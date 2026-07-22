"""VDS -> MatrixTable conversion must densify exactly ONCE.

Hail is lazy and does not cache. `to_dense_mt` builds a plan; every *eager* action on the
resulting MatrixTable (an aggregate, a count, a write) re-executes that plan from the top --
including the densify. `convert_vds_to_mt` used to issue two eager actions on the dense MT: an
`aggregate_entries` to validate biallelic entries, and then the `write`. So it densified the whole
cohort twice, and the validation pass cost ~42% of the stage's wall time on a 500-sample chr1
callset.

The fix moves the *audit* onto the sparse `variant_data` (where the defects can actually
originate) and applies the *repair* as a lazy, unconditional expression that folds into the single
write pass. These tests pin that: they count eager actions rather than measure time, so they run in
the fast suite with Hail stubbed out entirely.
"""

from unittest.mock import MagicMock, patch

import pytest

from hvantk.algorithms.hgc.converters import convert_vds_to_mt


class SpyMT:
    """A stand-in MatrixTable that records the calls made against it.

    Hail's MatrixTable methods are a mix of lazy (annotate_entries, key_cols_by -> return a new
    MT) and eager (aggregate_entries, write -> execute the plan). We return `self` from the lazy
    ones so a chain of them stays observable on one object, and record everything in order.
    """

    EAGER = {
        "aggregate_entries",
        "aggregate_rows",
        "aggregate_cols",
        "count",
        "count_rows",
        "write",
    }

    def __init__(
        self, name, entry_fields=("GT", "AD", "GQ", "DP"), aggregate_result=None
    ):
        self.name = name
        self.entry = dict.fromkeys(entry_fields)
        self.calls = []
        self.annotate_entries_kwargs = []
        self._aggregate_result = aggregate_result

    # -- lazy ops: record, return self so the chain stays on this spy --
    def annotate_entries(self, **kwargs):
        self.calls.append("annotate_entries")
        self.annotate_entries_kwargs.append(set(kwargs))
        return self

    def key_cols_by(self, *_a, **_k):
        self.calls.append("key_cols_by")
        return self

    # -- eager ops: these are the ones that force a densify --
    def aggregate_entries(self, *_a, **_k):
        self.calls.append("aggregate_entries")
        return self._aggregate_result

    def write(self, *_a, **_k):
        self.calls.append("write")

    def __getitem__(self, key):
        return MagicMock(name=f"{self.name}[{key}]")

    def __getattr__(self, item):
        # locus / alleles / GT / AD ... -> opaque expression objects
        return MagicMock(name=f"{self.name}.{item}")

    @property
    def n_eager(self):
        return sum(1 for c in self.calls if c in self.EAGER)


class SpyVDS:
    def __init__(self, variant_data, reference_data=None):
        self.variant_data = variant_data
        self.reference_data = reference_data or SpyMT("reference_data")


def _clean_audit():
    """What the audit aggregate returns when the data is fine (the expected case)."""
    return MagicMock(n_invalid_gt=0, gt_examples=[], n_invalid_ad=0, ad_examples=[])


@pytest.fixture
def spies():
    """Patch converters.hl so no Hail/Spark is needed; hand back the two spy MTs."""
    variant_data = SpyMT("variant_data", aggregate_result=_clean_audit())
    dense = SpyMT("dense_mt", aggregate_result=_clean_audit())
    vds = SpyVDS(variant_data)

    with patch("hvantk.algorithms.hgc.converters.hl") as hl:
        hl.vds.read_vds.return_value = vds
        hl.vds.split_multi.return_value = vds
        hl.vds.to_dense_mt.return_value = dense
        # hl.if_else / hl.missing / hl.any ... just need to be inert expression builders
        yield {"hl": hl, "vds": vds, "variant_data": variant_data, "dense": dense}


def _run(**overrides):
    kwargs = dict(
        vds_path="/in.vds",
        output_path="/out.mt",
        adjust_genotypes=False,  # keep gnomad out of this test
    )
    kwargs.update(overrides)
    convert_vds_to_mt(**kwargs)


def test_dense_mt_is_materialized_exactly_once(spies):
    """THE regression test: only the write may force the dense plan to execute.

    Fails on the old code, which called aggregate_entries on the dense MT *and* wrote it -> two
    eager actions -> two densifies.
    """
    _run()
    dense = spies["dense"]

    assert dense.calls.count("aggregate_entries") == 0, (
        "the biallelic audit must not run on the dense MatrixTable -- an eager aggregate there "
        f"forces a second full densify (calls seen: {dense.calls})"
    )
    assert dense.calls.count("write") == 1
    assert (
        dense.n_eager == 1
    ), f"expected exactly one densify-forcing action, got {dense.calls}"


def test_audit_runs_on_sparse_variant_data(spies):
    """The audit must move to variant_data -- that is the whole point of the fix."""
    _run()
    assert spies["variant_data"].calls.count("aggregate_entries") == 1


def test_gt_repair_is_applied_unconditionally_even_when_audit_is_clean(spies):
    """The repair must be a lazy expression, not gated on a Python-visible count.

    Gating it on `if n_invalid_gt > 0` is precisely what forced the second densify: reading that
    count is an eager action. So the repair has to be applied always -- on clean data it is the
    identity (`if_else(False, missing, GT) == GT`).
    """
    _run()
    dense = spies["dense"]

    gt_annotations = [kw for kw in dense.annotate_entries_kwargs if kw == {"GT"}]
    assert len(gt_annotations) == 1, (
        "expected exactly one lazy GT-repair annotate_entries on the dense MT, got "
        f"{dense.annotate_entries_kwargs}"
    )
    assert dense.calls.index("annotate_entries") < dense.calls.index("write")


def test_skip_validation_skips_both_audit_and_repair(spies):
    """skip_validation=True keeps its documented meaning: no audit, no repair."""
    _run(skip_validation=True)

    assert spies["variant_data"].calls.count("aggregate_entries") == 0
    assert {"GT"} not in spies["dense"].annotate_entries_kwargs
    assert spies["dense"].calls.count("write") == 1


def test_audit_tolerates_unsplit_variant_data(spies):
    """With skip_split_multi=True the entries are LGT/LAD, not GT/AD.

    The old validator dereferenced mt.GT / mt.AD unconditionally and blew up on this path. The
    audit must degrade to a warning instead of raising.
    """
    spies["variant_data"].entry = dict.fromkeys(("LGT", "LAD", "LA"))
    spies["dense"].entry = dict.fromkeys(("LGT", "LAD", "LA"))

    _run(skip_split_multi=True)  # must not raise

    assert spies["dense"].calls.count("write") == 1
