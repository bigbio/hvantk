"""The QC preflight aggregates were provably-zero counts, and must not come back.

`_prepare_qc_gt` writes a SANITIZED call into `__qc_gt`:

    invalid_gt = hl.is_defined(base_gt) & hl.any(lambda i: base_gt[i] >= hl.len(mt.alleles), ...)
    qc_gt      = hl.if_else(invalid_gt, hl.missing(hl.tcall), base_gt)

so `__qc_gt` is MISSING exactly where the genotype was invalid. compute_sample_qc and
compute_variant_qc then each ran a full `aggregate_entries` counting entries where `__qc_gt`
IS DEFINED *and* is invalid. That intersection is empty by construction: the count is the constant
0 and the `if invalid_count:` warning branch is unreachable.

Each of those aggregates was a full pass over the entry matrix (a 12 GiB MatrixTable on the
benchmark cohort) to compute a number that could not be anything but zero. These tests assert the
QC functions perform no eager entry aggregation at all.
"""

from unittest.mock import MagicMock, patch

import pytest

from hvantk.algorithms.hgc.qc import compute_sample_qc, compute_variant_qc


class SpyMT:
    """MatrixTable stand-in that records eager aggregations. Lazy ops return self."""

    def __init__(self):
        self.entry = dict.fromkeys(("GT", "AD", "GQ", "DP"))
        self.row = dict.fromkeys(("locus", "alleles"))
        self.calls = []

    def annotate_entries(self, **_kwargs):
        self.calls.append("annotate_entries")
        return self

    def drop(self, *_a):
        self.calls.append("drop")
        return self

    def aggregate_entries(self, *_a, **_k):
        self.calls.append("aggregate_entries")
        return 0

    def __getitem__(self, key):
        return MagicMock(name=f"mt[{key}]")

    def __getattr__(self, item):
        return MagicMock(name=f"mt.{item}")


@pytest.fixture
def spy_mt():
    mt = SpyMT()
    with patch("hvantk.algorithms.hgc.qc.hl") as hl:
        hl.sample_qc.side_effect = lambda m, **_k: m
        hl.variant_qc.side_effect = lambda m, **_k: m
        yield mt


@pytest.mark.parametrize("qc_fn", [compute_sample_qc, compute_variant_qc])
def test_qc_does_not_aggregate_entries(spy_mt, qc_fn):
    """No full entry scan may be issued just to compute a number that is always 0."""
    qc_fn(spy_mt)

    assert spy_mt.calls.count("aggregate_entries") == 0, (
        f"{qc_fn.__name__} issued a full entry aggregation; the preflight count it replaced was "
        f"tautologically zero (calls seen: {spy_mt.calls})"
    )
