"""Regression: `adj` must not be MISSING on reference-derived (hom-ref) entries.

DeepVariant gVCF reference blocks are ``GT:GQ:MIN_DP:PL`` -- they carry ``MIN_DP`` and
no ``DP``. Hail's VDS combiner keeps what it finds, so ``reference_data`` ends up as
``{END, GQ, MIN_DP}``; ``to_dense_mt`` then fills variant-only fields with missing on
the reference side, leaving ``DP`` undefined on every hom-ref entry.

``adj = (GQ >= 20) & (DP >= 10) & (allele-balance)``. With ``DP`` missing, ``DP >= 10``
is missing and ``True & missing`` is missing -- so ``adj`` is MISSING rather than False.
``filter_entries(adj)`` keeps only True, so every reference-block hom-ref genotype is
DELETED, silently. On the 1005-sample CHD WGS cohort that destroyed 96.5% of hom-ref
entries (call rate 0.033) and went unnoticed for ~18 months, because ``variant_qc``
recomputed AC/AF/AN afterwards and the output looked self-consistent.

Note a schema check cannot catch this: ``DP`` IS in the entry schema (it comes from
``variant_data``) while being undefined on the reference side.
"""

import pytest

import hail as hl

from hvantk.algorithms.hgc.adj import annotate_adj
from hvantk.algorithms.hgc.converters import (
    _assert_adj_is_computable,
    resolve_reference_depth,
)


def _densified_shape_mt():
    """A MatrixTable shaped like a densified DeepVariant VDS.

    One row, four columns, each column a case that matters:

    ==== ================= ==== ====== ==== ============================
    col  origin            GQ   DP     MIN_DP  expected adj
    ==== ================= ==== ====== ==== ============================
    0    reference block   50   MISSING  25   True  -- the regression
    1    variant record    50   20       MISSING True
    2    reference block   50   MISSING  3    False -- depth too low
    3    reference block   5    MISSING  25   False -- GQ too low
    ==== ================= ==== ====== ==== ============================

    Columns 2 and 3 matter as much as 0: they prove the fallback still *filters*
    rather than waving everything through.
    """
    mt = hl.utils.range_matrix_table(n_rows=1, n_cols=4)
    i = mt.col_idx
    is_variant_record = i == 1
    return mt.annotate_entries(
        GT=hl.if_else(is_variant_record, hl.call(0, 1), hl.call(0, 0)),
        # AD is a variant-only field: missing on every reference-derived entry.
        AD=hl.if_else(
            is_variant_record,
            hl.literal([10, 10]),
            hl.missing(hl.tarray(hl.tint32)),
        ),
        GQ=hl.if_else(i == 3, 5, 50),
        DP=hl.if_else(is_variant_record, 20, hl.missing(hl.tint32)),
        MIN_DP=hl.if_else(
            is_variant_record,
            hl.missing(hl.tint32),
            hl.if_else(i == 2, 3, 25),
        ),
    )


def _adj_by_col(mt):
    entries = annotate_adj(mt).entries()
    return {r.col_idx: r.adj for r in entries.select("adj").collect()}


def test_adj_is_never_missing_on_reference_derived_entries():
    """The regression itself: no entry may come out MISSING."""
    adj = _adj_by_col(_densified_shape_mt())

    missing = [c for c, v in adj.items() if v is None]
    assert not missing, (
        f"adj is MISSING for columns {missing}. filter_entries() would DELETE those "
        f"genotypes silently. Has the MIN_DP fallback in annotate_adj been removed?"
    )


def test_adj_still_discriminates_after_the_fallback():
    """The fallback must not turn adj into a rubber stamp."""
    adj = _adj_by_col(_densified_shape_mt())

    assert adj[0] is True, "ref block with MIN_DP=25, GQ=50 should pass adj"
    assert adj[1] is True, "variant record with DP=20, GQ=50, balanced AD should pass"
    assert adj[2] is False, "ref block with MIN_DP=3 should FAIL the depth threshold"
    assert adj[3] is False, "ref block with GQ=5 should FAIL the GQ threshold"


def test_filter_entries_retains_reference_derived_genotypes():
    """End-to-end: the genotypes survive the filter that used to delete them."""
    mt = annotate_adj(_densified_shape_mt())
    kept = mt.filter_entries(mt.adj, keep=True)

    assert kept.aggregate_entries(hl.agg.count()) == 2, (
        "expected the two passing entries (a reference block and a variant record) "
        "to survive filter_entries"
    )


def test_without_min_dp_adj_goes_missing():
    """Pins the mechanism: drop MIN_DP and the reference entries go MISSING again.

    This is what the shipped pipeline did. It documents *why* the fallback exists, so
    the fix is not mistaken for an incidental refactor and quietly reverted.
    """
    adj = _adj_by_col(_densified_shape_mt().drop("MIN_DP"))

    assert adj[0] is None, "without MIN_DP the reference-derived entry must go MISSING"
    assert adj[1] is True, "variant records carry their own DP and are unaffected"


def test_guard_rejects_a_systematically_missing_adj():
    """The guard must catch a matrix already built with the broken adj."""
    mt = annotate_adj(_densified_shape_mt().drop("MIN_DP"))

    with pytest.raises(ValueError, match="MISSING"):
        _assert_adj_is_computable(mt)


def test_guard_passes_a_healthy_matrix():
    _assert_adj_is_computable(annotate_adj(_densified_shape_mt()))


def test_resolve_reference_depth_fills_dp_from_min_dp():
    """Recovered reference-block genotypes must carry a usable DP.

    Keeping the genotype but exporting it with an empty DP hands the next tool the
    same missing-value trap: a `FORMAT/DP >= 10` filter deletes exactly the entries
    the MIN_DP fallback just rescued. Folding MIN_DP into DP at the joint-genotyping
    step is the standard convention (GLnexus `orig_names: [MIN_DP, DP] -> DP,
    combi_method: min`; DRAGEN/GATK print hom-ref MIN_DP as FORMAT/DP).
    """
    mt = resolve_reference_depth(_densified_shape_mt())
    dp = {r.col_idx: r.DP for r in mt.entries().select("DP").collect()}

    assert dp[0] == 25, "reference-block entry must take DP from MIN_DP"
    assert dp[1] == 20, "variant-record DP must be preserved, not overwritten"
    assert dp[2] == 3, "low-depth reference block keeps its true (low) MIN_DP"
    assert dp[3] == 25

    assert all(v is not None for v in dp.values()), (
        "no entry may be left without a DP -- that is what makes a downstream "
        "FORMAT/DP filter delete recovered genotypes"
    )


def test_resolve_reference_depth_is_a_noop_without_min_dp():
    """Non-DeepVariant callsets have no MIN_DP; the step must leave them alone."""
    mt = _densified_shape_mt().drop("MIN_DP")
    out = resolve_reference_depth(mt)

    dp = {r.col_idx: r.DP for r in out.entries().select("DP").collect()}
    assert dp[1] == 20, "variant-record DP untouched"
    assert dp[0] is None, "nothing to fall back to, so DP stays missing"
