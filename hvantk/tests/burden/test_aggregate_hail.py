import pytest
import hail as hl

from hvantk.algorithms.burden.aggregate import (
    build_per_gene_carrier_mt,
    assert_clean_mt,
    qualifies_expr,
)

pytestmark = pytest.mark.hail


def _toy_mt():
    """3 variants (2 in GENEA/lof, 1 in GENEB/mis) x 4 samples.

    Carrier pattern: variant0 (GENEA/lof) het in s0 (case); variant1
    (GENEA/lof) het in s1 (case); variant2 (GENEB/mis) het in s2 (control).
    Everyone else is hom-ref. s0/s1 are cases; s2/s3 are controls.

    Built on ``hl.utils.range_matrix_table`` (not ``from_rows_table``, which
    yields a zero-column MT that cannot be given samples via annotate_cols).
    Shared by later tasks' tests -- keep it a clean, reusable helper.
    """
    mt = hl.utils.range_matrix_table(n_rows=3, n_cols=4)

    genes = hl.literal(["GENEA", "GENEA", "GENEB"])
    routes = hl.literal(["lof", "lof", "mis"])
    positions = hl.literal([100, 200, 300])
    alt_alleles = hl.literal(["C", "G", "T"])
    mt = mt.annotate_rows(
        locus=hl.locus("1", positions[mt.row_idx], reference_genome="GRCh38"),
        alleles=hl.array(["A", alt_alleles[mt.row_idx]]),
        SYMBOL=genes[mt.row_idx],
        csq_group=routes[mt.row_idx],
    )
    mt = mt.key_rows_by("locus", "alleles")

    sample_ids = hl.literal(["s0", "s1", "s2", "s3"])
    mt = mt.annotate_cols(s=sample_ids[mt.col_idx], is_case=mt.col_idx < 2)
    mt = mt.key_cols_by("s")

    # variant_i is carried het by sample_i for i in {0, 1, 2}; everyone else hom-ref.
    mt = mt.annotate_entries(
        GT=hl.if_else(mt.row_idx == mt.col_idx, hl.call(0, 1), hl.call(0, 0))
    )
    return mt.drop("row_idx", "col_idx")


def test_carrier_primitive_schema_no_route():
    mt = _toy_mt()
    out = build_per_gene_carrier_mt(mt, gene_field="SYMBOL")
    assert set(out.entry) == {"hets", "homs", "multi_het"}
    assert list(out.row_key) == ["SYMBOL"]


def test_carrier_primitive_with_route_keys_on_gene_and_route():
    mt = _toy_mt()
    out = build_per_gene_carrier_mt(mt, gene_field="SYMBOL", route_field="csq_group")
    assert list(out.row_key) == ["SYMBOL", "csq_group"]
    # GENEA/lof exists, GENEB/mis exists
    rows = out.rows()
    genes = rows.aggregate(hl.agg.collect_as_set(rows.SYMBOL))
    assert genes == {"GENEA", "GENEB"}


def test_assert_clean_mt_rejects_missing_route():
    mt = _toy_mt().drop("csq_group")
    with pytest.raises(ValueError, match="route column 'csq_group' not found"):
        assert_clean_mt(
            mt,
            gene_col="SYMBOL",
            route_col="csq_group",
            arm_col="is_case",
            key="symbol",
        )


def test_assert_clean_mt_rejects_multiallelic():
    # A tiny, standalone 1x1 MT whose single row is multiallelic (3 alleles),
    # with the gene/route/arm/GT fields assert_clean_mt requires.
    mt = hl.utils.range_matrix_table(n_rows=1, n_cols=1)
    mt = mt.annotate_rows(
        locus=hl.locus("1", 400, reference_genome="GRCh38"),
        alleles=["A", "C", "G"],
        SYMBOL="GENEC",
        csq_group="lof",
    )
    mt = mt.key_rows_by("locus", "alleles")
    mt = mt.annotate_cols(s="s0", is_case=True)
    mt = mt.key_cols_by("s")
    mt = mt.annotate_entries(GT=hl.call(0, 0))
    mt = mt.drop("row_idx", "col_idx")

    with pytest.raises(ValueError, match="multiallelic"):
        assert_clean_mt(
            mt,
            gene_col="SYMBOL",
            route_col="csq_group",
            arm_col="is_case",
            key="symbol",
        )


def test_assert_clean_mt_rejects_non_binary_arm():
    mt = _toy_mt().annotate_cols(is_case=hl.missing(hl.tbool))
    with pytest.raises(
        ValueError, match="case/control column 'is_case' must be a defined boolean"
    ):
        assert_clean_mt(
            mt,
            gene_col="SYMBOL",
            route_col="csq_group",
            arm_col="is_case",
            key="symbol",
        )
