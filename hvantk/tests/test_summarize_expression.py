"""Tests for expression summarization pipeline (Phases 0–2).

Covers:
- annotate_column_summary / describe_expression_mt (Phase 0)
- summarize_expression (Phase 1)
- extract_marker_gene_sets (Phase 2)

Uses a small synthetic MatrixTable to keep tests fast.
"""

import shutil
import tempfile
from pathlib import Path

import hail as hl
import pandas as pd
import pytest

from hvantk.utils.matrix_utils import (
    annotate_column_summary,
    describe_expression_mt,
    summarize_expression,
)
from hvantk.utils.gene_sets import extract_marker_gene_sets, GeneSetCollection

pytestmark = [pytest.mark.hail]


@pytest.fixture(scope="module")
def temp_dir():
    d = tempfile.mkdtemp()
    yield d
    shutil.rmtree(d)


@pytest.fixture(scope="module")
def synthetic_mt():
    """Build a small synthetic expression MatrixTable.

    5 genes x 20 cells, 2 cell types (A, B with 12 and 8 cells),
    2 regions (R1, R2), one numeric field (n_UMI).
    Gene G1/G2 are markers for type A; G4/G5 for type B; G3 is housekeeping.

    Uses hl.utils.range_matrix_table to avoid numpy compatibility issues
    with hl.Table.from_pandas.
    """
    import numpy as np

    np.random.seed(42)
    n_genes, n_cells = 5, 20
    gene_ids = [f"ENSG{i:04d}" for i in range(n_genes)]
    gene_names = [f"G{i + 1}" for i in range(n_genes)]
    cell_ids = [f"cell_{i:03d}" for i in range(n_cells)]

    # Cell metadata
    cell_types = ["A"] * 12 + ["B"] * 8
    regions = (["R1", "R2"] * 10)[:n_cells]
    n_umis = [int(x) for x in np.random.randint(500, 5000, size=n_cells)]

    # Expression: markers skewed by cell type
    expr = np.zeros((n_genes, n_cells), dtype=float)
    for j in range(n_cells):
        if cell_types[j] == "A":
            expr[0, j] = np.random.uniform(5, 10)  # G1 high in A
            expr[1, j] = np.random.uniform(4, 8)  # G2 high in A
            expr[3, j] = np.random.uniform(0, 1)  # G4 low in A
            expr[4, j] = np.random.uniform(0, 1)  # G5 low in A
        else:
            expr[0, j] = np.random.uniform(0, 1)  # G1 low in B
            expr[1, j] = np.random.uniform(0, 1)  # G2 low in B
            expr[3, j] = np.random.uniform(5, 10)  # G4 high in B
            expr[4, j] = np.random.uniform(4, 8)  # G5 high in B
        expr[2, j] = np.random.uniform(3, 5)  # G3 housekeeping

    # Flatten expression for entry annotation (row-major)
    flat_expr = [float(v) for v in expr.flatten()]

    # Build MT from range_matrix_table (avoids np.bool issue in from_pandas)
    mt = hl.utils.range_matrix_table(n_rows=n_genes, n_cols=n_cells)

    # Annotate entries first (row_idx / col_idx still available)
    mt = mt.annotate_entries(
        x=hl.literal(flat_expr)[mt.row_idx * n_cells + mt.col_idx]
    )

    # Annotate rows
    mt = mt.annotate_rows(
        GeneID=hl.literal(gene_ids)[mt.row_idx],
    )
    mt = mt.annotate_rows(
        **{"Gene Name": hl.literal(gene_names)[mt.row_idx]}
    )

    # Annotate cols
    mt = mt.annotate_cols(
        sample_id=hl.literal(cell_ids)[mt.col_idx],
        metadata=hl.struct(
            cell_type=hl.literal(cell_types)[mt.col_idx],
            region=hl.literal(regions)[mt.col_idx],
            n_UMI=hl.literal(n_umis)[mt.col_idx],
        ),
    )

    # Re-key to match real expression MT layout
    mt = mt.key_rows_by("GeneID")
    mt = mt.key_cols_by("sample_id")
    mt = mt.drop("row_idx", "col_idx")

    return mt


# ── Phase 0 tests ──────────────────────────────────────────────────────────


class TestAnnotateColumnSummary:
    def test_basic(self, synthetic_mt):
        mt = annotate_column_summary(synthetic_mt)
        assert "column_summary" in mt.globals.dtype

        summary = hl.eval(mt.column_summary)
        assert "cell_type" in summary
        assert "region" in summary
        assert "n_UMI" in summary

    def test_categorical_fields(self, synthetic_mt):
        mt = annotate_column_summary(synthetic_mt)
        summary = hl.eval(mt.column_summary)

        ct = summary["cell_type"]
        assert ct.dtype == "categorical"
        assert ct.n_levels == 2
        assert set(ct.levels) == {"A", "B"}
        assert ct.truncated is False

    def test_numeric_fields(self, synthetic_mt):
        mt = annotate_column_summary(synthetic_mt)
        summary = hl.eval(mt.column_summary)

        n = summary["n_UMI"]
        assert n.dtype == "numeric"
        assert n.min_val >= 500
        assert n.max_val <= 5000
        assert n.mean_val > 0

    def test_high_cardinality_truncation(self, synthetic_mt):
        """With max_levels=1, all categorical fields should be truncated."""
        mt = annotate_column_summary(
            synthetic_mt, max_levels=1, top_n_levels=1
        )
        summary = hl.eval(mt.column_summary)
        ct = summary["cell_type"]
        assert ct.truncated is True
        assert len(ct.top_levels) == 1
        assert len(ct.levels) == 0


class TestDescribeExpressionMt:
    def test_with_summary(self, synthetic_mt):
        mt = annotate_column_summary(synthetic_mt)
        result = describe_expression_mt(mt)

        assert result["n_genes"] == 5
        assert result["n_cols"] == 20
        assert len(result["fields"]) == 3  # cell_type, region, n_UMI

    def test_fallback_without_summary(self, synthetic_mt):
        """Should compute summary on the fly when globals are missing."""
        result = describe_expression_mt(synthetic_mt)
        assert result["n_genes"] == 5
        assert len(result["fields"]) > 0


# ── Phase 1 tests ──────────────────────────────────────────────────────────


class TestSummarizeExpression:
    def test_single_group_by(self, synthetic_mt):
        tb = summarize_expression(
            synthetic_mt,
            group_by="cell_type",
            gene_id_field="GeneID",
            gene_name_field="Gene Name",
            min_cells_per_group=1,
        )
        assert tb.count() == 5
        assert "gene_id" in list(tb.key)
        assert "stats" in tb.row

        # Check stats structure
        row = tb.take(1)[0]
        assert "A" in row.stats
        assert "B" in row.stats
        assert hasattr(row.stats["A"], "mean")
        assert hasattr(row.stats["A"], "fraction_expressed")
        assert hasattr(row.stats["A"], "n_cells")

    def test_multi_group_by(self, synthetic_mt):
        tb = summarize_expression(
            synthetic_mt,
            group_by=["cell_type", "region"],
            gene_id_field="GeneID",
            gene_name_field="Gene Name",
            min_cells_per_group=1,
        )
        row = tb.take(1)[0]
        # Expect concatenated labels like "A_R1", "A_R2", "B_R1", "B_R2"
        assert any("_" in k for k in row.stats.keys())

    def test_filter_by(self, synthetic_mt):
        tb = summarize_expression(
            synthetic_mt,
            group_by="cell_type",
            filter_by={"region": "R1"},
            gene_id_field="GeneID",
            min_cells_per_group=1,
        )
        # Should still have groups A and B, but with fewer cells
        row = tb.take(1)[0]
        # n_cells should be less than full (12 A + 8 B)
        total = sum(s.n_cells for s in row.stats.values())
        assert total < 20

    def test_min_cells_filters_small_groups(self, synthetic_mt):
        # B has only 8 cells; setting min=10 should skip it
        tb = summarize_expression(
            synthetic_mt,
            group_by="cell_type",
            gene_id_field="GeneID",
            min_cells_per_group=10,
        )
        row = tb.take(1)[0]
        assert "A" in row.stats
        assert "B" not in row.stats

    def test_numeric_group_by_raises(self, synthetic_mt):
        """Grouping by a numeric field should raise ValueError."""
        mt = annotate_column_summary(synthetic_mt)
        with pytest.raises(ValueError, match="numeric"):
            summarize_expression(
                mt,
                group_by="n_UMI",
                gene_id_field="GeneID",
                min_cells_per_group=1,
            )

    def test_checkpoint(self, synthetic_mt, temp_dir):
        out = str(Path(temp_dir) / "summary.ht")
        tb = summarize_expression(
            synthetic_mt,
            group_by="cell_type",
            gene_id_field="GeneID",
            min_cells_per_group=1,
            output_path=out,
            overwrite=True,
        )
        # Should be readable from disk
        reloaded = hl.read_table(out)
        assert reloaded.count() == tb.count()


# ── Phase 2 tests ──────────────────────────────────────────────────────────


class TestExtractMarkerGeneSets:
    @pytest.fixture(scope="class")
    def summary_tb(self, synthetic_mt):
        return summarize_expression(
            synthetic_mt,
            group_by="cell_type",
            gene_id_field="GeneID",
            gene_name_field="Gene Name",
            min_cells_per_group=1,
        )

    def test_fold_change(self, summary_tb):
        collection = extract_marker_gene_sets(
            summary_tb,
            n_markers=5,
            min_fold_change=1.0,
            min_fraction_expressed=0.0,
            method="fold_change",
        )
        assert isinstance(collection, GeneSetCollection)
        assert len(collection) >= 1
        assert len(collection.background_genes) == 5

    def test_specificity(self, summary_tb):
        collection = extract_marker_gene_sets(
            summary_tb,
            n_markers=5,
            min_fold_change=1.0,
            min_fraction_expressed=0.0,
            method="specificity",
        )
        assert len(collection) >= 1

    def test_marker_genes_make_sense(self, summary_tb):
        """G1/G2 should be top markers for A; G4/G5 for B."""
        collection = extract_marker_gene_sets(
            summary_tb,
            n_markers=3,
            min_fold_change=1.0,
            min_fraction_expressed=0.0,
            method="fold_change",
        )
        # Check group A markers
        a_set = collection.get("A")
        if a_set:
            assert "G1" in a_set.genes or "G2" in a_set.genes
        # Check group B markers
        b_set = collection.get("B")
        if b_set:
            assert "G4" in b_set.genes or "G5" in b_set.genes

    def test_min_fraction_expressed_filter(self, summary_tb):
        # With very high min_fraction, fewer markers should pass
        strict = extract_marker_gene_sets(
            summary_tb,
            n_markers=5,
            min_fold_change=1.0,
            min_fraction_expressed=0.99,
        )
        lenient = extract_marker_gene_sets(
            summary_tb,
            n_markers=5,
            min_fold_change=1.0,
            min_fraction_expressed=0.0,
        )
        strict_total = sum(gs.n_genes for gs in strict)
        lenient_total = sum(gs.n_genes for gs in lenient)
        assert strict_total <= lenient_total

    def test_save_load_json(self, summary_tb, temp_dir):
        collection = extract_marker_gene_sets(
            summary_tb,
            n_markers=3,
            min_fold_change=1.0,
            min_fraction_expressed=0.0,
        )
        path = Path(temp_dir) / "markers.json"
        collection.save(path)
        reloaded = GeneSetCollection.load(path)
        assert len(reloaded) == len(collection)

    def test_save_gmt(self, summary_tb, temp_dir):
        collection = extract_marker_gene_sets(
            summary_tb,
            n_markers=3,
            min_fold_change=1.0,
            min_fraction_expressed=0.0,
        )
        path = Path(temp_dir) / "markers.gmt"
        collection.save_gmt(path)
        reloaded = GeneSetCollection.load_gmt(path)
        assert len(reloaded) == len(collection)
