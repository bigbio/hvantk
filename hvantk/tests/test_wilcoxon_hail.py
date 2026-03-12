"""Integration tests for Wilcoxon marker detection with Hail.

Uses the same synthetic MatrixTable pattern as test_summarize_expression.py.
Requires Hail; marked with @pytest.mark.hail.
"""

import shutil
import tempfile

import hail as hl
import numpy as np
import pytest

from hvantk.utils.wilcoxon import WilcoxonParams
from hvantk.utils.wilcoxon_hail import (
    extract_expression_for_wilcoxon,
    wilcoxon_markers_from_mt,
)
from hvantk.utils.gene_sets import GeneSetCollection

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
    2 regions (R1, R2).
    Gene G1/G2 are markers for type A; G4/G5 for type B; G3 is housekeeping.
    """
    rng = np.random.RandomState(42)
    n_genes, n_cells = 5, 20
    gene_ids = [f"ENSG{i:04d}" for i in range(n_genes)]
    gene_names = [f"G{i + 1}" for i in range(n_genes)]
    cell_ids = [f"cell_{i:03d}" for i in range(n_cells)]

    cell_types = ["A"] * 12 + ["B"] * 8
    regions = (["R1", "R2"] * 10)[:n_cells]
    n_umis = [int(x) for x in rng.randint(500, 5000, size=n_cells)]

    expr = np.zeros((n_genes, n_cells), dtype=float)
    for j in range(n_cells):
        if cell_types[j] == "A":
            expr[0, j] = rng.uniform(5, 10)
            expr[1, j] = rng.uniform(4, 8)
            expr[3, j] = rng.uniform(0, 1)
            expr[4, j] = rng.uniform(0, 1)
        else:
            expr[0, j] = rng.uniform(0, 1)
            expr[1, j] = rng.uniform(0, 1)
            expr[3, j] = rng.uniform(5, 10)
            expr[4, j] = rng.uniform(4, 8)
        expr[2, j] = rng.uniform(3, 5)

    flat_expr = [float(v) for v in expr.flatten()]

    mt = hl.utils.range_matrix_table(n_rows=n_genes, n_cols=n_cells)
    mt = mt.annotate_entries(x=hl.literal(flat_expr)[mt.row_idx * n_cells + mt.col_idx])
    mt = mt.annotate_rows(
        GeneID=hl.literal(gene_ids)[mt.row_idx],
    )
    mt = mt.annotate_rows(**{"Gene Name": hl.literal(gene_names)[mt.row_idx]})
    mt = mt.annotate_cols(
        sample_id=hl.literal(cell_ids)[mt.col_idx],
        metadata=hl.struct(
            cell_type=hl.literal(cell_types)[mt.col_idx],
            region=hl.literal(regions)[mt.col_idx],
            n_UMI=hl.literal(n_umis)[mt.col_idx],
        ),
    )
    mt = mt.key_rows_by("GeneID")
    mt = mt.key_cols_by("sample_id")
    mt = mt.drop("row_idx", "col_idx")

    return mt


class TestExtractExpressionForWilcoxon:
    def test_basic_extraction(self, synthetic_mt):
        expression, labels, gene_ids, gene_names = extract_expression_for_wilcoxon(
            synthetic_mt,
            group_by=["cell_type"],
            gene_id_field="GeneID",
            gene_name_field="Gene Name",
            min_cells_per_group=1,
        )
        assert expression.shape == (20, 5)
        assert len(labels) == 20
        assert set(labels) == {"A", "B"}
        assert len(gene_ids) == 5
        assert gene_names is not None
        assert len(gene_names) == 5

    def test_multi_group_by(self, synthetic_mt):
        expression, labels, gene_ids, _ = extract_expression_for_wilcoxon(
            synthetic_mt,
            group_by=["cell_type", "region"],
            gene_id_field="GeneID",
            gene_name_field="Gene Name",
            min_cells_per_group=1,
        )
        # Labels should be like "A_R1", "A_R2", "B_R1", "B_R2"
        assert any("_" in lbl for lbl in labels)
        unique_labels = set(labels)
        assert len(unique_labels) >= 3  # At least 3 combo groups

    def test_filter_by(self, synthetic_mt):
        expression, labels, _, _ = extract_expression_for_wilcoxon(
            synthetic_mt,
            group_by=["cell_type"],
            gene_id_field="GeneID",
            filter_by={"region": "R1"},
            min_cells_per_group=1,
        )
        assert expression.shape[0] < 20  # Fewer cells after filtering

    def test_candidate_filter(self, synthetic_mt):
        expression, labels, gene_ids, _ = extract_expression_for_wilcoxon(
            synthetic_mt,
            group_by=["cell_type"],
            gene_id_field="GeneID",
            candidate_gene_ids={"ENSG0000", "ENSG0001"},
            min_cells_per_group=1,
        )
        assert expression.shape[1] == 2
        assert set(gene_ids) == {"ENSG0000", "ENSG0001"}


class TestWilcoxonMarkersFromMt:
    def test_end_to_end(self, synthetic_mt):
        params = WilcoxonParams(
            min_fold_change=1.0,
            min_fraction_expressed=0.0,
            alpha=0.05,
            top_n=5,
        )
        results_df, collection = wilcoxon_markers_from_mt(
            synthetic_mt,
            group_by="cell_type",
            params=params,
            gene_id_field="GeneID",
            gene_name_field="Gene Name",
            min_cells_per_group=1,
        )

        assert len(results_df) > 0
        assert isinstance(collection, GeneSetCollection)
        assert "group" in results_df.columns
        assert "pvalue_adj" in results_df.columns

    def test_markers_make_sense(self, synthetic_mt):
        """G1/G2 should be markers for A; G4/G5 for B."""
        params = WilcoxonParams(
            min_fold_change=1.0,
            min_fraction_expressed=0.0,
            alpha=0.05,
            top_n=5,
        )
        results_df, collection = wilcoxon_markers_from_mt(
            synthetic_mt,
            group_by="cell_type",
            params=params,
            gene_id_field="GeneID",
            gene_name_field="Gene Name",
            min_cells_per_group=1,
        )

        a_set = collection.get("A")
        if a_set:
            assert "G1" in a_set.genes or "G2" in a_set.genes

        b_set = collection.get("B")
        if b_set:
            assert "G4" in b_set.genes or "G5" in b_set.genes

    def test_with_summary_prefilter(self, synthetic_mt, temp_dir):
        """Using a pre-computed summary should narrow candidates."""
        from pathlib import Path
        from hvantk.utils.matrix_utils import summarize_expression

        summary_path = str(Path(temp_dir) / "wilcoxon_summary.ht")
        summary_tb = summarize_expression(
            synthetic_mt,
            group_by="cell_type",
            gene_id_field="GeneID",
            gene_name_field="Gene Name",
            min_cells_per_group=1,
            output_path=summary_path,
            overwrite=True,
        )

        params = WilcoxonParams(
            min_fold_change=1.0,
            min_fraction_expressed=0.0,
            alpha=0.05,
        )
        results_df, collection = wilcoxon_markers_from_mt(
            synthetic_mt,
            group_by="cell_type",
            summary=summary_path,
            params=params,
            gene_id_field="GeneID",
            gene_name_field="Gene Name",
            min_cells_per_group=1,
        )
        assert len(results_df) > 0
        assert isinstance(collection, GeneSetCollection)

    def test_multi_group_by(self, synthetic_mt):
        params = WilcoxonParams(
            min_fold_change=1.0,
            min_fraction_expressed=0.0,
            alpha=0.5,  # lenient for small groups
        )
        results_df, collection = wilcoxon_markers_from_mt(
            synthetic_mt,
            group_by=["cell_type", "region"],
            params=params,
            gene_id_field="GeneID",
            gene_name_field="Gene Name",
            min_cells_per_group=1,
        )
        groups_found = set(results_df["group"].unique())
        assert len(groups_found) >= 3
