"""
Tests for overlap enrichment analysis.
"""

import pytest

from hvantk.enrichex.gene_sets import (
    GeneSet,
    GeneSetCollection,
    load_gene_sets_from_dict,
)
from hvantk.enrichex.overlap import (
    OverlapResult,
    compute_overlap_enrichment,
    compute_overlap_enrichment_pandas,
)


@pytest.mark.hail
class TestOverlapResult:
    """Tests for OverlapResult dataclass."""

    def test_overlap_result_creation(self):
        """Test basic OverlapResult creation."""
        result = OverlapResult(
            gene_set_name="test_set",
            n_query=10,
            n_gene_set=20,
            n_overlap=5,
            n_background=100,
            p_value=0.01,
            odds_ratio=2.5,
            ci_lower=1.2,
            ci_upper=4.0,
            overlap_genes=["BRCA1", "TP53"],
            p_adjusted=0.05,
        )

        assert result.gene_set_name == "test_set"
        assert result.n_overlap == 5
        assert result.p_value == 0.01
        assert len(result.overlap_genes) == 2

    def test_overlap_result_to_dict(self):
        """Test OverlapResult serialization."""
        result = OverlapResult(
            gene_set_name="test",
            n_query=10,
            n_gene_set=10,
            n_overlap=5,
            n_background=100,
            p_value=0.01,
            odds_ratio=2.0,
            ci_lower=1.0,
            ci_upper=3.0,
            overlap_genes=["A", "B"],
        )

        d = result.to_dict()

        assert d["gene_set_name"] == "test"
        assert d["n_overlap"] == 5
        assert d["p_value"] == 0.01
        assert d["overlap_genes"] == ["A", "B"]


@pytest.mark.hail
class TestComputeOverlapEnrichment:
    """Tests for compute_overlap_enrichment function."""

    def test_simple_enrichment(self, hail_session):
        """Test basic enrichment calculation."""
        # Create gene sets
        gene_sets_dict = {
            "enriched": ["A", "B", "C", "D", "E", "F"],  # 6 genes
            "not_enriched": ["X", "Y", "Z"],  # 3 genes
        }

        # Background of 26 genes (A-Z)
        background = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
        collection = load_gene_sets_from_dict(
            gene_sets_dict, background_genes=background
        )

        # Query with 5 genes, all in "enriched" set
        query = ["A", "B", "C", "D", "E"]

        results = compute_overlap_enrichment(
            query, collection, correction_method="none"
        )

        # Should have 2 results
        assert len(results) == 2

        # Find enriched set result
        enriched = next(r for r in results if r.gene_set_name == "enriched")
        not_enriched = next(r for r in results if r.gene_set_name == "not_enriched")

        # Enriched should have better p-value
        assert enriched.p_value < not_enriched.p_value

        # Enriched should have high overlap
        assert enriched.n_overlap == 5
        assert not_enriched.n_overlap == 0

        # Odds ratio for enriched should be > 1
        assert enriched.odds_ratio > 1

    def test_no_enrichment(self, hail_session):
        """Test case with no enrichment."""
        gene_sets_dict = {
            "set1": ["A", "B", "C"],
            "set2": ["D", "E", "F"],
        }

        background = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
        collection = load_gene_sets_from_dict(
            gene_sets_dict, background_genes=background
        )

        # Query genes not in any set
        query = ["X", "Y", "Z"]

        results = compute_overlap_enrichment(
            query, collection, correction_method="none"
        )

        # Should have results but no overlap
        assert len(results) == 2
        for r in results:
            assert r.n_overlap == 0

    def test_multiple_testing_correction(self, hail_session):
        """Test that multiple testing correction is applied."""
        gene_sets_dict = {
            "set1": ["A", "B", "C"],
            "set2": ["D", "E", "F"],
            "set3": ["G", "H", "I"],
        }

        background = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
        collection = load_gene_sets_from_dict(
            gene_sets_dict, background_genes=background
        )

        query = ["A", "B", "C"]

        # Test with Benjamini-Hochberg
        results_bh = compute_overlap_enrichment(
            query, collection, correction_method="benjamini-hochberg"
        )

        # All results should have p_adjusted
        assert all(r.p_adjusted is not None for r in results_bh)

        # Adjusted p-values should be >= raw p-values
        for r in results_bh:
            assert r.p_adjusted >= r.p_value

        # Test with Bonferroni
        results_bonf = compute_overlap_enrichment(
            query, collection, correction_method="bonferroni"
        )

        # Bonferroni should be more conservative
        for r_bh, r_bonf in zip(results_bh, results_bonf):
            assert (
                r_bonf.p_adjusted >= r_bh.p_adjusted
                or abs(r_bonf.p_adjusted - r_bh.p_adjusted) < 1e-10
            )

    def test_results_sorted_by_pvalue(self, hail_session):
        """Test that results are sorted by p-value."""
        gene_sets_dict = {
            "high_overlap": ["A", "B", "C", "D", "E"],  # 5 genes
            "low_overlap": ["A", "X", "Y", "Z"],  # 1 overlap
            "no_overlap": ["X", "Y", "Z"],  # 0 overlap
        }

        background = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
        collection = load_gene_sets_from_dict(
            gene_sets_dict, background_genes=background
        )

        query = ["A", "B", "C", "D", "E"]

        results = compute_overlap_enrichment(
            query, collection, correction_method="none"
        )

        # Results should be sorted by p-value (ascending)
        p_values = [r.p_value for r in results]
        assert p_values == sorted(p_values)

        # First result should be high_overlap
        assert results[0].gene_set_name == "high_overlap"

    def test_background_filtering(self, hail_session):
        """Test that genes not in background are filtered out."""
        gene_sets_dict = {
            "set1": ["A", "B", "C"],
        }

        # Small background
        background = {"A", "B", "C", "D", "E"}
        collection = load_gene_sets_from_dict(
            gene_sets_dict, background_genes=background
        )

        # Query includes gene not in background
        query = ["A", "B", "Z"]  # Z not in background

        results = compute_overlap_enrichment(
            query, collection, correction_method="none"
        )

        # Only A and B should be counted (2 genes, not 3)
        assert results[0].n_query == 2

    def test_empty_query(self, hail_session):
        """Test with empty query list."""
        gene_sets_dict = {"set1": ["A", "B", "C"]}
        background = set("ABCDEFGHIJ")
        collection = load_gene_sets_from_dict(
            gene_sets_dict, background_genes=background
        )

        query = []

        results = compute_overlap_enrichment(
            query, collection, correction_method="none"
        )

        # Should return empty results
        assert len(results) == 0

    def test_query_not_in_background(self, hail_session):
        """Test when no query genes are in background."""
        gene_sets_dict = {"set1": ["A", "B", "C"]}
        background = {"A", "B", "C", "D", "E"}
        collection = load_gene_sets_from_dict(
            gene_sets_dict, background_genes=background
        )

        # Query genes not in background
        query = ["X", "Y", "Z"]

        results = compute_overlap_enrichment(
            query, collection, correction_method="none"
        )

        # Should return empty (warning logged)
        assert len(results) == 0

    def test_overlap_genes_list(self, hail_session):
        """Test that overlap_genes contains correct genes."""
        gene_sets_dict = {
            "set1": ["A", "B", "C", "D"],
        }

        background = set("ABCDEFGHIJ")
        collection = load_gene_sets_from_dict(
            gene_sets_dict, background_genes=background
        )

        query = ["A", "B", "X", "Y"]

        results = compute_overlap_enrichment(
            query, collection, correction_method="none"
        )

        # Overlap should be A and B
        assert set(results[0].overlap_genes) == {"A", "B"}
        # Should be sorted
        assert results[0].overlap_genes == sorted(results[0].overlap_genes)


@pytest.mark.hail
class TestComputeOverlapEnrichmentPandas:
    """Tests for compute_overlap_enrichment_pandas function."""

    def test_pandas_output_format(self, hail_session):
        """Test that pandas function returns DataFrame."""
        gene_sets_dict = {
            "set1": ["A", "B", "C"],
            "set2": ["D", "E", "F"],
        }

        background = set("ABCDEFGHIJ")
        collection = load_gene_sets_from_dict(
            gene_sets_dict, background_genes=background
        )

        query = ["A", "B"]

        df = compute_overlap_enrichment_pandas(
            query, collection, correction_method="none"
        )

        # Check it's a DataFrame
        assert hasattr(df, "columns")

        # Check expected columns
        expected_columns = [
            "gene_set_name",
            "n_query",
            "n_gene_set",
            "n_overlap",
            "n_background",
            "p_value",
            "odds_ratio",
            "ci_lower",
            "ci_upper",
            "p_adjusted",
            "overlap_genes",
            "significant",
        ]

        for col in expected_columns:
            assert col in df.columns

    def test_pandas_significant_column(self, hail_session):
        """Test that significant column is correctly computed."""
        gene_sets_dict = {
            "enriched": ["A", "B", "C", "D", "E"],
            "not_enriched": ["X", "Y", "Z"],
        }

        background = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
        collection = load_gene_sets_from_dict(
            gene_sets_dict, background_genes=background
        )

        query = ["A", "B", "C", "D", "E"]

        df = compute_overlap_enrichment_pandas(
            query, collection, correction_method="benjamini-hochberg"
        )

        # Check significant column exists and is boolean
        assert "significant" in df.columns
        assert df["significant"].dtype == bool

        # Enriched set should be significant
        enriched_row = df[df["gene_set_name"] == "enriched"].iloc[0]
        assert enriched_row["significant"] == True

    def test_pandas_overlap_genes_string(self, hail_session):
        """Test that overlap_genes is comma-separated string in DataFrame."""
        gene_sets_dict = {"set1": ["A", "B", "C"]}
        background = set("ABCDEFGHIJ")
        collection = load_gene_sets_from_dict(
            gene_sets_dict, background_genes=background
        )

        query = ["A", "B"]

        df = compute_overlap_enrichment_pandas(
            query, collection, correction_method="none"
        )

        # overlap_genes should be string
        overlap_str = df.iloc[0]["overlap_genes"]
        assert isinstance(overlap_str, str)
        assert "," in overlap_str or len(overlap_str.split(",")) == 1

        # Should be able to split back to list
        genes = overlap_str.split(",")
        assert set(genes) == {"A", "B"}

    def test_pandas_empty_results(self, hail_session):
        """Test DataFrame structure with empty results."""
        gene_sets_dict = {"set1": ["A", "B", "C"]}
        background = {"A", "B", "C"}
        collection = load_gene_sets_from_dict(
            gene_sets_dict, background_genes=background
        )

        # Query not in background
        query = ["X", "Y", "Z"]

        df = compute_overlap_enrichment_pandas(
            query, collection, correction_method="none"
        )

        # Should return empty DataFrame with correct columns
        assert len(df) == 0
        assert "gene_set_name" in df.columns
        assert "p_value" in df.columns


@pytest.mark.hail
class TestOverlapEnrichmentIntegration:
    """Integration tests for overlap enrichment."""

    def test_realistic_scenario(self, hail_session):
        """Test realistic gene set enrichment scenario."""
        # Simulate cell-type marker genes
        gene_sets_dict = {
            "T_cell": ["CD3D", "CD3E", "CD3G", "CD8A", "CD4"],
            "B_cell": ["CD19", "CD79A", "CD79B", "MS4A1"],
            "Macrophage": ["CD14", "CD68", "FCGR3A", "CSF1R"],
        }

        # Realistic background: ~50 genes
        background = (
            set(gene_sets_dict["T_cell"])
            | set(gene_sets_dict["B_cell"])
            | set(gene_sets_dict["Macrophage"])
        )
        background.update([f"GENE{i}" for i in range(40)])  # Add more genes

        collection = load_gene_sets_from_dict(
            gene_sets_dict, background_genes=background
        )

        # Query: T cell specific genes
        query = ["CD3D", "CD3E", "CD8A", "GENE1", "GENE2"]

        results = compute_overlap_enrichment(
            query, collection, correction_method="benjamini-hochberg"
        )

        # Should have 3 results
        assert len(results) == 3

        # T_cell should be most enriched
        assert results[0].gene_set_name == "T_cell"
        assert results[0].n_overlap == 3  # CD3D, CD3E, CD8A
        assert results[0].p_adjusted < 0.05

    def test_with_test_data_file(self, hail_session, tmp_path):
        """Test with actual test data files."""
        from pathlib import Path

        # Create simple test gene sets
        gene_sets = {
            "cancer_genes": ["BRCA1", "TP53", "EGFR", "KRAS", "MYC", "PTEN"],
            "cell_cycle": ["CDK1", "CDK2", "CCNA1", "CCNB1"],
        }

        background = set(gene_sets["cancer_genes"]) | set(gene_sets["cell_cycle"])
        background.update([f"GENE{i}" for i in range(50)])

        collection = load_gene_sets_from_dict(gene_sets, background_genes=background)

        # Save collection
        json_path = tmp_path / "gene_sets.json"
        collection.save(json_path)

        # Load back
        loaded_collection = GeneSetCollection.load(json_path)

        # Use for enrichment
        query = ["BRCA1", "TP53", "EGFR", "KRAS"]

        df = compute_overlap_enrichment_pandas(
            query, loaded_collection, correction_method="benjamini-hochberg"
        )

        # Should find cancer_genes enriched
        cancer_row = df[df["gene_set_name"] == "cancer_genes"].iloc[0]
        assert cancer_row["n_overlap"] == 4
        assert cancer_row["significant"] == True
