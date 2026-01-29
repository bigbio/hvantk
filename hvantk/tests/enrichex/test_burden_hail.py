"""
Tests for burden analysis module (Hail-native regression).
"""

import pytest
import hail as hl
import pandas as pd

from hvantk.enrichex.burden import (
    VariantFilter,
    compute_geneset_burden_mt,
    logistic_burden_test,
    linear_burden_test,
    run_burden_analysis,
)


@pytest.mark.hail
class TestVariantFilter:
    """Tests for VariantFilter dataclass."""

    def test_variant_filter_creation(self):
        """Test basic VariantFilter creation."""
        vf = VariantFilter(
            max_af=0.01,
            min_cadd=20.0,
            consequences=["missense_variant", "frameshift_variant"],
        )

        assert vf.max_af == 0.01
        assert vf.min_cadd == 20.0
        assert len(vf.consequences) == 2

    def test_variant_filter_to_hail_expr(self, hail_session):
        """Test conversion to Hail expression."""
        # Create simple MT
        mt = hl.utils.range_matrix_table(10, 5)
        mt = mt.annotate_rows(
            gnomad_af=hl.rand_unif(0, 0.1),
            cadd_phred=hl.rand_unif(0, 40),
            consequence="missense_variant",
        )
        mt = mt.annotate_entries(GQ=30, DP=20)

        vf = VariantFilter(max_af=0.01, min_cadd=20.0)
        expr = vf.to_hail_expr(mt)

        # Should be a boolean expression
        assert isinstance(expr, hl.expr.BooleanExpression)

    def test_variant_filter_missing_fields(self, hail_session):
        """Test filter with missing fields."""
        mt = hl.utils.range_matrix_table(10, 5)

        vf = VariantFilter(max_af=0.01, min_cadd=20.0)

        # Should log warnings but not fail
        expr = vf.to_hail_expr(mt)
        assert expr is not None


@pytest.mark.hail
class TestComputeGenesetBurdenMt:
    """Tests for compute_geneset_burden_mt function."""

    def setup_test_mt(self):
        """Create a test MatrixTable with variants and genes."""
        # Create MT with 20 variants, 10 samples
        mt = hl.utils.range_matrix_table(20, 10)

        # Add gene annotations
        mt = mt.annotate_rows(
            SYMBOL=hl.if_else(
                mt.row_idx < 10, "GENE1", hl.if_else(mt.row_idx < 15, "GENE2", "GENE3")
            ),
            gnomad_af=hl.rand_unif(0, 0.05),
            cadd_phred=hl.rand_unif(15, 30),
        )

        # Add genotypes (random hets/homs)
        mt = mt.annotate_entries(
            GT=hl.call(
                hl.int(hl.rand_bool(0.3)), hl.int(hl.rand_bool(0.1))
            ),  # Random genotypes
            GQ=30,
            DP=20,
        )

        return mt

    def test_compute_geneset_burden_basic(self, hail_session):
        """Test basic burden computation."""
        mt = self.setup_test_mt()

        gene_sets = {"set1": ["GENE1", "GENE2"], "set2": ["GENE3"]}

        mt_burden = compute_geneset_burden_mt(
            mt, gene_sets, gene_field="SYMBOL", max_af=0.1, min_cadd=None
        )

        # Should have 2 gene sets (rows) and 10 samples (cols)
        n_rows, n_cols = mt_burden.count()
        assert n_rows == 2  # 2 gene sets
        assert n_cols == 10  # 10 samples

        # Should have burden field
        assert "burden" in mt_burden.entry

    def test_compute_geneset_burden_genotype_aggregation_methods(self, hail_session):
        """Test different genotype aggregation methods."""
        mt = self.setup_test_mt()
        gene_sets = {"set1": ["GENE1", "GENE2"]}

        for method in ["hets", "homs", "chets", "homs_chets"]:
            mt_burden = compute_geneset_burden_mt(
                mt,
                gene_sets,
                gene_field="SYMBOL",
                genotype_aggregation=method,
                max_af=0.1,
                min_cadd=None,
            )

            assert mt_burden.count_rows() == 1

    def test_compute_geneset_burden_no_variants(self, hail_session):
        """Test with very restrictive filters (no variants pass)."""
        mt = self.setup_test_mt()
        gene_sets = {"set1": ["GENE1"]}

        # Very restrictive AF filter
        with pytest.raises(ValueError, match="No variants in gene sets"):
            compute_geneset_burden_mt(
                mt,
                gene_sets,
                gene_field="SYMBOL",
                max_af=0.0001,  # Very low
                min_cadd=None,
            )

    def test_compute_geneset_burden_invalid_method(self, hail_session):
        """Test with invalid genotype aggregation method."""
        mt = self.setup_test_mt()
        gene_sets = {"set1": ["GENE1"]}

        with pytest.raises(ValueError, match="Invalid genotype_aggregation"):
            compute_geneset_burden_mt(
                mt,
                gene_sets,
                gene_field="SYMBOL",
                genotype_aggregation="invalid",
            )


@pytest.mark.hail
class TestLogisticBurdenTest:
    """Tests for logistic_burden_test function."""

    def setup_burden_mt(self):
        """Create a burden MatrixTable for testing."""
        # Create MT with 3 gene sets, 100 samples
        mt = hl.utils.range_matrix_table(3, 100)
        mt = mt.annotate_rows(gene_set_name=hl.str(mt.row_idx))

        # Add burden scores (random 0-5)
        mt = mt.annotate_entries(
            burden=hl.int(hl.rand_unif(0, 5)),
        )

        # Add phenotypes (50 cases, 50 controls)
        mt = mt.annotate_cols(
            is_case=mt.col_idx < 50,
            PC1=hl.rand_norm(),
            PC2=hl.rand_norm(),
            sex=mt.col_idx % 2,
        )

        return mt

    def test_logistic_burden_test_basic(self, hail_session):
        """Test basic logistic regression."""
        mt_burden = self.setup_burden_mt()

        result_ht = logistic_burden_test(mt_burden, phenotype_field="is_case")

        # Should have results for 3 gene sets
        assert result_ht.count() == 3

        # Check expected columns
        assert "beta" in result_ht.row
        assert "standard_error" in result_ht.row
        assert "p_value" in result_ht.row
        assert "odds_ratio" in result_ht.row
        assert "ci_lower" in result_ht.row
        assert "ci_upper" in result_ht.row

    def test_logistic_burden_test_with_covariates(self, hail_session):
        """Test logistic regression with covariates."""
        mt_burden = self.setup_burden_mt()

        result_ht = logistic_burden_test(
            mt_burden, phenotype_field="is_case", covariates=["PC1", "PC2", "sex"]
        )

        assert result_ht.count() == 3

    def test_logistic_burden_test_missing_phenotype(self, hail_session):
        """Test error with missing phenotype field."""
        mt_burden = self.setup_burden_mt()

        with pytest.raises(ValueError, match="Phenotype field"):
            logistic_burden_test(mt_burden, phenotype_field="nonexistent")

    def test_logistic_burden_test_missing_covariate(self, hail_session):
        """Test error with missing covariate."""
        mt_burden = self.setup_burden_mt()

        with pytest.raises(ValueError, match="Covariate"):
            logistic_burden_test(
                mt_burden, phenotype_field="is_case", covariates=["nonexistent"]
            )


@pytest.mark.hail
class TestLinearBurdenTest:
    """Tests for linear_burden_test function."""

    def setup_burden_mt_continuous(self):
        """Create burden MT with continuous phenotype."""
        mt = hl.utils.range_matrix_table(3, 100)
        mt = mt.annotate_rows(gene_set_name=hl.str(mt.row_idx))
        mt = mt.annotate_entries(burden=hl.int(hl.rand_unif(0, 5)))

        # Add continuous phenotype
        mt = mt.annotate_cols(
            cognitive_score=hl.rand_norm(mean=100, sd=15),
            age=hl.rand_unif(20, 80),
            sex=mt.col_idx % 2,
        )

        return mt

    def test_linear_burden_test_basic(self, hail_session):
        """Test basic linear regression."""
        mt_burden = self.setup_burden_mt_continuous()

        result_ht = linear_burden_test(mt_burden, phenotype_field="cognitive_score")

        assert result_ht.count() == 3

        # Check expected columns
        assert "beta" in result_ht.row
        assert "standard_error" in result_ht.row
        assert "p_value" in result_ht.row

    def test_linear_burden_test_with_covariates(self, hail_session):
        """Test linear regression with covariates."""
        mt_burden = self.setup_burden_mt_continuous()

        result_ht = linear_burden_test(
            mt_burden, phenotype_field="cognitive_score", covariates=["age", "sex"]
        )

        assert result_ht.count() == 3


@pytest.mark.hail
class TestRunBurdenAnalysis:
    """Tests for run_burden_analysis pipeline."""

    def setup_cohort_mt(self):
        """Create a realistic cohort MatrixTable."""
        # Create MT with 30 variants, 50 samples
        mt = hl.utils.range_matrix_table(30, 50)

        # Rename column key to 's' and make it string
        mt = mt.key_cols_by(s=hl.str(mt.col_idx))

        # Add variant annotations
        mt = mt.annotate_rows(
            SYMBOL=hl.if_else(
                mt.row_idx < 10,
                "BRCA1",
                hl.if_else(
                    mt.row_idx < 20, "TP53", hl.if_else(mt.row_idx < 25, "EGFR", "MYC")
                ),
            ),
            gnomad_af=hl.rand_unif(0, 0.05),
            cadd_phred=hl.rand_unif(15, 30),
        )

        # Add genotypes
        mt = mt.annotate_entries(
            GT=hl.call(hl.int(hl.rand_bool(0.3)), hl.int(hl.rand_bool(0.1))),
            GQ=30,
            DP=20,
        )

        return mt

    def setup_phenotype_ht(self):
        """Create a phenotype table."""
        # Create table with 50 samples
        ht = hl.utils.range_table(50)
        ht = ht.annotate(
            s=hl.str(ht.idx),
            is_case=ht.idx < 25,  # 25 cases, 25 controls
            PC1=hl.rand_norm(),
            PC2=hl.rand_norm(),
            sex=ht.idx % 2,
        )
        ht = ht.key_by("s")

        return ht

    def test_run_burden_analysis_binary(self, hail_session):
        """Test complete burden analysis pipeline with binary phenotype."""
        cohort_mt = self.setup_cohort_mt()
        phenotype_ht = self.setup_phenotype_ht()

        gene_sets = {
            "cancer_genes": ["BRCA1", "TP53"],
            "oncogenes": ["EGFR", "MYC"],
        }

        result_ht = run_burden_analysis(
            cohort_mt=cohort_mt,
            gene_sets=gene_sets,
            phenotype_ht=phenotype_ht,
            phenotype_field="is_case",
            phenotype_type="binary",
            max_af=0.1,
            min_cadd=None,
        )

        # Should have results for 2 gene sets
        assert result_ht.count() == 2

        # Check columns
        assert "gene_set_name" in result_ht.row
        assert "beta" in result_ht.row
        assert "p_value" in result_ht.row
        assert "odds_ratio" in result_ht.row

    def test_run_burden_analysis_with_covariates(self, hail_session):
        """Test burden analysis with covariates."""
        cohort_mt = self.setup_cohort_mt()
        phenotype_ht = self.setup_phenotype_ht()

        gene_sets = {"cancer_genes": ["BRCA1", "TP53"]}

        result_ht = run_burden_analysis(
            cohort_mt=cohort_mt,
            gene_sets=gene_sets,
            phenotype_ht=phenotype_ht,
            phenotype_field="is_case",
            covariate_fields=["PC1", "PC2"],
            phenotype_type="binary",
            max_af=0.1,
            min_cadd=None,
        )

        assert result_ht.count() == 1

    def test_run_burden_analysis_continuous(self, hail_session):
        """Test burden analysis with continuous phenotype."""
        cohort_mt = self.setup_cohort_mt()

        # Create continuous phenotype
        ht = hl.utils.range_table(50)
        ht = ht.annotate(
            s=hl.str(ht.idx),
            cognitive_score=hl.rand_norm(mean=100, sd=15),
            age=hl.rand_unif(20, 80),
        )
        ht = ht.key_by("s")

        gene_sets = {"cancer_genes": ["BRCA1", "TP53"]}

        result_ht = run_burden_analysis(
            cohort_mt=cohort_mt,
            gene_sets=gene_sets,
            phenotype_ht=ht,
            phenotype_field="cognitive_score",
            phenotype_type="continuous",
            max_af=0.1,
            min_cadd=None,
        )

        assert result_ht.count() == 1
        # Linear regression doesn't have odds_ratio
        assert "odds_ratio" not in result_ht.row

    def test_run_burden_analysis_invalid_phenotype_type(self, hail_session):
        """Test error with invalid phenotype type."""
        cohort_mt = self.setup_cohort_mt()
        phenotype_ht = self.setup_phenotype_ht()

        gene_sets = {"set1": ["BRCA1"]}

        with pytest.raises(ValueError, match="Invalid phenotype_type"):
            run_burden_analysis(
                cohort_mt=cohort_mt,
                gene_sets=gene_sets,
                phenotype_ht=phenotype_ht,
                phenotype_field="is_case",
                phenotype_type="invalid",
            )


@pytest.mark.hail
class TestBurdenAnalysisIntegration:
    """Integration tests for burden analysis."""

    def test_complete_workflow(self, hail_session):
        """Test complete burden analysis workflow."""
        # Create realistic data
        mt = hl.utils.range_matrix_table(40, 60)
        mt = mt.key_cols_by(s=hl.str(mt.col_idx))

        # Add annotations
        mt = mt.annotate_rows(
            SYMBOL=hl.if_else(
                mt.row_idx < 15,
                "TREM2",
                hl.if_else(
                    mt.row_idx < 25,
                    "CD33",
                    hl.if_else(mt.row_idx < 35, "ABI3", "OTHER"),
                ),
            ),
            gnomad_af=hl.rand_unif(0, 0.02),
            cadd_phred=hl.rand_unif(18, 35),
        )

        mt = mt.annotate_entries(
            GT=hl.call(hl.int(hl.rand_bool(0.2)), hl.int(hl.rand_bool(0.05))),
            GQ=25,
            DP=15,
        )

        # Phenotypes
        ht = hl.utils.range_table(60)
        ht = ht.annotate(
            s=hl.str(ht.idx),
            is_case=ht.idx < 30,
            PC1=hl.rand_norm(),
            PC2=hl.rand_norm(),
            PC3=hl.rand_norm(),
        )
        ht = ht.key_by("s")

        # Gene sets (simulate microglia genes)
        gene_sets = {
            "microglia": ["TREM2", "CD33", "ABI3"],
            "control_set": ["OTHER"],
        }

        # Run analysis
        result_ht = run_burden_analysis(
            cohort_mt=mt,
            gene_sets=gene_sets,
            phenotype_ht=ht,
            phenotype_field="is_case",
            covariate_fields=["PC1", "PC2", "PC3"],
            phenotype_type="binary",
            max_af=0.01,
            min_cadd=20.0,
            genotype_aggregation="hets",
        )

        # Convert to pandas for inspection
        df = result_ht.to_pandas()

        # Check structure
        assert len(df) == 2
        assert "gene_set_name" in df.columns
        assert "beta" in df.columns
        assert "p_value" in df.columns
        assert "odds_ratio" in df.columns

        # Check p-values are valid
        assert (df["p_value"] >= 0).all()
        assert (df["p_value"] <= 1).all()

    def test_genotype_aggregation_methods(self, hail_session):
        """Test all genotype aggregation methods work."""
        mt = hl.utils.range_matrix_table(20, 30)
        mt = mt.key_cols_by(s=hl.str(mt.col_idx))

        mt = mt.annotate_rows(
            SYMBOL=hl.if_else(mt.row_idx < 10, "GENE1", "GENE2"),
            gnomad_af=hl.rand_unif(0, 0.02),
        )

        mt = mt.annotate_entries(
            GT=hl.call(hl.int(hl.rand_bool(0.3)), hl.int(hl.rand_bool(0.1))),
            GQ=25,
            DP=15,
        )

        ht = hl.utils.range_table(30)
        ht = ht.annotate(
            s=hl.str(ht.idx),
            is_case=ht.idx < 15,
        )
        ht = ht.key_by("s")

        gene_sets = {"test_set": ["GENE1", "GENE2"]}

        for method in ["hets", "homs", "chets", "homs_chets"]:
            result_ht = run_burden_analysis(
                cohort_mt=mt,
                gene_sets=gene_sets,
                phenotype_ht=ht,
                phenotype_field="is_case",
                phenotype_type="binary",
                genotype_aggregation=method,
                max_af=0.1,
                min_cadd=None,
            )

            assert result_ht.count() == 1
