"""
Tests for burden analysis module (Hail-native regression).
"""

import pytest
import hail as hl
import pandas as pd

import numpy as np

from hvantk.enrichex.burden import (
    VariantFilter,
    _compute_test_statistic,
    _linear_t_statistic,
    _logistic_z_statistic,
    build_variant_classes_from_presets,
    permutation_burden_test,
    compute_geneset_burden_mt,
    compute_per_gene_burden_mt,
    linear_burden_test,
    logistic_burden_test,
    run_burden_analysis,
    run_stratified_burden_analysis,
)


@pytest.mark.hail
class TestVariantFilter:
    """Tests for VariantFilter dataclass."""

    def test_variant_filter_creation(self):
        """Test basic VariantFilter creation."""
        vf = VariantFilter(
            max_af=0.01,
            min_score=20.0,
            consequences=["missense_variant", "frameshift_variant"],
        )

        assert vf.max_af == 0.01
        assert vf.min_score == 20.0
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

        vf = VariantFilter(max_af=0.01, min_score=20.0)
        expr = vf.to_hail_expr(mt)

        # Should be a boolean expression
        assert isinstance(expr, hl.expr.BooleanExpression)

    def test_variant_filter_missing_fields(self, hail_session):
        """Test filter with missing fields."""
        mt = hl.utils.range_matrix_table(10, 5)

        vf = VariantFilter(max_af=0.01, min_score=20.0)

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
            mt, gene_sets, gene_field="SYMBOL"
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

        for method in ["hets", "homs", "multi_het", "homs_multi_het"]:
            mt_burden = compute_geneset_burden_mt(
                mt,
                gene_sets,
                gene_field="SYMBOL",
                genotype_aggregation=method,
            )

            assert mt_burden.count_rows() == 1

    def test_deprecated_chets_alias(self, hail_session):
        """Test that deprecated 'chets' alias works with warning."""
        import warnings

        mt = self.setup_test_mt()
        gene_sets = {"set1": ["GENE1", "GENE2"]}

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            mt_burden = compute_geneset_burden_mt(
                mt,
                gene_sets,
                gene_field="SYMBOL",
                genotype_aggregation="chets",
            )
            assert mt_burden.count_rows() == 1
            assert any("deprecated" in str(warning.message).lower() for warning in w)

    def test_deprecated_homs_chets_alias(self, hail_session):
        """Test that deprecated 'homs_chets' alias works with warning."""
        import warnings

        mt = self.setup_test_mt()
        gene_sets = {"set1": ["GENE1", "GENE2"]}

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            mt_burden = compute_geneset_burden_mt(
                mt,
                gene_sets,
                gene_field="SYMBOL",
                genotype_aggregation="homs_chets",
            )
            assert mt_burden.count_rows() == 1
            assert any("deprecated" in str(warning.message).lower() for warning in w)

    def test_compute_geneset_burden_no_filtering(self, hail_session):
        """Test burden computation with no filtering (pre-filtered MT)."""
        mt = self.setup_test_mt()
        gene_sets = {"set1": ["GENE1", "GENE2"], "set2": ["GENE3"]}

        mt_burden = compute_geneset_burden_mt(
            mt, gene_sets, gene_field="SYMBOL"
        )

        n_rows, n_cols = mt_burden.count()
        assert n_rows == 2
        assert n_cols == 10

    def test_compute_geneset_burden_with_variant_filter(self, hail_session):
        """Test burden computation with VariantFilter."""
        mt = self.setup_test_mt()
        gene_sets = {"set1": ["GENE1", "GENE2"]}

        vf = VariantFilter(max_af=0.1, min_score=None, pass_only=False, min_gq=0, min_dp=0)
        mt_burden = compute_geneset_burden_mt(
            mt, gene_sets, gene_field="SYMBOL", variant_filter=vf
        )

        assert mt_burden.count_rows() == 1

    def test_compute_geneset_burden_no_variants(self, hail_session):
        """Test with very restrictive filters (no variants pass) returns None."""
        mt = self.setup_test_mt()
        gene_sets = {"set1": ["GENE1"]}

        vf = VariantFilter(max_af=0.0001, min_score=None, pass_only=False, min_gq=0, min_dp=0)
        result = compute_geneset_burden_mt(
            mt,
            gene_sets,
            gene_field="SYMBOL",
            variant_filter=vf,
        )
        assert result is None

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

        # Run analysis with VariantFilter for variant-class selection
        vf = VariantFilter(
            max_af=0.01,
            min_score=20.0,
            pass_only=False,
            min_gq=0,
            min_dp=0,
        )
        result_ht = run_burden_analysis(
            cohort_mt=mt,
            gene_sets=gene_sets,
            phenotype_ht=ht,
            phenotype_field="is_case",
            covariate_fields=["PC1", "PC2", "PC3"],
            phenotype_type="binary",
            genotype_aggregation="hets",
            variant_filter=vf,
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

        for method in ["hets", "homs", "multi_het", "homs_multi_het"]:
            result_ht = run_burden_analysis(
                cohort_mt=mt,
                gene_sets=gene_sets,
                phenotype_ht=ht,
                phenotype_field="is_case",
                phenotype_type="binary",
                genotype_aggregation=method,
            )

            assert result_ht.count() == 1


@pytest.mark.hail
class TestBuildVariantClassesFromPresets:
    """Tests for build_variant_classes_from_presets."""

    def test_build_from_presets_default(self):
        """Test building variant classes from preset names."""
        classes = build_variant_classes_from_presets(
            ["lof", "missense_constrained", "synonymous"]
        )

        assert len(classes) == 3
        assert "lof" in classes
        assert "missense_constrained" in classes
        assert "synonymous" in classes

        # Check LoF consequences
        assert "stop_gained" in classes["lof"].consequences
        assert "frameshift_variant" in classes["lof"].consequences

        # Check missense has CADD threshold
        assert classes["missense_constrained"].min_score == 25.0
        assert "missense_variant" in classes["missense_constrained"].consequences

        # Check synonymous
        assert "synonymous_variant" in classes["synonymous"].consequences
        assert classes["synonymous"].min_score is None

    def test_build_from_presets_permissive_defaults(self):
        """Test that presets without base_filter use permissive defaults."""
        classes = build_variant_classes_from_presets(["lof"])

        vf = classes["lof"]
        assert vf.max_af == 1.0
        assert vf.pass_only is False
        assert vf.min_gq == 0
        assert vf.min_dp == 0

    def test_build_from_presets_with_base_filter(self):
        """Test that base_filter settings are inherited."""
        base = VariantFilter(
            max_af=0.001,
            min_score=15.0,
            pass_only=True,
            min_gq=20,
            min_dp=10,
            af_field="af",
            score_field="cadd",
            consequence_field="csq",
        )
        classes = build_variant_classes_from_presets(
            ["lof", "missense_constrained"], base_filter=base
        )

        # AF and field names inherited from base
        assert classes["lof"].max_af == 0.001
        assert classes["lof"].af_field == "af"
        assert classes["lof"].consequence_field == "csq"

        # LoF inherits base min_score (no preset override)
        assert classes["lof"].min_score == 15.0

        # Missense overrides min_score from preset
        assert classes["missense_constrained"].min_score == 25.0
        assert classes["missense_constrained"].max_af == 0.001

    def test_build_from_presets_unknown_class(self):
        """Test error for unknown variant class name."""
        with pytest.raises(ValueError, match="Unknown variant class preset"):
            build_variant_classes_from_presets(["nonexistent"])


@pytest.mark.hail
class TestRunStratifiedBurdenAnalysis:
    """Tests for run_stratified_burden_analysis."""

    def setup_cohort_mt(self):
        """Create a cohort MT with consequence annotations."""
        mt = hl.utils.range_matrix_table(40, 50)
        mt = mt.key_cols_by(s=hl.str(mt.col_idx))

        # Assign consequences: first 15 = stop_gained, next 10 = missense,
        # next 10 = synonymous, last 5 = intron
        mt = mt.annotate_rows(
            SYMBOL=hl.if_else(
                mt.row_idx < 20,
                "GENE1",
                hl.if_else(mt.row_idx < 30, "GENE2", "GENE3"),
            ),
            consequence=hl.if_else(
                mt.row_idx < 15,
                "stop_gained",
                hl.if_else(
                    mt.row_idx < 25,
                    "missense_variant",
                    hl.if_else(mt.row_idx < 35, "synonymous_variant", "intron_variant"),
                ),
            ),
            gnomad_af=hl.rand_unif(0, 0.02),
            cadd_phred=hl.rand_unif(10, 35),
        )

        mt = mt.annotate_entries(
            GT=hl.call(hl.int(hl.rand_bool(0.3)), hl.int(hl.rand_bool(0.1))),
            GQ=25,
            DP=15,
        )
        return mt

    def setup_phenotype_ht(self):
        """Create a phenotype table."""
        ht = hl.utils.range_table(50)
        ht = ht.annotate(
            s=hl.str(ht.idx),
            is_case=ht.idx < 25,
            PC1=hl.rand_norm(),
        )
        ht = ht.key_by("s")
        return ht

    def test_stratified_basic(self, hail_session):
        """Test stratified burden with two variant classes."""
        mt = self.setup_cohort_mt()
        ht = self.setup_phenotype_ht()

        gene_sets = {"set1": ["GENE1", "GENE2"], "set2": ["GENE3"]}

        variant_classes = {
            "lof": VariantFilter(
                consequences=["stop_gained"],
                pass_only=False,
                min_gq=0,
                min_dp=0,
                min_score=None,
            ),
            "synonymous": VariantFilter(
                consequences=["synonymous_variant"],
                pass_only=False,
                min_gq=0,
                min_dp=0,
                min_score=None,
            ),
        }

        results = run_stratified_burden_analysis(
            cohort_mt=mt,
            gene_sets=gene_sets,
            phenotype_ht=ht,
            variant_classes=variant_classes,
            phenotype_field="is_case",
        )

        # Should have results for both classes
        assert "lof" in results
        assert "synonymous" in results

        # Each should have gene set results
        for cls, result_ht in results.items():
            df = result_ht.to_pandas()
            assert "variant_class" in df.columns
            assert (df["variant_class"] == cls).all()
            assert "gene_set_name" in df.columns
            assert "p_value" in df.columns

    def test_stratified_with_presets(self, hail_session):
        """Test stratified burden using preset classes."""
        mt = self.setup_cohort_mt()
        ht = self.setup_phenotype_ht()

        gene_sets = {"set1": ["GENE1", "GENE2"]}

        # Use presets with field name override
        base = VariantFilter(
            max_af=1.0,
            min_score=None,
            pass_only=False,
            min_gq=0,
            min_dp=0,
            consequence_field="consequence",
        )
        variant_classes = build_variant_classes_from_presets(
            ["lof", "synonymous"], base_filter=base
        )

        results = run_stratified_burden_analysis(
            cohort_mt=mt,
            gene_sets=gene_sets,
            phenotype_ht=ht,
            variant_classes=variant_classes,
            phenotype_field="is_case",
        )

        assert len(results) >= 1  # at least one class should have results

    def test_stratified_skips_empty_class(self, hail_session):
        """Test that classes with no qualifying variants are skipped."""
        mt = self.setup_cohort_mt()
        ht = self.setup_phenotype_ht()

        gene_sets = {"set1": ["GENE1", "GENE2"]}

        variant_classes = {
            "lof": VariantFilter(
                consequences=["stop_gained"],
                pass_only=False,
                min_gq=0,
                min_dp=0,
                min_score=None,
            ),
            "nonexistent_class": VariantFilter(
                consequences=["this_consequence_does_not_exist"],
                pass_only=False,
                min_gq=0,
                min_dp=0,
                min_score=None,
            ),
        }

        results = run_stratified_burden_analysis(
            cohort_mt=mt,
            gene_sets=gene_sets,
            phenotype_ht=ht,
            variant_classes=variant_classes,
            phenotype_field="is_case",
        )

        # Should have lof but not nonexistent_class
        assert "lof" in results
        assert "nonexistent_class" not in results

    def test_stratified_variant_class_field_in_output(self, hail_session):
        """Test that variant_class column is present in results."""
        mt = self.setup_cohort_mt()
        ht = self.setup_phenotype_ht()

        gene_sets = {"set1": ["GENE1", "GENE2"]}

        variant_classes = {
            "lof": VariantFilter(
                consequences=["stop_gained"],
                pass_only=False,
                min_gq=0,
                min_dp=0,
                min_score=None,
            ),
        }

        results = run_stratified_burden_analysis(
            cohort_mt=mt,
            gene_sets=gene_sets,
            phenotype_ht=ht,
            variant_classes=variant_classes,
            phenotype_field="is_case",
        )

        df = results["lof"].to_pandas()
        assert "variant_class" in df.columns
        assert (df["variant_class"] == "lof").all()
        assert "p_value" in df.columns
        assert (df["p_value"] >= 0).all()
        assert (df["p_value"] <= 1).all()


@pytest.mark.hail
class TestNormalizeByLength:
    """Tests for gene-length normalization in burden computation."""

    def setup_test_mt(self):
        """Create a test MT with genes of different variant counts."""
        # 30 variants: GENE1 has 15, GENE2 has 10, GENE3 has 5
        mt = hl.utils.range_matrix_table(30, 20)
        mt = mt.key_cols_by(s=hl.str(mt.col_idx))

        mt = mt.annotate_rows(
            SYMBOL=hl.if_else(
                mt.row_idx < 15,
                "GENE1",
                hl.if_else(mt.row_idx < 25, "GENE2", "GENE3"),
            ),
        )

        mt = mt.annotate_entries(
            GT=hl.call(hl.int(hl.rand_bool(0.3)), hl.int(hl.rand_bool(0.1))),
        )
        return mt

    def test_normalize_with_provided_lengths(self, hail_session):
        """Test normalization using provided CDS lengths."""
        mt = self.setup_test_mt()
        gene_sets = {"set1": ["GENE1", "GENE2", "GENE3"]}
        gene_lengths = {"GENE1": 3000, "GENE2": 1000, "GENE3": 500}

        mt_burden = compute_geneset_burden_mt(
            mt,
            gene_sets,
            gene_field="SYMBOL",
            normalize_by_length=True,
            gene_lengths=gene_lengths,
        )

        n_rows, n_cols = mt_burden.count()
        assert n_rows == 1
        assert n_cols == 20

        # Burden should be float (not integer) when normalized
        burden_vals = mt_burden.burden.collect()
        assert all(isinstance(v, float) for v in burden_vals)

    def test_normalize_with_proxy(self, hail_session):
        """Test normalization using variant site count proxy."""
        mt = self.setup_test_mt()
        gene_sets = {"set1": ["GENE1", "GENE2", "GENE3"]}

        mt_burden = compute_geneset_burden_mt(
            mt,
            gene_sets,
            gene_field="SYMBOL",
            normalize_by_length=True,
            gene_lengths=None,
        )

        n_rows, n_cols = mt_burden.count()
        assert n_rows == 1
        assert n_cols == 20

        # Should still produce valid burden values
        burden_vals = mt_burden.burden.collect()
        assert all(v is not None for v in burden_vals)

    def test_normalize_changes_burden_values(self, hail_session):
        """Test that normalization produces different values than raw."""
        mt = self.setup_test_mt()
        gene_sets = {"set1": ["GENE1", "GENE2", "GENE3"]}
        gene_lengths = {"GENE1": 3000, "GENE2": 1000, "GENE3": 500}

        mt_raw = compute_geneset_burden_mt(mt, gene_sets, gene_field="SYMBOL")
        mt_norm = compute_geneset_burden_mt(
            mt,
            gene_sets,
            gene_field="SYMBOL",
            normalize_by_length=True,
            gene_lengths=gene_lengths,
        )

        raw_vals = mt_raw.burden.collect()
        norm_vals = mt_norm.burden.collect()

        # Dimensions should match
        assert len(raw_vals) == len(norm_vals)

        # At least some values should differ (genes have different lengths)
        raw_sum = sum(float(v) for v in raw_vals)
        norm_sum = sum(float(v) for v in norm_vals)
        assert raw_sum != norm_sum

    def test_normalize_threads_through_run_burden_analysis(self, hail_session):
        """Test normalization works through the full pipeline."""
        mt = self.setup_test_mt()
        ht = hl.utils.range_table(20)
        ht = ht.annotate(s=hl.str(ht.idx), is_case=ht.idx < 10)
        ht = ht.key_by("s")

        gene_sets = {"set1": ["GENE1", "GENE2"]}
        gene_lengths = {"GENE1": 3000, "GENE2": 1000, "GENE3": 500}

        result_ht = run_burden_analysis(
            cohort_mt=mt,
            gene_sets=gene_sets,
            phenotype_ht=ht,
            phenotype_field="is_case",
            phenotype_type="binary",
            normalize_by_length=True,
            gene_lengths=gene_lengths,
        )

        assert result_ht.count() == 1
        assert "p_value" in result_ht.row

    def test_normalize_in_permutation_test(self, hail_session):
        """Test normalization in permutation burden test."""
        mt = hl.utils.range_matrix_table(100, 80)
        mt = mt.key_cols_by(s=hl.str(mt.col_idx))
        gene_names = [f"GENE{i + 1}" for i in range(10)]
        mt = mt.annotate_rows(
            SYMBOL=hl.literal(gene_names)[mt.row_idx // 10],
        )
        mt = mt.annotate_entries(
            GT=hl.call(
                hl.int(hl.rand_bool(0.3)), hl.int(hl.rand_bool(0.1))
            ),
        )

        ht = hl.utils.range_table(80)
        ht = ht.annotate(s=hl.str(ht.idx), is_case=ht.idx < 40)
        ht = ht.key_by("s")

        gene_sets = {"set1": ["GENE1", "GENE2"]}
        gene_lengths = {f"GENE{i + 1}": (i + 1) * 1000 for i in range(10)}

        result_df = permutation_burden_test(
            cohort_mt=mt,
            gene_sets=gene_sets,
            phenotype_ht=ht,
            phenotype_field="is_case",
            n_permutations=50,
            seed=42,
            normalize_by_length=True,
            gene_lengths=gene_lengths,
        )

        assert len(result_df) == 1
        assert 0 < result_df.iloc[0]["empirical_p_value"] <= 1


class TestLocalRegressionHelpers:
    """Tests for local numpy-based regression helpers (no Hail needed)."""

    def test_linear_t_statistic_basic(self):
        """Test OLS t-statistic with known relationship."""
        rng = np.random.default_rng(42)
        n = 200
        x = rng.normal(size=n)
        y = 2.0 * x + rng.normal(size=n) * 0.5  # strong signal
        intercept = np.ones(n)
        X = np.column_stack([x, intercept])

        t_stat = _linear_t_statistic(X, y)
        assert not np.isnan(t_stat)
        assert abs(t_stat) > 5  # should be highly significant

    def test_linear_t_statistic_no_signal(self):
        """Test OLS t-statistic with no relationship."""
        rng = np.random.default_rng(42)
        n = 200
        x = rng.normal(size=n)
        y = rng.normal(size=n)  # no signal
        intercept = np.ones(n)
        X = np.column_stack([x, intercept])

        t_stat = _linear_t_statistic(X, y)
        assert not np.isnan(t_stat)
        assert abs(t_stat) < 5  # should not be significant

    def test_logistic_z_statistic_basic(self):
        """Test logistic z-statistic with a separable predictor."""
        rng = np.random.default_rng(42)
        n = 200
        x = rng.normal(size=n)
        prob = 1 / (1 + np.exp(-2 * x))
        y = (rng.random(n) < prob).astype(float)
        intercept = np.ones(n)
        X = np.column_stack([x, intercept])

        z_stat = _logistic_z_statistic(X, y)
        assert not np.isnan(z_stat)
        assert abs(z_stat) > 2  # should detect the signal

    def test_logistic_z_statistic_no_signal(self):
        """Test logistic z-statistic with no relationship."""
        rng = np.random.default_rng(42)
        n = 200
        x = rng.normal(size=n)
        y = (rng.random(n) < 0.5).astype(float)  # random labels
        intercept = np.ones(n)
        X = np.column_stack([x, intercept])

        z_stat = _logistic_z_statistic(X, y)
        assert not np.isnan(z_stat)

    def test_compute_test_statistic_dispatches(self):
        """Test that _compute_test_statistic dispatches correctly."""
        rng = np.random.default_rng(42)
        n = 100
        x = rng.normal(size=n)
        intercept = np.ones(n)
        X = np.column_stack([x, intercept])

        y_continuous = rng.normal(size=n)
        stat = _compute_test_statistic(X, y_continuous, "continuous")
        assert not np.isnan(stat)

        y_binary = (rng.random(n) < 0.5).astype(float)
        stat = _compute_test_statistic(X, y_binary, "binary")
        assert not np.isnan(stat)


@pytest.mark.hail
class TestComputePerGeneBurdenMt:
    """Tests for compute_per_gene_burden_mt."""

    def test_basic(self, hail_session):
        """Test per-gene burden computation."""
        mt = hl.utils.range_matrix_table(20, 10)
        mt = mt.annotate_rows(
            SYMBOL=hl.if_else(mt.row_idx < 10, "GENE1", "GENE2"),
        )
        mt = mt.annotate_entries(
            GT=hl.call(hl.int(hl.rand_bool(0.3)), hl.int(hl.rand_bool(0.1))),
        )

        mt_genes = compute_per_gene_burden_mt(mt, gene_field="SYMBOL")

        # Should have 2 genes (rows) and 10 samples (cols)
        n_rows, n_cols = mt_genes.count()
        assert n_rows == 2
        assert n_cols == 10

        # Should have expected entry fields
        assert "hets" in mt_genes.entry
        assert "homs" in mt_genes.entry
        assert "multi_het" in mt_genes.entry

    def test_with_variant_filter(self, hail_session):
        """Test per-gene burden with variant filter."""
        mt = hl.utils.range_matrix_table(20, 10)
        mt = mt.annotate_rows(
            SYMBOL=hl.if_else(mt.row_idx < 10, "GENE1", "GENE2"),
            consequence=hl.if_else(
                mt.row_idx < 10, "stop_gained", "synonymous_variant"
            ),
        )
        mt = mt.annotate_entries(
            GT=hl.call(hl.int(hl.rand_bool(0.3)), hl.int(hl.rand_bool(0.1))),
        )

        vf = VariantFilter(
            consequences=["stop_gained"],
            pass_only=False,
            min_gq=0,
            min_dp=0,
            min_score=None,
        )
        mt_genes = compute_per_gene_burden_mt(
            mt, gene_field="SYMBOL", variant_filter=vf
        )

        # Only GENE1 has stop_gained variants
        assert mt_genes.count_rows() == 1


@pytest.mark.hail
class TestPermutationBurdenTest:
    """Tests for permutation_burden_test."""

    def setup_cohort_mt(self):
        """Create a cohort MT with multiple genes.

        Uses 100 variants across 10 genes with 80 samples to provide
        enough background diversity for competitive permutation testing.
        """
        mt = hl.utils.range_matrix_table(100, 80)
        mt = mt.key_cols_by(s=hl.str(mt.col_idx))

        # 10 genes, 10 variants each
        gene_names = [f"GENE{i + 1}" for i in range(10)]
        mt = mt.annotate_rows(
            SYMBOL=hl.literal(gene_names)[mt.row_idx // 10],
        )

        mt = mt.annotate_entries(
            GT=hl.call(
                hl.int(hl.rand_bool(0.3)), hl.int(hl.rand_bool(0.1))
            ),
        )
        return mt

    def setup_phenotype_ht(self):
        """Create phenotype table."""
        ht = hl.utils.range_table(80)
        ht = ht.annotate(
            s=hl.str(ht.idx),
            is_case=ht.idx < 40,
            PC1=hl.rand_norm(),
        )
        ht = ht.key_by("s")
        return ht

    def test_permutation_basic(self, hail_session):
        """Test basic permutation burden test."""
        mt = self.setup_cohort_mt()
        ht = self.setup_phenotype_ht()

        gene_sets = {"set1": ["GENE1", "GENE2"]}

        result_df = permutation_burden_test(
            cohort_mt=mt,
            gene_sets=gene_sets,
            phenotype_ht=ht,
            phenotype_field="is_case",
            n_permutations=100,
            seed=42,
        )

        # Check result structure
        assert len(result_df) == 1
        assert "gene_set_name" in result_df.columns
        assert "observed_statistic" in result_df.columns
        assert "empirical_p_value" in result_df.columns
        assert "n_permutations" in result_df.columns
        assert "n_exceeded" in result_df.columns
        assert "n_genes_tested" in result_df.columns

        # Check values are valid
        row = result_df.iloc[0]
        assert row["gene_set_name"] == "set1"
        assert row["n_genes_tested"] == 2
        assert 0 < row["empirical_p_value"] <= 1
        assert row["n_permutations"] == 100

    def test_permutation_multiple_gene_sets(self, hail_session):
        """Test permutation test with multiple gene sets."""
        mt = self.setup_cohort_mt()
        ht = self.setup_phenotype_ht()

        gene_sets = {
            "set1": ["GENE1", "GENE2"],
            "set2": ["GENE3", "GENE4"],
        }

        result_df = permutation_burden_test(
            cohort_mt=mt,
            gene_sets=gene_sets,
            phenotype_ht=ht,
            phenotype_field="is_case",
            n_permutations=50,
            seed=42,
        )

        assert len(result_df) == 2
        assert set(result_df["gene_set_name"]) == {"set1", "set2"}

    def test_permutation_with_covariates(self, hail_session):
        """Test permutation test with covariates."""
        mt = self.setup_cohort_mt()
        ht = self.setup_phenotype_ht()

        gene_sets = {"set1": ["GENE1", "GENE2"]}

        result_df = permutation_burden_test(
            cohort_mt=mt,
            gene_sets=gene_sets,
            phenotype_ht=ht,
            phenotype_field="is_case",
            covariate_fields=["PC1"],
            n_permutations=50,
            seed=42,
        )

        assert len(result_df) == 1
        assert not np.isnan(result_df.iloc[0]["observed_statistic"])

    def test_permutation_reproducible_with_seed(self, hail_session):
        """Test that same seed gives same results."""
        mt = self.setup_cohort_mt()
        ht = self.setup_phenotype_ht()

        gene_sets = {"set1": ["GENE1", "GENE2"]}
        kwargs = dict(
            cohort_mt=mt,
            gene_sets=gene_sets,
            phenotype_ht=ht,
            phenotype_field="is_case",
            n_permutations=50,
            seed=123,
        )

        result1 = permutation_burden_test(**kwargs)
        result2 = permutation_burden_test(**kwargs)

        assert result1.iloc[0]["empirical_p_value"] == result2.iloc[0]["empirical_p_value"]
        assert result1.iloc[0]["observed_statistic"] == result2.iloc[0]["observed_statistic"]

    def test_permutation_length_matched(self, hail_session):
        """Test permutation test with length-matched sampling."""
        mt = self.setup_cohort_mt()
        ht = self.setup_phenotype_ht()

        gene_sets = {"set1": ["GENE1", "GENE2"]}

        # With length matching
        result_lm = permutation_burden_test(
            cohort_mt=mt,
            gene_sets=gene_sets,
            phenotype_ht=ht,
            phenotype_field="is_case",
            n_permutations=50,
            length_matched=True,
            seed=42,
        )

        # Without length matching
        result_no_lm = permutation_burden_test(
            cohort_mt=mt,
            gene_sets=gene_sets,
            phenotype_ht=ht,
            phenotype_field="is_case",
            n_permutations=50,
            length_matched=False,
            seed=42,
        )

        # Both should produce valid results
        assert len(result_lm) == 1
        assert len(result_no_lm) == 1
        assert 0 < result_lm.iloc[0]["empirical_p_value"] <= 1
        assert 0 < result_no_lm.iloc[0]["empirical_p_value"] <= 1

    def test_permutation_skips_missing_genes(self, hail_session):
        """Test that gene sets with no genes in MT are skipped."""
        mt = self.setup_cohort_mt()
        ht = self.setup_phenotype_ht()

        gene_sets = {"missing_set": ["NONEXISTENT1", "NONEXISTENT2"]}

        result_df = permutation_burden_test(
            cohort_mt=mt,
            gene_sets=gene_sets,
            phenotype_ht=ht,
            phenotype_field="is_case",
            n_permutations=50,
            seed=42,
        )

        # Should skip the gene set, returning empty DataFrame
        assert len(result_df) == 0


@pytest.mark.hail
class TestGracefulEmptyHandling:
    """Tests for Phase 3.1: graceful handling of empty/degenerate results."""

    def setup_cohort_mt(self):
        """Create a test cohort MatrixTable."""
        mt = hl.utils.range_matrix_table(20, 50)
        mt = mt.key_cols_by(s=hl.str(mt.col_idx))
        mt = mt.annotate_rows(
            SYMBOL=hl.if_else(
                mt.row_idx < 10,
                "GENE1",
                hl.if_else(mt.row_idx < 15, "GENE2", "GENE3"),
            ),
            gnomad_af=hl.rand_unif(0, 0.05),
            cadd_phred=hl.rand_unif(15, 30),
            consequence=hl.if_else(
                mt.row_idx % 3 == 0, "stop_gained", "missense_variant"
            ),
        )
        mt = mt.annotate_entries(
            GT=hl.call(hl.int(hl.rand_bool(0.3)), hl.int(hl.rand_bool(0.1))),
            GQ=30,
            DP=20,
        )
        return mt

    def setup_phenotype_ht(self):
        """Create a phenotype Hail Table matching the cohort."""
        ht = hl.utils.range_table(50)
        ht = ht.annotate(s=hl.str(ht.idx))
        ht = ht.key_by("s")
        ht = ht.annotate(
            is_case=hl.if_else(ht.idx < 25, 1, 0),
            PC1=hl.rand_norm(0, 1),
        )
        return ht

    def test_run_burden_no_variants_returns_none(self, hail_session):
        """run_burden_analysis returns None when no variants match gene sets."""
        mt = self.setup_cohort_mt()
        ht = self.setup_phenotype_ht()

        # Gene sets with genes not in MT
        gene_sets = {"set1": ["NONEXISTENT_GENE1", "NONEXISTENT_GENE2"]}

        result = run_burden_analysis(
            cohort_mt=mt,
            gene_sets=gene_sets,
            phenotype_ht=ht,
            phenotype_field="is_case",
        )
        assert result is None

    def test_run_burden_no_phenotype_overlap_returns_none(self, hail_session):
        """run_burden_analysis returns None when no samples have phenotype data."""
        mt = self.setup_cohort_mt()

        # Create phenotype table with non-matching sample IDs
        ht = hl.utils.range_table(10)
        ht = ht.annotate(s=hl.format("NOMATCH_%d", ht.idx))
        ht = ht.key_by("s")
        ht = ht.annotate(is_case=hl.int32(hl.if_else(ht.idx < 5, 1, 0)))

        gene_sets = {"set1": ["GENE1", "GENE2"]}

        result = run_burden_analysis(
            cohort_mt=mt,
            gene_sets=gene_sets,
            phenotype_ht=ht,
            phenotype_field="is_case",
        )
        assert result is None

    def test_run_burden_min_carriers_filters_gene_sets(self, hail_session):
        """min_carriers filters out gene sets with too few carriers."""
        mt = self.setup_cohort_mt()
        ht = self.setup_phenotype_ht()

        gene_sets = {"set1": ["GENE1", "GENE2"], "set2": ["GENE3"]}

        # With min_carriers=0 (default), all gene sets should be tested
        result_all = run_burden_analysis(
            cohort_mt=mt,
            gene_sets=gene_sets,
            phenotype_ht=ht,
            phenotype_field="is_case",
            min_carriers=0,
        )
        assert result_all is not None

        # With very high min_carriers, all gene sets should be filtered
        result_none = run_burden_analysis(
            cohort_mt=mt,
            gene_sets=gene_sets,
            phenotype_ht=ht,
            phenotype_field="is_case",
            min_carriers=999999,
        )
        assert result_none is None

    def test_stratified_all_empty_returns_empty_dict(self, hail_session):
        """run_stratified_burden_analysis returns {} when all classes are empty."""
        mt = self.setup_cohort_mt()
        ht = self.setup_phenotype_ht()

        gene_sets = {"set1": ["NONEXISTENT1", "NONEXISTENT2"]}

        variant_classes = {
            "lof": VariantFilter(
                consequences=["stop_gained"],
                pass_only=False,
                min_gq=0,
                min_dp=0,
                min_score=None,
            ),
        }

        results = run_stratified_burden_analysis(
            cohort_mt=mt,
            gene_sets=gene_sets,
            phenotype_ht=ht,
            variant_classes=variant_classes,
            phenotype_field="is_case",
        )
        assert len(results) == 0

    def test_permutation_no_phenotype_overlap_returns_empty_df(self, hail_session):
        """permutation_burden_test returns empty DataFrame with no sample overlap."""
        mt = self.setup_cohort_mt()

        # Non-matching phenotype table
        ht = hl.utils.range_table(10)
        ht = ht.annotate(s=hl.format("NOMATCH_%d", ht.idx))
        ht = ht.key_by("s")
        ht = ht.annotate(is_case=hl.int32(hl.if_else(ht.idx < 5, 1, 0)))

        gene_sets = {"set1": ["GENE1", "GENE2"]}

        result_df = permutation_burden_test(
            cohort_mt=mt,
            gene_sets=gene_sets,
            phenotype_ht=ht,
            phenotype_field="is_case",
            n_permutations=10,
            seed=42,
        )
        assert isinstance(result_df, pd.DataFrame)
        assert len(result_df) == 0
        assert "gene_set_name" in result_df.columns

    def test_compute_geneset_burden_mt_no_gene_overlap_returns_none(self, hail_session):
        """compute_geneset_burden_mt returns None when no genes match."""
        mt = self.setup_cohort_mt()
        gene_sets = {"set1": ["NO_SUCH_GENE"]}

        result = compute_geneset_burden_mt(mt, gene_sets, gene_field="SYMBOL")
        assert result is None
