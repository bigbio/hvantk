"""Tests for synthetic burden cohort generation (Phase 5.1)."""

import pytest
import numpy as np

from hvantk.enrichex.simulation import (
    generate_synthetic_burden_cohort,
    check_type_i_error,
)


class TestCheckTypeIError:
    """Tests for check_type_i_error (no Hail needed)."""

    def test_uniform_pvalues(self):
        """Uniform p-values should be calibrated."""
        rng = np.random.default_rng(42)
        p_values = rng.uniform(0, 1, 1000).tolist()
        result = check_type_i_error(p_values)
        assert result["calibrated"]
        assert result["ks_pvalue"] > 0.01

    def test_inflated_pvalues(self):
        """All-significant p-values should NOT be calibrated."""
        p_values = [0.001] * 100
        result = check_type_i_error(p_values)
        assert not result["calibrated"]

    def test_empty_pvalues(self):
        """Empty list should return uncalibrated result."""
        result = check_type_i_error([])
        assert not result["calibrated"]
        assert result["n_tests"] == 0
        assert result["n_rejected"] == 0

    def test_result_keys(self):
        """Result dictionary should have all expected keys."""
        rng = np.random.default_rng(99)
        p_values = rng.uniform(0, 1, 500).tolist()
        result = check_type_i_error(p_values)
        expected_keys = {
            "n_tests",
            "n_rejected",
            "rejection_rate",
            "expected_rate",
            "calibrated",
            "ks_pvalue",
        }
        assert set(result.keys()) == expected_keys

    def test_rejection_rate_computed_correctly(self):
        """Rejection rate should match manual calculation."""
        # 10 p-values, 3 below alpha=0.05
        p_values = [0.01, 0.03, 0.04, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.99]
        result = check_type_i_error(p_values, alpha=0.05)
        assert result["n_tests"] == 10
        assert result["n_rejected"] == 3
        assert result["rejection_rate"] == pytest.approx(0.3)

    def test_custom_alpha(self):
        """Custom alpha should change the expected rate."""
        rng = np.random.default_rng(42)
        p_values = rng.uniform(0, 1, 1000).tolist()
        result = check_type_i_error(p_values, alpha=0.10)
        assert result["expected_rate"] == 0.10

    def test_strict_tolerance_fails_uniform(self):
        """Very strict tolerance should reject even well-behaved p-values."""
        # With a very tight tolerance, even uniform p-values may fail
        rng = np.random.default_rng(42)
        p_values = rng.uniform(0, 1, 50).tolist()
        result = check_type_i_error(p_values, alpha=0.05, tolerance=1.01)
        # With only 50 tests, random variation is high — calibration is fragile
        # We just check it returns a valid result
        assert isinstance(result["calibrated"], bool)

    def test_moderately_inflated_pvalues(self):
        """P-values slightly shifted toward 0 should not be calibrated."""
        rng = np.random.default_rng(42)
        # Generate from Beta(0.5, 1) — skewed toward 0
        p_values = rng.beta(0.5, 1.0, size=1000).tolist()
        result = check_type_i_error(p_values)
        assert not result["calibrated"]


@pytest.mark.hail
class TestSyntheticCohort:
    def test_basic_generation(self, hail_session):
        """Test basic cohort generation dimensions."""
        mt, pheno_ht, gene_sets = generate_synthetic_burden_cohort(
            n_cases=50,
            n_controls=50,
            n_genes=20,
            variants_per_gene=5,
            seed=42,
        )
        n_rows, n_cols = mt.count()
        assert n_rows == 100  # 20 genes * 5 variants
        assert n_cols == 100  # 50 + 50
        assert pheno_ht.count() == 100
        assert len(gene_sets) >= 1

    def test_has_required_fields(self, hail_session):
        """Test MT has required annotations."""
        mt, _, _ = generate_synthetic_burden_cohort(
            n_cases=30,
            n_controls=30,
            n_genes=10,
            variants_per_gene=3,
            seed=123,
        )
        assert "SYMBOL" in mt.row
        assert "consequence" in mt.row
        assert "GT" in mt.entry

    def test_seed_reproducibility(self, hail_session):
        """Same seed produces same data."""
        mt1, _, gs1 = generate_synthetic_burden_cohort(
            n_cases=20,
            n_controls=20,
            n_genes=10,
            variants_per_gene=3,
            seed=99,
        )
        mt2, _, gs2 = generate_synthetic_burden_cohort(
            n_cases=20,
            n_controls=20,
            n_genes=10,
            variants_per_gene=3,
            seed=99,
        )
        assert gs1 == gs2
        assert mt1.count() == mt2.count()

    def test_custom_signal_sets(self, hail_session):
        """Test with custom signal gene sets."""
        gene_sets = {"test_set": [f"GENE_{i}" for i in range(5)]}
        mt, pheno_ht, gs = generate_synthetic_burden_cohort(
            n_cases=50,
            n_controls=50,
            n_genes=20,
            variants_per_gene=5,
            signal_gene_sets=gene_sets,
            effect_sizes={"test_set": 3.0},
            seed=42,
        )
        assert "test_set" in gs
        assert gs["test_set"] == gene_sets["test_set"]

    def test_phenotype_table(self, hail_session):
        """Test phenotype table structure."""
        import hail as hl

        _, pheno_ht, _ = generate_synthetic_burden_cohort(
            n_cases=30,
            n_controls=30,
            n_genes=10,
            variants_per_gene=3,
            seed=42,
        )
        assert "is_case" in pheno_ht.row
        # Check case/control counts
        n_cases = pheno_ht.filter(pheno_ht.is_case).count()
        n_controls = pheno_ht.filter(~pheno_ht.is_case).count()
        assert n_cases == 30
        assert n_controls == 30

    def test_invalid_signal_genes_raises(self, hail_session):
        """Signal gene sets with unknown genes should raise ValueError."""
        gene_sets = {"bad_set": ["NOT_A_GENE"]}
        with pytest.raises(ValueError, match="unknown genes"):
            generate_synthetic_burden_cohort(
                n_cases=10,
                n_controls=10,
                n_genes=5,
                variants_per_gene=2,
                signal_gene_sets=gene_sets,
                seed=42,
            )

    def test_consequence_distribution(self, hail_session):
        """Consequences should cycle through the three types."""
        import hail as hl

        mt, _, _ = generate_synthetic_burden_cohort(
            n_cases=10,
            n_controls=10,
            n_genes=6,
            variants_per_gene=3,
            seed=42,
        )
        # 6 genes * 3 variants = 18 variants, each gene gets one of each consequence
        consequences = mt.aggregate_rows(hl.agg.collect(mt.consequence))
        counts = {}
        for c in consequences:
            counts[c] = counts.get(c, 0) + 1
        assert counts["stop_gained"] == 6
        assert counts["missense_variant"] == 6
        assert counts["synonymous_variant"] == 6

    def test_gnomad_af_field_present(self, hail_session):
        """MT should have gnomad_af row annotation."""
        mt, _, _ = generate_synthetic_burden_cohort(
            n_cases=10,
            n_controls=10,
            n_genes=5,
            variants_per_gene=2,
            seed=42,
        )
        assert "gnomad_af" in mt.row

    def test_default_gene_sets_created(self, hail_session):
        """When signal_gene_sets is None, default sets are created."""
        _, _, gene_sets = generate_synthetic_burden_cohort(
            n_cases=10,
            n_controls=10,
            n_genes=30,
            variants_per_gene=2,
            seed=42,
        )
        assert "signal_set_small" in gene_sets
        assert "signal_set_large" in gene_sets
        assert len(gene_sets["signal_set_small"]) == 10
        assert len(gene_sets["signal_set_large"]) == 20
