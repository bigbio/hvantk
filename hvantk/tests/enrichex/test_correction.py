"""
Tests for multiple testing correction methods.
"""

import pytest
import numpy as np

from hvantk.utils.correction import apply_correction, fdr_threshold


class TestApplyCorrection:
    """Tests for apply_correction function."""

    def test_no_correction(self):
        """Test that 'none' method returns original p-values."""
        p_values = [0.001, 0.01, 0.05, 0.1]
        adjusted = apply_correction(p_values, method="none")

        assert adjusted == p_values

    def test_bonferroni_correction(self):
        """Test Bonferroni correction."""
        p_values = [0.001, 0.01, 0.05, 0.1]
        adjusted = apply_correction(p_values, method="bonferroni")

        # Each p-value multiplied by n_tests (4)
        expected = [0.004, 0.04, 0.2, 0.4]

        assert len(adjusted) == len(p_values)
        for adj, exp in zip(adjusted, expected):
            assert abs(adj - exp) < 1e-10

    def test_bonferroni_capped_at_one(self):
        """Test that Bonferroni correction caps values at 1.0."""
        p_values = [0.5, 0.6, 0.7]
        adjusted = apply_correction(p_values, method="bonferroni")

        # All should be capped at 1.0
        assert all(adj == 1.0 for adj in adjusted)

    def test_benjamini_hochberg_correction(self):
        """Test Benjamini-Hochberg FDR correction."""
        p_values = [0.001, 0.01, 0.05, 0.1]
        adjusted = apply_correction(p_values, method="benjamini-hochberg")

        # Check basic properties
        assert len(adjusted) == len(p_values)

        # Adjusted p-values should be >= original
        for orig, adj in zip(p_values, adjusted):
            assert adj >= orig

        # Should maintain order (if p[i] < p[j], then adj_p[i] <= adj_p[j])
        sorted_indices = np.argsort(p_values)
        sorted_adjusted = [adjusted[i] for i in sorted_indices]
        assert sorted_adjusted == sorted(sorted_adjusted)

    def test_benjamini_hochberg_known_values(self):
        """Test BH correction with known values."""
        # Example from literature
        p_values = [0.001, 0.008, 0.039, 0.041, 0.042]
        adjusted = apply_correction(p_values, method="benjamini-hochberg")

        # At alpha=0.05, first 3 should be significant
        assert adjusted[0] < 0.05
        assert adjusted[1] < 0.05
        assert adjusted[2] < 0.05

    def test_empty_p_values(self):
        """Test with empty list."""
        adjusted = apply_correction([], method="bonferroni")
        assert adjusted == []

    def test_single_p_value(self):
        """Test with single p-value."""
        adjusted = apply_correction([0.05], method="benjamini-hochberg")
        assert len(adjusted) == 1
        assert adjusted[0] == 0.05

    def test_all_zeros(self):
        """Test with all zero p-values."""
        p_values = [0.0, 0.0, 0.0]
        adjusted = apply_correction(p_values, method="benjamini-hochberg")

        assert all(adj == 0.0 for adj in adjusted)

    def test_all_ones(self):
        """Test with all p-values = 1.0."""
        p_values = [1.0, 1.0, 1.0]
        adjusted = apply_correction(p_values, method="bonferroni")

        assert all(adj == 1.0 for adj in adjusted)

    def test_invalid_method(self):
        """Test error handling for invalid method."""
        with pytest.raises(ValueError, match="Unknown correction method"):
            apply_correction([0.01], method="invalid")

    def test_n_total_bonferroni(self):
        """n_total makes Bonferroni use the given denominator."""
        p_values = [0.001, 0.01]
        # Default: n_tests = 2
        adj_default = apply_correction(p_values, method="bonferroni")
        # n_total = 20000 (simulating pre-filtered subset of 20k genes)
        adj_total = apply_correction(p_values, method="bonferroni", n_total=20000)
        assert adj_default == [0.002, 0.02]
        assert adj_total == [min(1.0, 0.001 * 20000), min(1.0, 0.01 * 20000)]

    def test_n_total_bh(self):
        """n_total makes BH correction more conservative."""
        p_values = [0.001, 0.01, 0.05]
        adj_default = apply_correction(p_values, method="benjamini-hochberg")
        adj_total = apply_correction(
            p_values, method="benjamini-hochberg", n_total=30000
        )
        # With larger n_total, adjusted p-values should be >= default
        for d, t in zip(adj_default, adj_total):
            assert t >= d - 1e-12

    def test_n_total_less_than_len_raises(self):
        """n_total < len(p_values) should raise ValueError."""
        with pytest.raises(ValueError, match="n_total"):
            apply_correction([0.01, 0.02, 0.03], method="bonferroni", n_total=2)


class TestFDRThreshold:
    """Tests for fdr_threshold function."""

    def test_fdr_threshold_basic(self):
        """Test basic FDR threshold calculation."""
        p_values = [0.001, 0.01, 0.05, 0.1, 0.2]
        threshold = fdr_threshold(p_values, alpha=0.05)

        # Threshold should be somewhere in the list
        assert threshold in p_values
        assert threshold > 0

    def test_fdr_threshold_no_significant(self):
        """Test when no p-values are significant."""
        p_values = [0.5, 0.6, 0.7, 0.8]
        threshold = fdr_threshold(p_values, alpha=0.05)

        # Should be 0 if none significant
        assert threshold == 0.0

    def test_fdr_threshold_all_significant(self):
        """Test when all p-values are significant."""
        p_values = [0.001, 0.002, 0.003, 0.004]
        threshold = fdr_threshold(p_values, alpha=0.05)

        # Threshold should be the largest p-value
        assert threshold == max(p_values)

    def test_fdr_threshold_empty(self):
        """Test with empty p-values."""
        threshold = fdr_threshold([], alpha=0.05)
        assert threshold == 0.0

    def test_fdr_threshold_alpha_variation(self):
        """Test that stricter alpha gives lower threshold."""
        p_values = [0.001, 0.01, 0.05, 0.1]

        threshold_05 = fdr_threshold(p_values, alpha=0.05)
        threshold_01 = fdr_threshold(p_values, alpha=0.01)

        # Stricter alpha should give lower or equal threshold
        assert threshold_01 <= threshold_05


class TestCorrectionProperties:
    """Tests for mathematical properties of corrections."""

    def test_monotonicity(self):
        """Test that correction preserves order."""
        p_values = [0.001, 0.005, 0.01, 0.05, 0.1]

        for method in ["bonferroni", "benjamini-hochberg"]:
            adjusted = apply_correction(p_values, method=method)

            # Get sort order
            orig_order = np.argsort(p_values)
            adj_order = np.argsort(adjusted)

            # Should have same order
            assert list(orig_order) == list(adj_order)

    def test_conservative_property(self):
        """Test that adjusted p-values are >= original."""
        p_values = [0.001, 0.01, 0.05, 0.1]

        for method in ["bonferroni", "benjamini-hochberg"]:
            adjusted = apply_correction(p_values, method=method)

            for orig, adj in zip(p_values, adjusted):
                assert adj >= orig, f"Adjusted {adj} < original {orig} for {method}"

    def test_bonferroni_more_conservative_than_bh(self):
        """Test that Bonferroni is more conservative than BH."""
        p_values = [0.001, 0.01, 0.05, 0.1]

        bonf = apply_correction(p_values, method="bonferroni")
        bh = apply_correction(p_values, method="benjamini-hochberg")

        for b, fdr in zip(bonf, bh):
            # Bonferroni should be >= BH (more conservative)
            assert b >= fdr or abs(b - fdr) < 1e-10
