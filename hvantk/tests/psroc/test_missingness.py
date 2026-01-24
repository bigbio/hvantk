"""
Tests for PSROC score missingness handling.

This module tests the missingness computation, filtering, and reporting
functionality of the PSROC module.
"""

import numpy as np
import pytest

from hvantk.psroc.roc import (
    ScoreMissingness,
    compute_score_missingness,
    filter_scores_by_missingness,
    compute_all_missingness,
)


class TestComputeScoreMissingness:
    """Test compute_score_missingness function."""

    def test_no_missing_values(self):
        """Test with no missing values."""
        values = np.array([0.1, 0.2, 0.3, 0.4, 0.5])

        result = compute_score_missingness(values, "test_score")

        assert result.score_name == "test_score"
        assert result.n_total == 5
        assert result.n_present == 5
        assert result.n_missing == 0
        assert result.missingness_rate == 0.0
        assert result.included_in_analysis is True
        assert result.exclusion_reason is None

    def test_all_missing_values(self):
        """Test with all missing values."""
        values = np.array([np.nan, np.nan, np.nan, np.nan])

        result = compute_score_missingness(values, "test_score", max_missingness=0.3)

        assert result.n_total == 4
        assert result.n_present == 0
        assert result.n_missing == 4
        assert result.missingness_rate == 1.0
        assert result.included_in_analysis is False
        assert "exceeds" in result.exclusion_reason

    def test_partial_missing_below_threshold(self):
        """Test with partial missing values below threshold."""
        values = np.array([0.1, 0.2, np.nan, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

        result = compute_score_missingness(values, "test_score", max_missingness=0.3)

        assert result.n_total == 10
        assert result.n_present == 9
        assert result.n_missing == 1
        assert result.missingness_rate == pytest.approx(0.1)
        assert result.included_in_analysis is True

    def test_partial_missing_above_threshold(self):
        """Test with partial missing values above threshold."""
        values = np.array([0.1, np.nan, np.nan, np.nan, np.nan, 0.6, 0.7, 0.8, 0.9, 1.0])

        result = compute_score_missingness(values, "test_score", max_missingness=0.3)

        assert result.n_total == 10
        assert result.n_present == 6
        assert result.n_missing == 4
        assert result.missingness_rate == pytest.approx(0.4)
        assert result.included_in_analysis is False
        assert "0.40" in result.exclusion_reason
        assert "0.30" in result.exclusion_reason

    def test_exactly_at_threshold(self):
        """Test score exactly at missingness threshold is included."""
        # 30% missing with threshold of 0.3 should be included (<=)
        values = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, np.nan, np.nan, np.nan])

        result = compute_score_missingness(values, "test_score", max_missingness=0.3)

        assert result.missingness_rate == pytest.approx(0.3)
        assert result.included_in_analysis is True

    def test_just_above_threshold(self):
        """Test score just above missingness threshold is excluded."""
        # 31% missing with threshold of 0.3 should be excluded (>)
        values = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, np.nan, np.nan, np.nan, np.nan])

        result = compute_score_missingness(values, "test_score", max_missingness=0.3)

        assert result.missingness_rate == pytest.approx(0.4)
        assert result.included_in_analysis is False

    def test_empty_array_raises_error(self):
        """Test that empty array raises ValueError."""
        with pytest.raises(ValueError, match="Cannot compute missingness for empty array"):
            compute_score_missingness(np.array([]), "test_score")

    def test_custom_threshold(self):
        """Test with custom missingness threshold."""
        values = np.array([0.1, 0.2, np.nan, np.nan, np.nan])  # 60% missing

        # With 50% threshold, should be excluded
        result1 = compute_score_missingness(values, "test", max_missingness=0.5)
        assert result1.included_in_analysis is False

        # With 70% threshold, should be included
        result2 = compute_score_missingness(values, "test", max_missingness=0.7)
        assert result2.included_in_analysis is True


class TestFilterScoresByMissingness:
    """Test filter_scores_by_missingness function."""

    def test_all_included(self):
        """Test when all scores pass threshold."""
        missingness = {
            "score_a": ScoreMissingness(
                score_name="score_a",
                n_total=100,
                n_present=95,
                n_missing=5,
                missingness_rate=0.05,
                included_in_analysis=True,
            ),
            "score_b": ScoreMissingness(
                score_name="score_b",
                n_total=100,
                n_present=80,
                n_missing=20,
                missingness_rate=0.20,
                included_in_analysis=True,
            ),
        }

        included, excluded = filter_scores_by_missingness(missingness, max_missingness=0.3)

        assert included == ["score_a", "score_b"]
        assert excluded == []

    def test_all_excluded(self):
        """Test when all scores fail threshold."""
        missingness = {
            "score_a": ScoreMissingness(
                score_name="score_a",
                n_total=100,
                n_present=50,
                n_missing=50,
                missingness_rate=0.50,
                included_in_analysis=False,
            ),
            "score_b": ScoreMissingness(
                score_name="score_b",
                n_total=100,
                n_present=30,
                n_missing=70,
                missingness_rate=0.70,
                included_in_analysis=False,
            ),
        }

        included, excluded = filter_scores_by_missingness(missingness, max_missingness=0.3)

        assert included == []
        assert sorted(excluded) == ["score_a", "score_b"]

    def test_mixed_inclusion(self):
        """Test with mix of included and excluded scores."""
        missingness = {
            "good_score": ScoreMissingness(
                score_name="good_score",
                n_total=100,
                n_present=95,
                n_missing=5,
                missingness_rate=0.05,
                included_in_analysis=True,
            ),
            "bad_score": ScoreMissingness(
                score_name="bad_score",
                n_total=100,
                n_present=50,
                n_missing=50,
                missingness_rate=0.50,
                included_in_analysis=False,
            ),
            "borderline_score": ScoreMissingness(
                score_name="borderline_score",
                n_total=100,
                n_present=70,
                n_missing=30,
                missingness_rate=0.30,
                included_in_analysis=True,
            ),
        }

        included, excluded = filter_scores_by_missingness(missingness, max_missingness=0.3)

        assert sorted(included) == ["borderline_score", "good_score"]
        assert excluded == ["bad_score"]

    def test_empty_input(self):
        """Test with empty missingness dict."""
        included, excluded = filter_scores_by_missingness({}, max_missingness=0.3)

        assert included == []
        assert excluded == []


class TestComputeAllMissingness:
    """Test compute_all_missingness function."""

    def test_multiple_scores(self):
        """Test computing missingness for multiple scores."""
        scores = {
            "score_a": np.array([0.1, 0.2, 0.3, 0.4, 0.5]),
            "score_b": np.array([0.1, np.nan, 0.3, np.nan, 0.5]),
            "score_c": np.array([np.nan, np.nan, np.nan, np.nan, np.nan]),
        }

        result = compute_all_missingness(scores, max_missingness=0.3)

        assert len(result) == 3
        assert "score_a" in result
        assert "score_b" in result
        assert "score_c" in result

        # score_a: 0% missing, included
        assert result["score_a"].missingness_rate == 0.0
        assert result["score_a"].included_in_analysis is True

        # score_b: 40% missing, excluded
        assert result["score_b"].missingness_rate == pytest.approx(0.4)
        assert result["score_b"].included_in_analysis is False

        # score_c: 100% missing, excluded
        assert result["score_c"].missingness_rate == 1.0
        assert result["score_c"].included_in_analysis is False

    def test_empty_scores_dict(self):
        """Test with empty scores dict."""
        result = compute_all_missingness({}, max_missingness=0.3)
        assert result == {}


class TestMissingnessReportFormat:
    """Test that missingness data can be serialized to expected format."""

    def test_missingness_report_json_structure(self):
        """Test that missingness report matches expected JSON schema."""
        missingness = ScoreMissingness(
            score_name="CADD_phred",
            n_total=1250,
            n_present=1200,
            n_missing=50,
            missingness_rate=0.04,
            included_in_analysis=True,
        )

        report = missingness.to_dict()

        # Verify all expected fields are present
        assert "n_present" in report
        assert "n_missing" in report
        assert "missingness_rate" in report
        assert "included_in_analysis" in report

        # Verify values
        assert report["n_present"] == 1200
        assert report["n_missing"] == 50
        assert report["missingness_rate"] == 0.04
        assert report["included_in_analysis"] is True

    def test_excluded_score_report_includes_reason(self):
        """Test that excluded scores include exclusion reason."""
        missingness = ScoreMissingness(
            score_name="REVEL_score",
            n_total=1250,
            n_present=800,
            n_missing=450,
            missingness_rate=0.36,
            included_in_analysis=False,
            exclusion_reason="missingness_rate (0.36) exceeds max_missingness (0.30)",
        )

        report = missingness.to_dict()

        assert report["included_in_analysis"] is False
        assert report["exclusion_reason"] is not None
        assert "0.36" in report["exclusion_reason"]


class TestClassStratifiedMissingness:
    """Test missingness computation stratified by class label."""

    def test_class_stratified_missingness_concept(self):
        """Test computing missingness separately for pathogenic and benign.

        This test demonstrates the pattern for stratified missingness,
        which will be used in the pipeline for bias detection.
        """
        # Create data with different missingness patterns by class
        labels = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1])
        scores = np.array([0.1, 0.2, np.nan, 0.4, 0.5, np.nan, np.nan, np.nan, 0.9, 1.0])

        # Compute stratified missingness
        benign_mask = labels == 0
        pathogenic_mask = labels == 1

        benign_scores = scores[benign_mask]
        pathogenic_scores = scores[pathogenic_mask]

        benign_missingness = compute_score_missingness(
            benign_scores, "score_benign", max_missingness=1.0
        )
        pathogenic_missingness = compute_score_missingness(
            pathogenic_scores, "score_pathogenic", max_missingness=1.0
        )

        # Benign: 1/5 = 20% missing
        assert benign_missingness.n_missing == 1
        assert benign_missingness.missingness_rate == pytest.approx(0.2)

        # Pathogenic: 3/5 = 60% missing
        assert pathogenic_missingness.n_missing == 3
        assert pathogenic_missingness.missingness_rate == pytest.approx(0.6)

        # This pattern would be flagged as potentially biased
        # (different missingness rates between classes)


if __name__ == "__main__":
    print("Running PSROC Missingness tests...")

    print("\n1. Testing compute_score_missingness")
    test_compute = TestComputeScoreMissingness()
    test_compute.test_no_missing_values()
    print("  ✓ No missing values handled correctly")
    test_compute.test_all_missing_values()
    print("  ✓ All missing values handled correctly")
    test_compute.test_partial_missing_below_threshold()
    print("  ✓ Partial missing below threshold works")
    test_compute.test_partial_missing_above_threshold()
    print("  ✓ Partial missing above threshold works")
    test_compute.test_exactly_at_threshold()
    print("  ✓ Exactly at threshold is included")
    test_compute.test_just_above_threshold()
    print("  ✓ Just above threshold is excluded")

    print("\n2. Testing filter_scores_by_missingness")
    test_filter = TestFilterScoresByMissingness()
    test_filter.test_all_included()
    print("  ✓ All included case works")
    test_filter.test_all_excluded()
    print("  ✓ All excluded case works")
    test_filter.test_mixed_inclusion()
    print("  ✓ Mixed inclusion works")

    print("\n3. Testing compute_all_missingness")
    test_all = TestComputeAllMissingness()
    test_all.test_multiple_scores()
    print("  ✓ Multiple scores work")

    print("\n4. Testing report format")
    test_report = TestMissingnessReportFormat()
    test_report.test_missingness_report_json_structure()
    print("  ✓ Report JSON structure correct")
    test_report.test_excluded_score_report_includes_reason()
    print("  ✓ Exclusion reason included")

    print("\n5. Testing class-stratified missingness")
    test_stratified = TestClassStratifiedMissingness()
    test_stratified.test_class_stratified_missingness_concept()
    print("  ✓ Stratified missingness works")

    print("\n✅ All missingness tests passed!")
