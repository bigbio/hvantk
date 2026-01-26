"""
Tests for PSROC ROC analysis functions.

This module tests the core ROC computation functions including AUC calculation,
optimal threshold finding, and handling of edge cases.
"""

import numpy as np
import pytest

from hvantk.psroc.roc import (
    ROCResult,
    ScoreMissingness,
    compute_roc_metrics,
    find_optimal_threshold,
)


class TestROCResult:
    """Test ROCResult dataclass."""

    def test_roc_result_creation(self):
        """Test creating a ROCResult instance."""
        missingness = ScoreMissingness(
            score_name="test_score",
            n_total=100,
            n_present=95,
            n_missing=5,
            missingness_rate=0.05,
            included_in_analysis=True,
        )

        result = ROCResult(
            score_name="test_score",
            fpr=np.array([0.0, 0.5, 1.0]),
            tpr=np.array([0.0, 0.8, 1.0]),
            thresholds=np.array([1.0, 0.5, 0.0]),
            auc=0.85,
            optimal_threshold=0.5,
            sensitivity_at_optimal=0.8,
            specificity_at_optimal=0.5,
            n_variants_used=95,
            missingness=missingness,
        )

        assert result.score_name == "test_score"
        assert result.auc == 0.85
        assert result.n_variants_used == 95

    def test_roc_result_to_dict(self):
        """Test ROCResult serialization to dict."""
        missingness = ScoreMissingness(
            score_name="test_score",
            n_total=100,
            n_present=95,
            n_missing=5,
            missingness_rate=0.05,
            included_in_analysis=True,
        )

        result = ROCResult(
            score_name="test_score",
            fpr=np.array([0.0, 0.5, 1.0]),
            tpr=np.array([0.0, 0.8, 1.0]),
            thresholds=np.array([1.0, 0.5, 0.0]),
            auc=0.85,
            optimal_threshold=0.5,
            sensitivity_at_optimal=0.8,
            specificity_at_optimal=0.5,
            n_variants_used=95,
            missingness=missingness,
        )

        d = result.to_dict()
        assert d["score_name"] == "test_score"
        assert d["auc"] == 0.85
        assert "fpr" not in d  # Arrays excluded from dict
        assert "missingness" in d
        assert d["missingness"]["n_total"] == 100


class TestFindOptimalThreshold:
    """Test optimal threshold finding functions."""

    def test_youden_method(self):
        """Test Youden's J statistic method."""
        # Simple case: threshold at index 1 maximizes TPR - FPR
        fpr = np.array([0.0, 0.2, 0.5, 1.0])
        tpr = np.array([0.0, 0.9, 0.95, 1.0])
        thresholds = np.array([1.0, 0.7, 0.4, 0.0])

        # Youden's J = TPR - FPR
        # Index 0: 0.0 - 0.0 = 0.0
        # Index 1: 0.9 - 0.2 = 0.7  <- max
        # Index 2: 0.95 - 0.5 = 0.45
        # Index 3: 1.0 - 1.0 = 0.0

        threshold, sens, spec = find_optimal_threshold(
            fpr, tpr, thresholds, method="youden"
        )

        assert threshold == 0.7
        assert sens == 0.9
        assert spec == pytest.approx(0.8)  # 1 - 0.2

    def test_closest_to_corner_method(self):
        """Test closest to corner (0,1) method."""
        # Point closest to (0, 1) should be selected
        fpr = np.array([0.0, 0.1, 0.5, 1.0])
        tpr = np.array([0.0, 0.95, 0.98, 1.0])
        thresholds = np.array([1.0, 0.6, 0.3, 0.0])

        # Distance to (0, 1):
        # Index 0: sqrt(0^2 + 1^2) = 1.0
        # Index 1: sqrt(0.1^2 + 0.05^2) ≈ 0.112  <- min
        # Index 2: sqrt(0.5^2 + 0.02^2) ≈ 0.5
        # Index 3: sqrt(1^2 + 0^2) = 1.0

        threshold, sens, spec = find_optimal_threshold(
            fpr, tpr, thresholds, method="closest_to_corner"
        )

        assert threshold == 0.6
        assert sens == 0.95
        assert spec == pytest.approx(0.9)

    def test_invalid_method_raises_error(self):
        """Test that invalid method raises ValueError."""
        fpr = np.array([0.0, 0.5, 1.0])
        tpr = np.array([0.0, 0.8, 1.0])
        thresholds = np.array([1.0, 0.5, 0.0])

        with pytest.raises(ValueError, match="Unknown optimization method"):
            find_optimal_threshold(fpr, tpr, thresholds, method="invalid")

    def test_empty_arrays_raise_error(self):
        """Test that empty arrays raise ValueError."""
        with pytest.raises(ValueError, match="empty arrays"):
            find_optimal_threshold(
                np.array([]), np.array([]), np.array([]), method="youden"
            )


class TestComputeROCMetrics:
    """Test ROC metrics computation."""

    def test_perfect_separator(self):
        """Test AUC=1.0 for perfect separation."""
        # Perfect classifier: all pathogenic have higher scores than benign
        labels = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1])
        scores = {
            "perfect_score": np.array(
                [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
            )
        }

        results = compute_roc_metrics(labels, scores)

        assert "perfect_score" in results
        assert results["perfect_score"].auc == pytest.approx(1.0)

    def test_random_classifier(self):
        """Test AUC≈0.5 for random/uninformative scores."""
        np.random.seed(42)
        n = 1000
        labels = np.array([0] * (n // 2) + [1] * (n // 2))
        scores = {"random_score": np.random.rand(n)}

        results = compute_roc_metrics(labels, scores)

        # Random classifier should have AUC close to 0.5
        assert "random_score" in results
        assert 0.4 < results["random_score"].auc < 0.6

    def test_inverse_classifier(self):
        """Test AUC≈0.0 for inverse classifier (benign scores higher)."""
        # Inverse classifier: all benign have higher scores than pathogenic
        labels = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1])
        scores = {
            "inverse_score": np.array(
                [0.6, 0.7, 0.8, 0.9, 1.0, 0.1, 0.2, 0.3, 0.4, 0.5]
            )
        }

        results = compute_roc_metrics(labels, scores)

        assert "inverse_score" in results
        assert results["inverse_score"].auc == pytest.approx(0.0)

    def test_multiple_scores(self):
        """Test computing metrics for multiple scores simultaneously."""
        labels = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1])
        scores = {
            "good_score": np.array(
                [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
            ),
            "medium_score": np.array(
                [0.2, 0.3, 0.5, 0.6, 0.7, 0.4, 0.5, 0.8, 0.9, 0.95]
            ),
        }

        results = compute_roc_metrics(labels, scores)

        assert len(results) == 2
        assert "good_score" in results
        assert "medium_score" in results
        assert results["good_score"].auc > results["medium_score"].auc

    def test_missing_values_excluded(self):
        """Test that NaN score values are excluded from ROC computation."""
        labels = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1])
        scores = {
            "score_with_nans": np.array(
                [0.1, np.nan, 0.3, 0.4, 0.5, 0.6, 0.7, np.nan, 0.9, 1.0]
            )
        }

        results = compute_roc_metrics(labels, scores)

        assert "score_with_nans" in results
        assert results["score_with_nans"].n_variants_used == 8
        assert results["score_with_nans"].missingness.n_missing == 2

    def test_high_missingness_excludes_score(self):
        """Test that scores with high missingness are excluded."""
        labels = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1])
        scores = {
            "good_score": np.array(
                [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
            ),
            "bad_score": np.array(
                [np.nan, np.nan, np.nan, np.nan, 0.5, np.nan, np.nan, np.nan, np.nan, 1.0]
            ),
        }

        results = compute_roc_metrics(labels, scores, max_missingness=0.3)

        # good_score should be included (0% missing)
        assert "good_score" in results
        # bad_score should be excluded (80% missing > 30% threshold)
        assert "bad_score" not in results

    def test_empty_labels_raises_error(self):
        """Test that empty labels array raises ValueError."""
        with pytest.raises(ValueError, match="Labels array is empty"):
            compute_roc_metrics(np.array([]), {"score": np.array([])})

    def test_mismatched_lengths_raises_error(self):
        """Test that mismatched array lengths raise ValueError."""
        labels = np.array([0, 0, 1, 1])
        scores = {"score": np.array([0.1, 0.2, 0.3])}  # Wrong length

        with pytest.raises(ValueError, match="has 3 values but labels has 4"):
            compute_roc_metrics(labels, scores)

    def test_single_class_raises_error(self):
        """Test that single class after filtering raises ValueError."""
        labels = np.array([0, 0, 0, 0, 1, 1, 1, 1])
        # NaN for all pathogenic variants
        scores = {
            "problematic": np.array(
                [0.1, 0.2, 0.3, 0.4, np.nan, np.nan, np.nan, np.nan]
            )
        }

        with pytest.raises(ValueError, match="only one class"):
            compute_roc_metrics(labels, scores, max_missingness=1.0)

    def test_all_scores_excluded_raises_error(self):
        """Test error when all scores exceed missingness threshold."""
        labels = np.array([0, 0, 1, 1])
        scores = {
            "score1": np.array([np.nan, np.nan, np.nan, 0.5]),  # 75% missing
            "score2": np.array([np.nan, np.nan, 0.5, np.nan]),  # 75% missing
        }

        with pytest.raises(ValueError, match="All scores were excluded"):
            compute_roc_metrics(labels, scores, max_missingness=0.3)

    def test_optimal_threshold_included(self):
        """Test that optimal threshold is computed and included."""
        labels = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1])
        scores = {
            "test_score": np.array(
                [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
            )
        }

        results = compute_roc_metrics(labels, scores, threshold_method="youden")

        result = results["test_score"]
        assert result.optimal_threshold is not None
        assert 0.0 <= result.sensitivity_at_optimal <= 1.0
        assert 0.0 <= result.specificity_at_optimal <= 1.0

    def test_custom_pos_label(self):
        """Test using custom positive label."""
        # Use 2 as positive label instead of 1
        labels = np.array([0, 0, 0, 0, 0, 2, 2, 2, 2, 2])
        scores = {
            "score": np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
        }

        results = compute_roc_metrics(labels, scores, pos_label=2)

        assert "score" in results
        assert results["score"].auc == pytest.approx(1.0)


class TestScoreMissingness:
    """Test ScoreMissingness dataclass."""

    def test_missingness_to_dict(self):
        """Test ScoreMissingness serialization."""
        missingness = ScoreMissingness(
            score_name="test",
            n_total=100,
            n_present=70,
            n_missing=30,
            missingness_rate=0.3,
            included_in_analysis=False,
            exclusion_reason="exceeded threshold",
        )

        d = missingness.to_dict()
        assert d["n_total"] == 100
        assert d["n_present"] == 70
        assert d["n_missing"] == 30
        assert d["missingness_rate"] == 0.3
        assert d["included_in_analysis"] is False
        assert d["exclusion_reason"] == "exceeded threshold"


if __name__ == "__main__":
    print("Running PSROC ROC tests...")

    print("\n1. Testing ROCResult dataclass")
    test_result = TestROCResult()
    test_result.test_roc_result_creation()
    print("  ✓ ROCResult creation works")
    test_result.test_roc_result_to_dict()
    print("  ✓ ROCResult to_dict works")

    print("\n2. Testing optimal threshold finding")
    test_threshold = TestFindOptimalThreshold()
    test_threshold.test_youden_method()
    print("  ✓ Youden method works")
    test_threshold.test_closest_to_corner_method()
    print("  ✓ Closest to corner method works")

    print("\n3. Testing ROC metrics computation")
    test_roc = TestComputeROCMetrics()
    test_roc.test_perfect_separator()
    print("  ✓ Perfect separator has AUC=1.0")
    test_roc.test_random_classifier()
    print("  ✓ Random classifier has AUC≈0.5")
    test_roc.test_multiple_scores()
    print("  ✓ Multiple scores work")
    test_roc.test_missing_values_excluded()
    print("  ✓ Missing values excluded correctly")

    print("\n✅ All tests passed!")
