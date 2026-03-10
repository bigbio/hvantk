"""
ROC analysis utilities for PSROC pipeline.

This module provides functions for computing ROC metrics, handling score
missingness, and finding optimal classification thresholds.
"""

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
import numpy as np

from sklearn.metrics import roc_curve, roc_auc_score


@dataclass
class ScoreMissingness:
    """Missingness statistics for a single prediction score.

    Attributes:
        score_name: Name of the prediction score (e.g., "CADD_phred").
        n_total: Total number of variants evaluated.
        n_present: Number of variants with non-null score values.
        n_missing: Number of variants with null/NaN score values.
        missingness_rate: Fraction of missing values (n_missing / n_total).
        included_in_analysis: Whether the score passed the missingness threshold.
        exclusion_reason: If excluded, the reason why (e.g., exceeded threshold).
    """

    score_name: str
    n_total: int
    n_present: int
    n_missing: int
    missingness_rate: float
    included_in_analysis: bool
    exclusion_reason: Optional[str] = None

    def to_dict(self) -> Dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "n_total": self.n_total,
            "n_present": self.n_present,
            "n_missing": self.n_missing,
            "missingness_rate": self.missingness_rate,
            "included_in_analysis": self.included_in_analysis,
            "exclusion_reason": self.exclusion_reason,
        }


@dataclass
class ROCResult:
    """Results from ROC analysis for a single prediction score.

    Attributes:
        score_name: Name of the prediction score.
        fpr: Array of false positive rates at each threshold.
        tpr: Array of true positive rates at each threshold.
        thresholds: Array of score thresholds corresponding to fpr/tpr.
        auc: Area under the ROC curve (0.0 to 1.0).
        optimal_threshold: Optimal classification threshold (e.g., Youden's J).
        sensitivity_at_optimal: Sensitivity (TPR) at the optimal threshold.
        specificity_at_optimal: Specificity (1-FPR) at the optimal threshold.
        n_variants_used: Number of variants with non-missing scores used in analysis.
        missingness: Missingness statistics for this score.
    """

    score_name: str
    fpr: np.ndarray
    tpr: np.ndarray
    thresholds: np.ndarray
    auc: float
    optimal_threshold: float
    sensitivity_at_optimal: float
    specificity_at_optimal: float
    n_variants_used: int
    missingness: ScoreMissingness

    def to_dict(self) -> Dict:
        """Convert to dictionary for JSON serialization (excludes arrays)."""
        return {
            "score_name": self.score_name,
            "auc": float(self.auc),
            "optimal_threshold": float(self.optimal_threshold),
            "sensitivity_at_optimal": float(self.sensitivity_at_optimal),
            "specificity_at_optimal": float(self.specificity_at_optimal),
            "n_variants_used": self.n_variants_used,
            "missingness": self.missingness.to_dict(),
        }


def compute_score_missingness(
    values: np.ndarray,
    score_name: str,
    max_missingness: float = 0.3,
) -> ScoreMissingness:
    """Compute missingness statistics for a single prediction score.

    Args:
        values: Array of score values (may contain NaN/None).
        score_name: Name of the score for reporting.
        max_missingness: Maximum allowed missingness rate for inclusion (default: 0.3).

    Returns:
        ScoreMissingness with computed statistics and inclusion decision.

    Raises:
        ValueError: If values array is empty.
    """
    if len(values) == 0:
        return ScoreMissingness(
            score_name=score_name,
            n_total=0,
            n_present=0,
            n_missing=0,
            missingness_rate=1.0,
            included_in_analysis=False,
            exclusion_reason="no variants available for this score",
        )

    n_total = len(values)
    n_missing = int(np.sum(np.isnan(values.astype(float))))
    n_present = n_total - n_missing
    missingness_rate = n_missing / n_total if n_total > 0 else 0.0

    # Scores with missingness > threshold are excluded
    included = missingness_rate <= max_missingness
    exclusion_reason = None
    if not included:
        exclusion_reason = (
            f"missingness_rate ({missingness_rate:.2f}) exceeds "
            f"max_missingness ({max_missingness:.2f})"
        )

    return ScoreMissingness(
        score_name=score_name,
        n_total=n_total,
        n_present=n_present,
        n_missing=n_missing,
        missingness_rate=missingness_rate,
        included_in_analysis=included,
        exclusion_reason=exclusion_reason,
    )


def filter_scores_by_missingness(
    missingness: Dict[str, ScoreMissingness],
    max_missingness: float = 0.3,
) -> Tuple[List[str], List[str]]:
    """Partition scores into included/excluded based on missingness threshold.

    Args:
        missingness: Dict mapping score name to ScoreMissingness.
        max_missingness: Maximum allowed missingness rate for inclusion.

    Returns:
        Tuple of (scores_included, scores_excluded) as lists of score names.
    """
    scores_included = []
    scores_excluded = []

    for score_name, stats in missingness.items():
        if stats.missingness_rate <= max_missingness:
            scores_included.append(score_name)
        else:
            scores_excluded.append(score_name)

    return scores_included, scores_excluded


def find_optimal_threshold(
    fpr: np.ndarray,
    tpr: np.ndarray,
    thresholds: np.ndarray,
    method: str = "youden",
) -> Tuple[float, float, float]:
    """Find optimal classification threshold using specified method.

    Args:
        fpr: Array of false positive rates.
        tpr: Array of true positive rates.
        thresholds: Array of score thresholds.
        method: Optimization method. One of:
            - "youden": Maximize Youden's J statistic (TPR - FPR)
            - "closest_to_corner": Minimize distance to (0, 1) corner
            - "f1": Maximize F1 score approximation

    Returns:
        Tuple of (optimal_threshold, sensitivity, specificity).

    Raises:
        ValueError: If method is not recognized or arrays are empty.
    """
    if len(fpr) == 0 or len(tpr) == 0 or len(thresholds) == 0:
        raise ValueError("Cannot find optimal threshold with empty arrays")

    if method == "youden":
        # Youden's J statistic: maximize (sensitivity + specificity - 1) = TPR - FPR
        j_scores = tpr - fpr
        optimal_idx = np.argmax(j_scores)

    elif method == "closest_to_corner":
        # Minimize Euclidean distance to (0, 1) - perfect classifier corner
        distances = np.sqrt((fpr - 0) ** 2 + (tpr - 1) ** 2)
        optimal_idx = np.argmin(distances)

    elif method == "f1":
        # Approximate F1 using precision = TP / (TP + FP) and recall = TPR
        # This is an approximation since we don't have actual counts
        # Precision ≈ TPR / (TPR + FPR) for balanced classes
        with np.errstate(divide="ignore", invalid="ignore"):
            precision = np.where(
                (tpr + fpr) > 0,
                tpr / (tpr + fpr),
                0,
            )
            f1_scores = np.where(
                (precision + tpr) > 0,
                2 * (precision * tpr) / (precision + tpr),
                0,
            )
        optimal_idx = np.argmax(f1_scores)

    else:
        raise ValueError(
            f"Unknown optimization method: {method}. "
            f"Expected one of: 'youden', 'closest_to_corner', 'f1'"
        )

    optimal_threshold = thresholds[optimal_idx]
    sensitivity = tpr[optimal_idx]
    specificity = 1 - fpr[optimal_idx]

    return float(optimal_threshold), float(sensitivity), float(specificity)


def compute_roc_metrics(
    labels: np.ndarray,
    scores: Dict[str, np.ndarray],
    max_missingness: float = 0.3,
    pos_label: int = 1,
    threshold_method: str = "youden",
) -> Dict[str, ROCResult]:
    """Compute ROC metrics for multiple prediction scores.

    For each score, variants with missing values are excluded from that score's
    ROC analysis. Scores with missingness exceeding max_missingness are excluded
    entirely from the results.

    Args:
        labels: Binary labels (1=pathogenic, 0=benign).
        scores: Dict mapping score name to array of score values.
        max_missingness: Maximum allowed missingness rate per score (default: 0.3).
        pos_label: Label value considered positive (default: 1).
        threshold_method: Method for finding optimal threshold (default: "youden").

    Returns:
        Dict mapping score name to ROCResult. Only scores that pass the
        missingness threshold are included.

    Raises:
        ValueError: If labels and score arrays have different lengths,
                   if no valid labels exist, or if all scores are excluded.
    """
    if len(labels) == 0:
        raise ValueError("Labels array is empty")

    results = {}

    for score_name, score_values in scores.items():
        if len(score_values) != len(labels):
            raise ValueError(
                f"Score '{score_name}' has {len(score_values)} values but "
                f"labels has {len(labels)} values"
            )

        # Compute missingness
        missingness = compute_score_missingness(
            score_values, score_name, max_missingness
        )

        # Skip scores with too much missingness
        if not missingness.included_in_analysis:
            continue

        # Filter to non-missing values
        score_array = np.array(score_values, dtype=float)
        valid_mask = ~np.isnan(score_array)
        valid_labels = labels[valid_mask]
        valid_scores = score_array[valid_mask]

        n_variants_used = len(valid_labels)

        # Check for valid class distribution
        unique_labels = np.unique(valid_labels)
        if len(unique_labels) < 2:
            raise ValueError(
                f"Score '{score_name}' has only one class after filtering "
                f"missing values. Need both pathogenic and benign variants."
            )

        if pos_label not in unique_labels:
            raise ValueError(
                f"Positive label {pos_label} not found in labels for "
                f"score '{score_name}'"
            )

        # Compute ROC curve
        fpr, tpr, thresholds = roc_curve(
            valid_labels, valid_scores, pos_label=pos_label
        )

        # Compute AUC - convert labels to binary (0/1) based on pos_label
        # roc_auc_score doesn't support pos_label parameter, so we need to
        # ensure labels are in binary format where 1 represents the positive class
        binary_labels = (valid_labels == pos_label).astype(int)
        auc = roc_auc_score(binary_labels, valid_scores)

        # Find optimal threshold
        optimal_threshold, sensitivity, specificity = find_optimal_threshold(
            fpr, tpr, thresholds, method=threshold_method
        )

        results[score_name] = ROCResult(
            score_name=score_name,
            fpr=fpr,
            tpr=tpr,
            thresholds=thresholds,
            auc=auc,
            optimal_threshold=optimal_threshold,
            sensitivity_at_optimal=sensitivity,
            specificity_at_optimal=specificity,
            n_variants_used=n_variants_used,
            missingness=missingness,
        )

    if len(results) == 0 and len(scores) > 0:
        import logging as _logging

        excluded_scores = list(scores.keys())
        _logging.getLogger(__name__).warning(
            "All scores were excluded due to high missingness: %s. "
            "Consider increasing --max-missingness threshold.",
            excluded_scores,
        )

    return results


def compute_all_missingness(
    scores: Dict[str, np.ndarray],
    max_missingness: float = 0.3,
) -> Dict[str, ScoreMissingness]:
    """Compute missingness statistics for all scores.

    Unlike compute_roc_metrics, this function returns missingness stats for
    ALL scores, including those that would be excluded from analysis.

    Args:
        scores: Dict mapping score name to array of score values.
        max_missingness: Maximum allowed missingness rate for inclusion.

    Returns:
        Dict mapping score name to ScoreMissingness for all input scores.
    """
    return {
        score_name: compute_score_missingness(values, score_name, max_missingness)
        for score_name, values in scores.items()
    }
