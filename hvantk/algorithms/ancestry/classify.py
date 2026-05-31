"""Ancestry classification using Random Forest.

This module provides functions for training a Random Forest classifier
on PC scores from a reference panel and predicting ancestry for query
samples.
"""

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
# Make scikit-learn import optional - only required for classifier training.
try:
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.metrics import classification_report, confusion_matrix
    from sklearn.model_selection import StratifiedKFold, cross_val_predict

    SKLEARN_AVAILABLE = True
except ImportError:  # pragma: no cover
    SKLEARN_AVAILABLE = False
    # Defined to satisfy linters; guarded by SKLEARN_AVAILABLE.
    RandomForestClassifier = None
    classification_report = None
    confusion_matrix = None
    StratifiedKFold = None
    cross_val_predict = None

from hvantk.algorithms.ancestry.constants import (
    ANCESTRY_PROB_COL,
    DEFAULT_MIN_PROB,
    DEFAULT_N_CV_FOLDS,
    DEFAULT_N_ESTIMATORS,
    DEFAULT_N_PCS_CLASSIFY,
    DEFAULT_RANDOM_SEED,
    MIN_SAMPLES_PER_POP,
    PREDICTED_ANCESTRY_COL,
)

logger = logging.getLogger(__name__)


@dataclass
class ClassificationResult:
    """Container for classification results.

    Attributes
    ----------
    predictions_df : pd.DataFrame
        DataFrame with predictions for all samples.
    model : RandomForestClassifier
        Trained classifier model.
    classes : List[str]
        List of class labels (ancestry populations).
    validation_metrics : Dict[str, Any], optional
        Cross-validation metrics if validation was performed.
    confusion_matrix : np.ndarray, optional
        Confusion matrix from validation.
    confusion_matrix_labels : Tuple[np.ndarray, np.ndarray], optional
        Tuple of (y_true, y_pred) arrays used to compute confusion matrix.
        Stored for visualization purposes.
    feature_importances : Dict[str, float], optional
        Importance scores for each PC feature.
    """

    predictions_df: pd.DataFrame
    model: RandomForestClassifier
    classes: List[str]
    validation_metrics: Optional[Dict[str, Any]] = None
    confusion_matrix: Optional[np.ndarray] = field(default=None, repr=False)
    confusion_matrix_labels: Optional[Tuple[np.ndarray, np.ndarray]] = field(
        default=None, repr=False
    )
    feature_importances: Optional[Dict[str, float]] = None

    def get_accuracy(self) -> Optional[float]:
        """Get overall accuracy from validation metrics.

        Returns
        -------
        float or None
            Overall accuracy if validation was performed.

        """
        if self.validation_metrics and "accuracy" in self.validation_metrics:
            return self.validation_metrics["accuracy"]
        return None


def _get_pc_columns(n_pcs: int) -> List[str]:
    """Generate list of PC column names.

    Parameters
    ----------
    n_pcs : int
        Number of PCs.

    Returns
    -------
    List[str]
        Column names ['PC1', 'PC2', ..., 'PCn'].
    """
    return [f"PC{i + 1}" for i in range(n_pcs)]


def validate_training_data(
    labels: pd.Series,
    min_samples_per_pop: int = MIN_SAMPLES_PER_POP,
) -> Dict[str, int]:
    """Validate that training data has sufficient samples per population.

    Parameters
    ----------
    labels : pd.Series
        Ancestry labels for training samples.
    min_samples_per_pop : int, optional
        Minimum required samples per population. Default: 10.

    Returns
    -------
    Dict[str, int]
        Dictionary of population -> sample count.

    Raises
    ------
    ValueError
        If any population has fewer than min_samples_per_pop samples.
        If fewer than 2 populations are present.
    """
    pop_counts = labels.value_counts().to_dict()

    n_pops = len(pop_counts)
    if n_pops < 2:
        raise ValueError(
            f"At least 2 populations required for classification, "
            f"but only {n_pops} found: {list(pop_counts.keys())}"
        )

    insufficient = {
        pop: count for pop, count in pop_counts.items() if count < min_samples_per_pop
    }

    if insufficient:
        msg = "; ".join(
            f"'{pop}' has {count} samples" for pop, count in insufficient.items()
        )
        raise ValueError(
            f"Populations with insufficient samples (need at least "
            f"{min_samples_per_pop}): {msg}"
        )

    logger.info(f"Training data validated: {n_pops} populations")
    for pop, count in sorted(pop_counts.items()):
        logger.debug(f"  {pop}: {count} samples")

    return pop_counts


def train_classifier(
    scores_df: pd.DataFrame,
    labels: pd.Series,
    n_pcs: int = DEFAULT_N_PCS_CLASSIFY,
    n_estimators: int = DEFAULT_N_ESTIMATORS,
    random_state: int = DEFAULT_RANDOM_SEED,
    validate: bool = True,
    n_cv_folds: int = DEFAULT_N_CV_FOLDS,
    min_samples_per_pop: int = MIN_SAMPLES_PER_POP,
) -> ClassificationResult:
    """Train Random Forest classifier on labeled samples.

    Parameters
    ----------
    scores_df : pd.DataFrame
        PC scores with columns PC1, PC2, ..., PCn.
    labels : pd.Series
        Ancestry labels aligned with scores_df index.
    n_pcs : int, optional
        Number of PCs to use as features. Default: 10.
    n_estimators : int, optional
        Number of trees in the forest. Default: 100.
    random_state : int, optional
        Random seed for reproducibility. Default: 42.
    validate : bool, optional
        Whether to perform cross-validation. Default: True.
    n_cv_folds : int, optional
        Number of cross-validation folds. Default: 5.
    min_samples_per_pop : int, optional
        Minimum samples per population. Default: 10.

    Returns
    -------
    ClassificationResult
        Trained model with optional validation metrics.

    Raises
    ------
    ValueError
        If training data validation fails or if required columns missing.

    Example
    -------
    >>> result = train_classifier(
    ...     scores_df=reference_scores,
    ...     labels=reference_labels,
    ...     n_pcs=10,
    ... )
    >>> print(f"Accuracy: {result.get_accuracy():.2%}")
    """
    if not SKLEARN_AVAILABLE:
        raise RuntimeError(
            "Ancestry classification requires scikit-learn. Install it with "
            "'poetry install --extras ml' (or --extras ancestry / --extras psroc)."
        )

    # Get PC columns
    pc_cols = _get_pc_columns(n_pcs)

    # Validate PC columns exist
    missing_cols = set(pc_cols) - set(scores_df.columns)
    if missing_cols:
        raise ValueError(
            f"Missing PC columns in scores_df: {sorted(missing_cols)}. "
            f"Available columns: {list(scores_df.columns)}"
        )

    # Validate training data
    pop_counts = validate_training_data(labels, min_samples_per_pop)

    # Prepare features
    X_train = scores_df[pc_cols].values
    y_train = labels.values

    logger.info(
        f"Training Random Forest with {n_estimators} trees "
        f"on {len(y_train)} samples using {n_pcs} PCs"
    )

    # Create classifier with class balancing
    rf = RandomForestClassifier(
        n_estimators=n_estimators,
        max_features="sqrt",
        min_samples_leaf=5,
        class_weight="balanced",
        random_state=random_state,
        n_jobs=-1,
    )

    # Perform cross-validation if requested
    validation_metrics = None
    conf_matrix = None
    conf_matrix_labels = None

    if validate:
        # Verify n_cv_folds is valid for the smallest class
        min_count = min(pop_counts.values())
        if n_cv_folds > min_count:
            raise ValueError(
                f"n_cv_folds ({n_cv_folds}) exceeds the smallest class count ({min_count}). "
                f"StratifiedKFold requires n_splits <= min_class_count. "
                f"Population counts: {pop_counts}"
            )

        logger.info(f"Performing {n_cv_folds}-fold stratified cross-validation")

        cv = StratifiedKFold(
            n_splits=n_cv_folds,
            shuffle=True,
            random_state=random_state,
        )

        # Get cross-validated predictions
        y_pred_cv = cross_val_predict(rf, X_train, y_train, cv=cv)

        # Compute metrics
        validation_metrics = classification_report(
            y_train, y_pred_cv, output_dict=True, zero_division=0
        )
        conf_matrix = confusion_matrix(y_train, y_pred_cv)
        # Store labels for visualization
        conf_matrix_labels = (y_train, y_pred_cv)

        accuracy = validation_metrics.get("accuracy", 0)
        logger.info(f"Cross-validation accuracy: {accuracy:.2%}")

        # Log per-class metrics
        classes = sorted(set(y_train))
        for cls in classes:
            if cls in validation_metrics:
                metrics = validation_metrics[cls]
                logger.debug(
                    f"  {cls}: precision={metrics['precision']:.2f}, "
                    f"recall={metrics['recall']:.2f}, "
                    f"f1={metrics['f1-score']:.2f}"
                )

    # Train final model on all data
    rf.fit(X_train, y_train)

    # Get feature importances
    feature_importances = dict(zip(pc_cols, rf.feature_importances_))

    # Create predictions DataFrame for training samples
    predictions_df = _create_predictions_df(
        scores_df=scores_df,
        model=rf,
        pc_cols=pc_cols,
        min_prob=0.0,  # No threshold for training samples
    )

    return ClassificationResult(
        predictions_df=predictions_df,
        model=rf,
        classes=list(rf.classes_),
        validation_metrics=validation_metrics,
        confusion_matrix=conf_matrix,
        confusion_matrix_labels=conf_matrix_labels,
        feature_importances=feature_importances,
    )


def predict_ancestry(
    model: RandomForestClassifier,
    scores_df: pd.DataFrame,
    n_pcs: int = DEFAULT_N_PCS_CLASSIFY,
    min_prob: float = DEFAULT_MIN_PROB,
) -> pd.DataFrame:
    """Predict ancestry for samples using trained model.

    Parameters
    ----------
    model : RandomForestClassifier
        Trained classifier from train_classifier().
    scores_df : pd.DataFrame
        PC scores for samples to classify.
    n_pcs : int, optional
        Number of PCs to use (must match training). Default: 10.
    min_prob : float, optional
        Minimum probability threshold for assignment. Samples with
        max probability below this are labeled 'unassigned'. Default: 0.75.

    Returns
    -------
    pd.DataFrame
        Predictions with columns:
        - predicted_ancestry: Assigned label or 'unassigned'
        - ancestry_probability: Maximum class probability
        - prob_<CLASS>: Per-class probabilities for each class

    Raises
    ------
    ValueError
        If required columns are missing or min_prob is invalid.

    Example
    -------
    >>> predictions = predict_ancestry(
    ...     model=result.model,
    ...     scores_df=query_scores,
    ...     min_prob=0.75,
    ... )
    >>> print(predictions['predicted_ancestry'].value_counts())
    """
    if not 0 <= min_prob <= 1:
        raise ValueError(f"min_prob must be between 0 and 1, got {min_prob}")

    pc_cols = _get_pc_columns(n_pcs)

    # Validate PC columns exist
    missing_cols = set(pc_cols) - set(scores_df.columns)
    if missing_cols:
        raise ValueError(f"Missing PC columns in scores_df: {sorted(missing_cols)}")

    logger.info(
        f"Predicting ancestry for {len(scores_df)} samples " f"with min_prob={min_prob}"
    )

    return _create_predictions_df(
        scores_df=scores_df,
        model=model,
        pc_cols=pc_cols,
        min_prob=min_prob,
    )


def _create_predictions_df(
    scores_df: pd.DataFrame,
    model: RandomForestClassifier,
    pc_cols: List[str],
    min_prob: float,
) -> pd.DataFrame:
    """Create predictions DataFrame from model and scores.

    Parameters
    ----------
    scores_df : pd.DataFrame
        PC scores DataFrame.
    model : RandomForestClassifier
        Trained classifier.
    pc_cols : List[str]
        List of PC column names to use.
    min_prob : float
        Minimum probability threshold.

    Returns
    -------
    pd.DataFrame
        Predictions DataFrame.
    """
    X = scores_df[pc_cols].values
    classes = model.classes_

    # Get probabilities
    probs = model.predict_proba(X)
    max_probs = probs.max(axis=1)
    pred_indices = probs.argmax(axis=1)

    # Apply threshold
    predictions = np.where(
        max_probs >= min_prob,
        classes[pred_indices],
        "unassigned",
    )

    # Build result DataFrame
    result = pd.DataFrame(index=scores_df.index)
    result[PREDICTED_ANCESTRY_COL] = predictions
    result[ANCESTRY_PROB_COL] = max_probs

    # Add per-class probabilities
    for i, cls in enumerate(classes):
        result[f"prob_{cls}"] = probs[:, i]

    # Log summary
    pred_counts = pd.Series(predictions).value_counts()
    logger.info("Prediction summary:")
    for pop, count in pred_counts.items():
        pct = 100 * count / len(predictions)
        logger.info(f"  {pop}: {count} ({pct:.1f}%)")

    return result


def get_training_samples(
    scores_df: pd.DataFrame,
    source_col: str,
    ancestry_col: str,
    source_value: str = "reference",
) -> Tuple[pd.DataFrame, pd.Series]:
    """Extract training samples from combined scores DataFrame.

    Parameters
    ----------
    scores_df : pd.DataFrame
        Combined scores DataFrame with all samples.
    source_col : str
        Column indicating sample source (query vs reference).
    ancestry_col : str
        Column containing known ancestry labels.
    source_value : str, optional
        Value indicating reference samples. Default: 'reference'.

    Returns
    -------
    Tuple[pd.DataFrame, pd.Series]
        (training_scores, training_labels) tuple.

    Raises
    ------
    ValueError
        If no training samples found or required columns missing.
    """
    # Validate columns exist
    if source_col not in scores_df.columns:
        raise ValueError(f"source column '{source_col}' not found in DataFrame")
    if ancestry_col not in scores_df.columns:
        raise ValueError(f"ancestry column '{ancestry_col}' not found in DataFrame")

    # Filter to reference samples
    mask = scores_df[source_col] == source_value
    train_df = scores_df[mask].copy()

    if len(train_df) == 0:
        raise ValueError(
            f"No training samples found with {source_col}='{source_value}'"
        )

    # Get labels, dropping any missing
    labels = train_df[ancestry_col].dropna()

    # Filter scores to match
    train_scores = train_df.loc[labels.index]

    logger.info(f"Extracted {len(train_scores)} training samples with known ancestry")

    return train_scores, labels


def get_query_samples(
    scores_df: pd.DataFrame,
    source_col: str,
    source_value: str = "query",
) -> pd.DataFrame:
    """Extract query samples from combined scores DataFrame.

    Parameters
    ----------
    scores_df : pd.DataFrame
        Combined scores DataFrame with all samples.
    source_col : str
        Column indicating sample source (query vs reference).
    source_value : str, optional
        Value indicating query samples. Default: 'query'.

    Returns
    -------
    pd.DataFrame
        Query samples scores.

    Raises
    ------
    ValueError
        If no query samples found.
    """
    if source_col not in scores_df.columns:
        raise ValueError(f"source column '{source_col}' not found in DataFrame")

    mask = scores_df[source_col] == source_value
    query_df = scores_df[mask].copy()

    if len(query_df) == 0:
        raise ValueError(f"No query samples found with {source_col}='{source_value}'")

    logger.info(f"Extracted {len(query_df)} query samples for prediction")

    return query_df
