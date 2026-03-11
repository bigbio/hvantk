"""Unit tests for ancestry classification module.

These tests do not require Hail and use synthetic numpy/pandas data.
"""

import numpy as np
import pandas as pd
import pytest
from sklearn.ensemble import RandomForestClassifier

from hvantk.ancestry.classify import (
    ClassificationResult,
    _get_pc_columns,
    get_query_samples,
    get_training_samples,
    predict_ancestry,
    train_classifier,
    validate_training_data,
)
from hvantk.ancestry.constants import (
    ANCESTRY_PROB_COL,
    PREDICTED_ANCESTRY_COL,
)


@pytest.fixture
def synthetic_scores():
    """Create synthetic PC scores for testing."""
    np.random.seed(42)

    # Create 3 distinct population clusters
    n_per_pop = 30

    # EUR: centered around (2, 1)
    eur_pc1 = np.random.normal(2, 0.5, n_per_pop)
    eur_pc2 = np.random.normal(1, 0.5, n_per_pop)

    # AFR: centered around (-2, 0)
    afr_pc1 = np.random.normal(-2, 0.5, n_per_pop)
    afr_pc2 = np.random.normal(0, 0.5, n_per_pop)

    # EAS: centered around (0, -2)
    eas_pc1 = np.random.normal(0, 0.5, n_per_pop)
    eas_pc2 = np.random.normal(-2, 0.5, n_per_pop)

    # Combine
    pc1 = np.concatenate([eur_pc1, afr_pc1, eas_pc1])
    pc2 = np.concatenate([eur_pc2, afr_pc2, eas_pc2])

    df = pd.DataFrame(
        {
            "PC1": pc1,
            "PC2": pc2,
            "PC3": np.random.randn(n_per_pop * 3),
            "PC4": np.random.randn(n_per_pop * 3),
            "PC5": np.random.randn(n_per_pop * 3),
        },
        index=[f"sample_{i}" for i in range(n_per_pop * 3)],
    )

    return df


@pytest.fixture
def synthetic_labels():
    """Create labels matching synthetic_scores."""
    n_per_pop = 30
    return pd.Series(
        ["EUR"] * n_per_pop + ["AFR"] * n_per_pop + ["EAS"] * n_per_pop,
        index=[f"sample_{i}" for i in range(n_per_pop * 3)],
    )


@pytest.fixture
def trained_model(synthetic_scores, synthetic_labels):
    """Create a pre-trained classifier for testing predictions."""
    result = train_classifier(
        scores_df=synthetic_scores,
        labels=synthetic_labels,
        n_pcs=2,
        validate=False,
    )
    return result.model


class TestGetPcColumns:
    """Tests for _get_pc_columns helper function."""

    def test_generates_correct_columns(self):
        """Should generate PC1, PC2, etc."""
        cols = _get_pc_columns(5)
        assert cols == ["PC1", "PC2", "PC3", "PC4", "PC5"]

    def test_single_pc(self):
        """Should work with single PC."""
        cols = _get_pc_columns(1)
        assert cols == ["PC1"]


class TestValidateTrainingData:
    """Tests for validate_training_data function."""

    def test_valid_data_passes(self, synthetic_labels):
        """Valid training data should pass validation."""
        pop_counts = validate_training_data(synthetic_labels)
        assert "EUR" in pop_counts
        assert "AFR" in pop_counts
        assert "EAS" in pop_counts
        assert all(count >= 10 for count in pop_counts.values())

    def test_single_population_raises(self):
        """Single population should raise ValueError."""
        labels = pd.Series(["EUR"] * 50)
        with pytest.raises(ValueError, match="At least 2 populations"):
            validate_training_data(labels)

    def test_insufficient_samples_raises(self):
        """Too few samples per population should raise ValueError."""
        labels = pd.Series(["EUR"] * 5 + ["AFR"] * 30)
        with pytest.raises(ValueError, match="insufficient samples"):
            validate_training_data(labels, min_samples_per_pop=10)

    def test_custom_min_samples(self):
        """Should respect custom min_samples_per_pop."""
        labels = pd.Series(["EUR"] * 5 + ["AFR"] * 5)
        # Should pass with lower threshold
        pop_counts = validate_training_data(labels, min_samples_per_pop=5)
        assert pop_counts["EUR"] == 5
        assert pop_counts["AFR"] == 5


class TestTrainClassifier:
    """Tests for train_classifier function."""

    def test_basic_training(self, synthetic_scores, synthetic_labels):
        """Basic training should succeed."""
        result = train_classifier(
            scores_df=synthetic_scores,
            labels=synthetic_labels,
            n_pcs=2,
            validate=False,
        )

        assert isinstance(result, ClassificationResult)
        assert isinstance(result.model, RandomForestClassifier)
        assert len(result.classes) == 3
        assert set(result.classes) == {"EUR", "AFR", "EAS"}

    def test_training_with_validation(self, synthetic_scores, synthetic_labels):
        """Training with cross-validation should produce metrics."""
        result = train_classifier(
            scores_df=synthetic_scores,
            labels=synthetic_labels,
            n_pcs=2,
            validate=True,
            n_cv_folds=3,
        )

        assert result.validation_metrics is not None
        assert "accuracy" in result.validation_metrics
        assert result.confusion_matrix is not None
        assert result.confusion_matrix.shape == (3, 3)

    def test_feature_importances(self, synthetic_scores, synthetic_labels):
        """Should compute feature importances."""
        result = train_classifier(
            scores_df=synthetic_scores,
            labels=synthetic_labels,
            n_pcs=2,
            validate=False,
        )

        assert result.feature_importances is not None
        assert "PC1" in result.feature_importances
        assert "PC2" in result.feature_importances
        assert sum(result.feature_importances.values()) == pytest.approx(1.0)

    def test_missing_pc_columns_raises(self, synthetic_labels):
        """Missing PC columns should raise ValueError."""
        scores = pd.DataFrame({"X": [1, 2, 3], "Y": [4, 5, 6]})
        labels = synthetic_labels[:3]

        with pytest.raises(ValueError, match="Missing PC columns"):
            train_classifier(scores, labels, n_pcs=2)

    def test_predictions_df_created(self, synthetic_scores, synthetic_labels):
        """Should create predictions DataFrame for training samples."""
        result = train_classifier(
            scores_df=synthetic_scores,
            labels=synthetic_labels,
            n_pcs=2,
            validate=False,
        )

        assert PREDICTED_ANCESTRY_COL in result.predictions_df.columns
        assert ANCESTRY_PROB_COL in result.predictions_df.columns
        assert len(result.predictions_df) == len(synthetic_scores)


class TestPredictAncestry:
    """Tests for predict_ancestry function."""

    def test_basic_prediction(self, trained_model, synthetic_scores):
        """Basic prediction should succeed."""
        predictions = predict_ancestry(
            model=trained_model,
            scores_df=synthetic_scores,
            n_pcs=2,
            min_prob=0.5,
        )

        assert PREDICTED_ANCESTRY_COL in predictions.columns
        assert ANCESTRY_PROB_COL in predictions.columns
        assert len(predictions) == len(synthetic_scores)

    def test_probability_threshold_respected(self, trained_model):
        """Samples below min_prob should be 'unassigned'."""
        # Create samples at boundary between clusters (low confidence)
        ambiguous = pd.DataFrame(
            {
                "PC1": [0.0, 0.0, 0.0],  # Center point, equidistant from all
                "PC2": [0.0, -0.5, 0.5],
            }
        )

        # High threshold should cause unassigned
        predictions_high = predict_ancestry(
            model=trained_model,
            scores_df=ambiguous,
            n_pcs=2,
            min_prob=0.99,
        )

        # Deterministically check: rows with max_prob < 0.99 should be unassigned
        prob_cols = [c for c in predictions_high.columns if c.startswith("prob_")]
        max_probs = predictions_high[prob_cols].max(axis=1)

        # Assert rows below threshold are unassigned
        low_conf_mask = max_probs < 0.99
        assert (
            predictions_high.loc[low_conf_mask, PREDICTED_ANCESTRY_COL] == "unassigned"
        ).all()

        # Assert rows at/above threshold are NOT unassigned (if any exist)
        high_conf_mask = max_probs >= 0.99
        if high_conf_mask.any():
            assert (
                predictions_high.loc[high_conf_mask, PREDICTED_ANCESTRY_COL]
                != "unassigned"
            ).all()

    def test_probability_columns_created(self, trained_model, synthetic_scores):
        """Should create per-class probability columns."""
        predictions = predict_ancestry(
            model=trained_model,
            scores_df=synthetic_scores,
            n_pcs=2,
        )

        # Check prob columns exist for each class
        assert "prob_EUR" in predictions.columns
        assert "prob_AFR" in predictions.columns
        assert "prob_EAS" in predictions.columns

        # Probabilities should sum to 1
        prob_cols = [c for c in predictions.columns if c.startswith("prob_")]
        prob_sums = predictions[prob_cols].sum(axis=1)
        assert all(prob_sums.apply(lambda x: x == pytest.approx(1.0)))

    def test_invalid_min_prob_raises(self, trained_model, synthetic_scores):
        """Invalid min_prob should raise ValueError."""
        with pytest.raises(ValueError, match="min_prob must be between"):
            predict_ancestry(trained_model, synthetic_scores, n_pcs=2, min_prob=-0.1)

        with pytest.raises(ValueError, match="min_prob must be between"):
            predict_ancestry(trained_model, synthetic_scores, n_pcs=2, min_prob=1.5)

    def test_missing_columns_raises(self, trained_model):
        """Missing PC columns should raise ValueError."""
        scores = pd.DataFrame({"X": [1, 2], "Y": [3, 4]})

        with pytest.raises(ValueError, match="Missing PC columns"):
            predict_ancestry(trained_model, scores, n_pcs=2)


class TestGetTrainingSamples:
    """Tests for get_training_samples function."""

    def test_extracts_reference_samples(self):
        """Should extract reference samples correctly."""
        df = pd.DataFrame(
            {
                "PC1": [1, 2, 3, 4],
                "PC2": [1, 2, 3, 4],
                "_source": ["reference", "reference", "query", "query"],
                "_ancestry": ["EUR", "AFR", None, None],
            },
            index=["s1", "s2", "s3", "s4"],
        )

        scores, labels = get_training_samples(
            df, source_col="_source", ancestry_col="_ancestry"
        )

        assert len(scores) == 2
        assert len(labels) == 2
        assert set(labels.values) == {"EUR", "AFR"}

    def test_missing_source_col_raises(self):
        """Missing source column should raise ValueError."""
        df = pd.DataFrame({"PC1": [1, 2], "ancestry": ["EUR", "AFR"]})

        with pytest.raises(ValueError, match="source column"):
            get_training_samples(df, source_col="_source", ancestry_col="ancestry")

    def test_no_training_samples_raises(self):
        """No training samples should raise ValueError."""
        df = pd.DataFrame(
            {
                "PC1": [1, 2],
                "_source": ["query", "query"],
                "_ancestry": [None, None],
            }
        )

        with pytest.raises(ValueError, match="No training samples"):
            get_training_samples(df, source_col="_source", ancestry_col="_ancestry")


class TestGetQuerySamples:
    """Tests for get_query_samples function."""

    def test_extracts_query_samples(self):
        """Should extract query samples correctly."""
        df = pd.DataFrame(
            {
                "PC1": [1, 2, 3, 4],
                "PC2": [1, 2, 3, 4],
                "_source": ["reference", "reference", "query", "query"],
            },
            index=["s1", "s2", "s3", "s4"],
        )

        query = get_query_samples(df, source_col="_source")

        assert len(query) == 2
        assert set(query.index) == {"s3", "s4"}

    def test_no_query_samples_raises(self):
        """No query samples should raise ValueError."""
        df = pd.DataFrame(
            {
                "PC1": [1, 2],
                "_source": ["reference", "reference"],
            }
        )

        with pytest.raises(ValueError, match="No query samples"):
            get_query_samples(df, source_col="_source")


class TestClassificationResult:
    """Tests for ClassificationResult dataclass."""

    def test_get_accuracy(self, synthetic_scores, synthetic_labels):
        """Should return accuracy from validation metrics."""
        result = train_classifier(
            scores_df=synthetic_scores,
            labels=synthetic_labels,
            n_pcs=2,
            validate=True,
            n_cv_folds=3,
        )

        accuracy = result.get_accuracy()
        assert accuracy is not None
        assert 0 <= accuracy <= 1

    def test_get_accuracy_no_validation(self, synthetic_scores, synthetic_labels):
        """Should return None when no validation performed."""
        result = train_classifier(
            scores_df=synthetic_scores,
            labels=synthetic_labels,
            n_pcs=2,
            validate=False,
        )

        assert result.get_accuracy() is None
