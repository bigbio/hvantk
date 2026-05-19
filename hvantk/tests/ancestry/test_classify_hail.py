"""Integration tests for ancestry classification with Hail.

These tests combine PCA computation with classification to verify
end-to-end functionality with synthetic genetic data.
"""

import pytest

from hvantk.algorithms.ancestry.classify import (
    get_query_samples,
    get_training_samples,
    predict_ancestry,
    train_classifier,
)
from hvantk.algorithms.ancestry.constants import (
    ANCESTRY_PROB_COL,
    KNOWN_ANCESTRY_COL,
    PREDICTED_ANCESTRY_COL,
    SOURCE_COL,
)
from hvantk.algorithms.ancestry.filter import filter_variants_for_ancestry
from hvantk.algorithms.ancestry.pca import compute_pca


@pytest.mark.hail
class TestClassificationWithPCA:
    """Integration tests combining PCA with classification."""

    def test_train_on_reference_pca(self, hail_session, synthetic_reference_mt):
        """Should train classifier on PCA scores from reference panel."""
        # Filter and compute PCA
        mt = filter_variants_for_ancestry(synthetic_reference_mt)
        pca_result = compute_pca(mt, n_pcs=10)

        # Get scores and labels
        scores_df = pca_result.get_scores_df()
        samples_df = synthetic_reference_mt.cols().to_pandas()

        # Merge to get labels
        scores_with_labels = scores_df.merge(samples_df[["s", "ancestry"]], on="s")
        scores_with_labels = scores_with_labels.set_index("s")

        # Train classifier
        result = train_classifier(
            scores_df=scores_with_labels,
            labels=scores_with_labels["ancestry"],
            n_pcs=5,
            validate=True,
            n_cv_folds=3,
        )

        # Should achieve reasonable accuracy on well-separated synthetic data
        accuracy = result.get_accuracy()
        assert accuracy is not None
        assert accuracy > 0.7, f"Expected accuracy > 0.7, got {accuracy}"

    def test_predict_query_samples(self, hail_session, small_merged_mt):
        """Should predict ancestry for query samples in merged MT."""
        # Filter and compute PCA
        mt = filter_variants_for_ancestry(small_merged_mt)
        pca_result = compute_pca(mt, n_pcs=10)

        # Get scores and annotate with source info
        scores_df = pca_result.get_scores_df()
        samples_df = small_merged_mt.cols().to_pandas()

        # Merge source and ancestry info
        scores_with_meta = scores_df.merge(
            samples_df[["s", SOURCE_COL, KNOWN_ANCESTRY_COL]],
            on="s",
        )
        scores_with_meta = scores_with_meta.set_index("s")

        # Split into training and query
        train_scores, train_labels = get_training_samples(
            scores_with_meta,
            source_col=SOURCE_COL,
            ancestry_col=KNOWN_ANCESTRY_COL,
        )

        # Train classifier
        train_result = train_classifier(
            scores_df=train_scores,
            labels=train_labels,
            n_pcs=5,
            validate=False,
        )

        # Get query samples
        query_scores = get_query_samples(scores_with_meta, source_col=SOURCE_COL)

        # Predict
        predictions = predict_ancestry(
            model=train_result.model,
            scores_df=query_scores,
            n_pcs=5,
            min_prob=0.5,
        )

        # Should have predictions for all query samples
        assert len(predictions) == len(query_scores)
        assert PREDICTED_ANCESTRY_COL in predictions.columns
        assert ANCESTRY_PROB_COL in predictions.columns

        # Check that predictions are valid population labels or 'unassigned'
        valid_labels = set(train_result.classes) | {"unassigned"}
        assert all(p in valid_labels for p in predictions[PREDICTED_ANCESTRY_COL])


@pytest.mark.hail
class TestClassificationEdgeCases:
    """Tests for edge cases in classification with real Hail data."""

    def test_single_population_raises(
        self, hail_session, synthetic_reference_mt_single_pop
    ):
        """Training with single population should raise error."""
        # Filter and compute PCA
        mt = filter_variants_for_ancestry(synthetic_reference_mt_single_pop)
        pca_result = compute_pca(mt, n_pcs=5)

        scores_df = pca_result.get_scores_df()
        samples_df = synthetic_reference_mt_single_pop.cols().to_pandas()

        scores_with_labels = scores_df.merge(samples_df[["s", "ancestry"]], on="s")
        scores_with_labels = scores_with_labels.set_index("s")

        # Should raise because only one population
        with pytest.raises(ValueError, match="At least 2 populations"):
            train_classifier(
                scores_df=scores_with_labels,
                labels=scores_with_labels["ancestry"],
                n_pcs=5,
            )


@pytest.mark.hail
@pytest.mark.slow
class TestFullAncestryWorkflow:
    """Integration tests for the full ancestry inference workflow."""

    def test_full_workflow(
        self, hail_session, synthetic_query_mt, synthetic_reference_mt
    ):
        """Test complete workflow: merge -> filter -> PCA -> train -> predict."""
        from hvantk.algorithms.ancestry.merge import merge_matrixtables

        # Step 1: Merge
        merged_mt = merge_matrixtables(
            query_mt=synthetic_query_mt,
            reference_mt=synthetic_reference_mt,
            ancestry_col="ancestry",
            min_shared_variants=100,  # Lower threshold for testing
        )

        # Step 2: Filter
        filtered_mt = filter_variants_for_ancestry(merged_mt)
        n_variants = filtered_mt.count_rows()
        assert n_variants > 0

        # Step 3: PCA
        pca_result = compute_pca(filtered_mt, n_pcs=10)
        assert len(pca_result.eigenvalues) == 10

        # Step 4: Prepare training data
        scores_df = pca_result.get_scores_df()
        samples_df = merged_mt.cols().to_pandas()

        scores_with_meta = scores_df.merge(
            samples_df[["s", SOURCE_COL, KNOWN_ANCESTRY_COL]],
            on="s",
        )
        scores_with_meta = scores_with_meta.set_index("s")

        train_scores, train_labels = get_training_samples(
            scores_with_meta,
            source_col=SOURCE_COL,
            ancestry_col=KNOWN_ANCESTRY_COL,
        )

        # Step 5: Train
        train_result = train_classifier(
            scores_df=train_scores,
            labels=train_labels,
            n_pcs=5,
            validate=True,
            n_cv_folds=3,
        )

        # Verify training
        assert train_result.validation_metrics is not None
        accuracy = train_result.get_accuracy()
        assert accuracy > 0.5  # Should be better than random

        # Step 6: Predict query samples
        query_scores = get_query_samples(scores_with_meta, source_col=SOURCE_COL)

        predictions = predict_ancestry(
            model=train_result.model,
            scores_df=query_scores,
            n_pcs=5,
            min_prob=0.5,
        )

        # Verify predictions
        n_query = synthetic_query_mt.count_cols()
        assert len(predictions) == n_query

        # Count predictions by ancestry
        pred_counts = predictions[PREDICTED_ANCESTRY_COL].value_counts()
        assert len(pred_counts) > 0

    def test_workflow_with_different_thresholds(self, hail_session, small_merged_mt):
        """Test prediction behavior with different probability thresholds."""
        # Filter and compute PCA
        mt = filter_variants_for_ancestry(small_merged_mt)
        pca_result = compute_pca(mt, n_pcs=10)

        scores_df = pca_result.get_scores_df()
        samples_df = small_merged_mt.cols().to_pandas()

        scores_with_meta = scores_df.merge(
            samples_df[["s", SOURCE_COL, KNOWN_ANCESTRY_COL]],
            on="s",
        )
        scores_with_meta = scores_with_meta.set_index("s")

        train_scores, train_labels = get_training_samples(
            scores_with_meta,
            source_col=SOURCE_COL,
            ancestry_col=KNOWN_ANCESTRY_COL,
        )

        train_result = train_classifier(
            scores_df=train_scores,
            labels=train_labels,
            n_pcs=5,
            validate=False,
        )

        query_scores = get_query_samples(scores_with_meta, source_col=SOURCE_COL)

        # Low threshold - more assignments
        pred_low = predict_ancestry(
            train_result.model, query_scores, n_pcs=5, min_prob=0.5
        )
        n_assigned_low = (pred_low[PREDICTED_ANCESTRY_COL] != "unassigned").sum()

        # High threshold - fewer assignments
        pred_high = predict_ancestry(
            train_result.model, query_scores, n_pcs=5, min_prob=0.9
        )
        n_assigned_high = (pred_high[PREDICTED_ANCESTRY_COL] != "unassigned").sum()

        # Higher threshold should have fewer or equal assignments
        assert n_assigned_high <= n_assigned_low
