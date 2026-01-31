"""Tests for ancestry plotting functions.

Tests the visualization functions in hvantk.ancestry.plot without requiring
Hail or actual genetic data. Uses mock DataFrames to test plotting logic.
"""

import pytest
import numpy as np
import pandas as pd

from hvantk.ancestry.constants import (
    PREDICTED_ANCESTRY_COL,
    ANCESTRY_PROB_COL,
    SOURCE_COL,
    KNOWN_ANCESTRY_COL,
)


@pytest.fixture
def sample_predictions_df():
    """Create sample predictions DataFrame for testing plots."""
    np.random.seed(42)
    n_query = 50
    n_ref = 90

    # Generate PC coordinates - 3 clusters for EUR, AFR, EAS
    pc1_means = {"EUR": -5, "AFR": 5, "EAS": 0}
    pc2_means = {"EUR": 0, "AFR": 2, "EAS": -5}

    data = []

    # Reference samples (30 each from EUR, AFR, EAS)
    for pop in ["EUR", "AFR", "EAS"]:
        for i in range(30):
            data.append({
                "s": f"ref_{pop}_{i}",
                "PC1": pc1_means[pop] + np.random.normal(0, 1),
                "PC2": pc2_means[pop] + np.random.normal(0, 1),
                "PC3": np.random.normal(0, 1),
                PREDICTED_ANCESTRY_COL: pop,
                ANCESTRY_PROB_COL: np.random.uniform(0.85, 0.99),
                KNOWN_ANCESTRY_COL: pop,
                SOURCE_COL: "reference",
            })

    # Query samples (randomly distributed)
    pops = ["EUR", "AFR", "EAS"]
    probs = [0.6, 0.3, 0.1]
    for i in range(n_query):
        assigned_pop = np.random.choice(pops, p=probs)
        prob = np.random.uniform(0.5, 0.95)
        predicted = assigned_pop if prob >= 0.75 else "unassigned"

        data.append({
            "s": f"query_{i}",
            "PC1": pc1_means[assigned_pop] + np.random.normal(0, 2),
            "PC2": pc2_means[assigned_pop] + np.random.normal(0, 2),
            "PC3": np.random.normal(0, 1),
            PREDICTED_ANCESTRY_COL: predicted,
            ANCESTRY_PROB_COL: prob,
            KNOWN_ANCESTRY_COL: None,
            SOURCE_COL: "query",
        })

    return pd.DataFrame(data)


@pytest.fixture
def sample_eigenvalues():
    """Create sample eigenvalues for variance explained plot."""
    # Eigenvalues that decrease exponentially
    eigenvalues = [10 * (0.7 ** i) for i in range(20)]
    return eigenvalues


@pytest.fixture
def sample_cv_labels():
    """Create sample cross-validation labels for confusion matrix."""
    np.random.seed(42)
    classes = ["EUR", "AFR", "EAS"]
    n_samples = 90

    y_true = np.array(classes * 30)
    # Simulate 90% accuracy
    y_pred = y_true.copy()
    mistakes = np.random.choice(n_samples, size=9, replace=False)
    for idx in mistakes:
        wrong_classes = [c for c in classes if c != y_true[idx]]
        y_pred[idx] = np.random.choice(wrong_classes)

    return y_true, y_pred, classes


class TestPlotPCAScatter:
    """Tests for plot_pca_scatter function."""

    def test_basic_scatter_plot(self, sample_predictions_df):
        """Test basic PCA scatter plot creation."""
        from hvantk.ancestry.plot import plot_pca_scatter

        fig = plot_pca_scatter(sample_predictions_df)
        assert fig is not None

        # Check figure has axes
        assert len(fig.axes) >= 1

    def test_scatter_plot_different_pcs(self, sample_predictions_df):
        """Test PCA scatter with different PC combinations."""
        from hvantk.ancestry.plot import plot_pca_scatter

        fig = plot_pca_scatter(sample_predictions_df, pc_x=1, pc_y=3)
        assert fig is not None

    def test_scatter_plot_query_as_undefined(self, sample_predictions_df):
        """Test showing query samples as undefined."""
        from hvantk.ancestry.plot import plot_pca_scatter

        fig = plot_pca_scatter(
            sample_predictions_df,
            show_query_as_undefined=True,
        )
        assert fig is not None

    def test_scatter_plot_custom_colors(self, sample_predictions_df):
        """Test scatter plot with custom color palette."""
        from hvantk.ancestry.plot import plot_pca_scatter

        custom_colors = {
            "EUR": "#FF0000",
            "AFR": "#00FF00",
            "EAS": "#0000FF",
            "unassigned": "#999999",
        }
        fig = plot_pca_scatter(sample_predictions_df, colors=custom_colors)
        assert fig is not None

    def test_scatter_plot_filter_populations(self, sample_predictions_df):
        """Test filtering to specific populations."""
        from hvantk.ancestry.plot import plot_pca_scatter

        fig = plot_pca_scatter(
            sample_predictions_df,
            filter_populations=["EUR", "AFR"],
        )
        assert fig is not None

    def test_scatter_plot_filter_source(self, sample_predictions_df):
        """Test filtering to query or reference only."""
        from hvantk.ancestry.plot import plot_pca_scatter

        fig = plot_pca_scatter(sample_predictions_df, filter_source="query")
        assert fig is not None

    def test_scatter_plot_no_legend(self, sample_predictions_df):
        """Test scatter plot without legend."""
        from hvantk.ancestry.plot import plot_pca_scatter

        fig = plot_pca_scatter(sample_predictions_df, show_legend=False)
        assert fig is not None


class TestPlotPCAPanel:
    """Tests for plot_pca_panel function."""

    def test_two_panel_plot(self, sample_predictions_df):
        """Test two-panel PCA plot creation."""
        from hvantk.ancestry.plot import plot_pca_panel

        fig = plot_pca_panel(sample_predictions_df)
        assert fig is not None

        # Should have 2 subplots
        assert len(fig.axes) == 2

    def test_panel_plot_query_as_undefined(self, sample_predictions_df):
        """Test two-panel plot with query as undefined."""
        from hvantk.ancestry.plot import plot_pca_panel

        fig = plot_pca_panel(
            sample_predictions_df,
            show_query_as_undefined=True,
        )
        assert fig is not None


class TestPlotVarianceExplained:
    """Tests for plot_variance_explained function."""

    def test_variance_explained_plot(self, sample_eigenvalues):
        """Test variance explained bar plot."""
        from hvantk.ancestry.plot import plot_variance_explained

        fig = plot_variance_explained(sample_eigenvalues)
        assert fig is not None

    def test_variance_explained_limited_pcs(self, sample_eigenvalues):
        """Test limiting number of PCs shown."""
        from hvantk.ancestry.plot import plot_variance_explained

        fig = plot_variance_explained(sample_eigenvalues, n_pcs=5)
        assert fig is not None

    def test_variance_explained_with_cumulative(self, sample_eigenvalues):
        """Test variance plot with cumulative line."""
        from hvantk.ancestry.plot import plot_variance_explained

        fig = plot_variance_explained(sample_eigenvalues, cumulative=True)
        assert fig is not None


class TestPlotAncestryProportions:
    """Tests for plot_ancestry_proportions function."""

    def test_ancestry_proportions_plot(self, sample_predictions_df):
        """Test ancestry proportions bar chart."""
        from hvantk.ancestry.plot import plot_ancestry_proportions

        fig = plot_ancestry_proportions(sample_predictions_df)
        assert fig is not None

    def test_ancestry_proportions_custom_colors(self, sample_predictions_df):
        """Test ancestry proportions with custom colors."""
        from hvantk.ancestry.plot import plot_ancestry_proportions

        custom_colors = {
            "EUR": "#FF0000",
            "AFR": "#00FF00",
            "EAS": "#0000FF",
            "unassigned": "#999999",
        }
        fig = plot_ancestry_proportions(sample_predictions_df, colors=custom_colors)
        assert fig is not None


class TestPlotProbabilityDistribution:
    """Tests for plot_probability_distribution function."""

    def test_probability_distribution_plot(self, sample_predictions_df):
        """Test probability distribution histogram."""
        from hvantk.ancestry.plot import plot_probability_distribution

        fig = plot_probability_distribution(sample_predictions_df)
        assert fig is not None

    def test_probability_distribution_custom_bins(self, sample_predictions_df):
        """Test probability histogram with custom bins."""
        from hvantk.ancestry.plot import plot_probability_distribution

        fig = plot_probability_distribution(sample_predictions_df, bins=20)
        assert fig is not None


class TestPlotConfusionMatrix:
    """Tests for plot_confusion_matrix function."""

    def test_confusion_matrix_plot(self, sample_cv_labels):
        """Test confusion matrix heatmap."""
        from hvantk.ancestry.plot import plot_confusion_matrix

        y_true, y_pred, labels = sample_cv_labels
        fig = plot_confusion_matrix(y_true, y_pred, labels=labels)
        assert fig is not None

    def test_confusion_matrix_normalized(self, sample_cv_labels):
        """Test normalized confusion matrix."""
        from hvantk.ancestry.plot import plot_confusion_matrix

        y_true, y_pred, labels = sample_cv_labels
        fig = plot_confusion_matrix(y_true, y_pred, labels=labels, normalize=True)
        assert fig is not None

    def test_confusion_matrix_unnormalized(self, sample_cv_labels):
        """Test unnormalized confusion matrix."""
        from hvantk.ancestry.plot import plot_confusion_matrix

        y_true, y_pred, labels = sample_cv_labels
        fig = plot_confusion_matrix(y_true, y_pred, labels=labels, normalize=False)
        assert fig is not None


class TestEncodeFigure:
    """Tests for figure encoding functions."""

    def test_encode_figure_to_base64(self, sample_predictions_df):
        """Test encoding figure to base64 string."""
        from hvantk.ancestry.plot import plot_pca_scatter, encode_figure_to_base64

        fig = plot_pca_scatter(sample_predictions_df)
        encoded = encode_figure_to_base64(fig)

        assert encoded is not None
        assert encoded.startswith("data:image/png;base64,")
        assert len(encoded) > 100

    def test_encode_figure_different_dpi(self, sample_predictions_df):
        """Test encoding with different DPI settings."""
        from hvantk.ancestry.plot import plot_pca_scatter, encode_figure_to_base64

        fig = plot_pca_scatter(sample_predictions_df)
        encoded_low = encode_figure_to_base64(fig, dpi=50)
        encoded_high = encode_figure_to_base64(fig, dpi=150)

        # Higher DPI should produce larger output
        assert len(encoded_high) > len(encoded_low)


class TestCloseFigure:
    """Tests for close_figure function."""

    def test_close_figure(self, sample_predictions_df):
        """Test closing figure releases memory."""
        from hvantk.ancestry.plot import plot_pca_scatter, close_figure
        import matplotlib.pyplot as plt

        initial_figs = len(plt.get_fignums())
        fig = plot_pca_scatter(sample_predictions_df)
        assert len(plt.get_fignums()) == initial_figs + 1

        close_figure(fig)
        assert len(plt.get_fignums()) == initial_figs


class TestPlotEdgeCases:
    """Test edge cases and error handling."""

    def test_empty_dataframe(self):
        """Test plotting with empty DataFrame."""
        from hvantk.ancestry.plot import plot_pca_scatter

        empty_df = pd.DataFrame(columns=[
            "s", "PC1", "PC2", PREDICTED_ANCESTRY_COL, SOURCE_COL, KNOWN_ANCESTRY_COL
        ])

        # Should not raise, may produce empty plot
        fig = plot_pca_scatter(empty_df)
        assert fig is not None

    def test_single_population(self, sample_predictions_df):
        """Test plotting with single population."""
        from hvantk.ancestry.plot import plot_pca_scatter

        # Filter to single population
        single_pop = sample_predictions_df[
            sample_predictions_df[PREDICTED_ANCESTRY_COL] == "EUR"
        ].copy()

        fig = plot_pca_scatter(single_pop)
        assert fig is not None

    def test_missing_probability_column(self, sample_predictions_df):
        """Test probability plot with missing column."""
        from hvantk.ancestry.plot import plot_probability_distribution

        df_no_prob = sample_predictions_df.drop(columns=[ANCESTRY_PROB_COL])

        # Should handle gracefully
        with pytest.raises((KeyError, ValueError)):
            plot_probability_distribution(df_no_prob)
