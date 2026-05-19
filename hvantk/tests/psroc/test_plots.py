"""
Tests for PSROC plotting functions.

This module tests the ROC curve visualization functions to ensure
they create valid figures and handle edge cases correctly.
"""

import tempfile
from pathlib import Path

import numpy as np
import pytest
import matplotlib

matplotlib.use("Agg")  # Use non-interactive backend for tests
import matplotlib.pyplot as plt

from hvantk.algorithms.psroc.roc import (
    ScoreMissingness,
    compute_roc_metrics,
    compute_all_missingness,
)
from hvantk.algorithms.psroc.plots import (
    plot_roc_curves,
    plot_roc_curve_single,
    plot_auc_comparison,
    plot_missingness_summary,
    plot_psroc_summary_dashboard,
    plot_collection_heatmap,
)


@pytest.fixture
def sample_roc_results():
    """Create sample ROC results for testing."""
    labels = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1])
    scores = {
        "good_score": np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]),
        "medium_score": np.array([0.2, 0.3, 0.5, 0.6, 0.7, 0.4, 0.5, 0.8, 0.9, 0.95]),
        "poor_score": np.array([0.5, 0.4, 0.6, 0.5, 0.55, 0.45, 0.5, 0.55, 0.6, 0.5]),
    }
    return compute_roc_metrics(labels, scores)


@pytest.fixture
def sample_missingness():
    """Create sample missingness data for testing."""
    return {
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
            n_present=75,
            n_missing=25,
            missingness_rate=0.25,
            included_in_analysis=True,
        ),
        "score_c": ScoreMissingness(
            score_name="score_c",
            n_total=100,
            n_present=50,
            n_missing=50,
            missingness_rate=0.50,
            included_in_analysis=False,
            exclusion_reason="exceeds threshold",
        ),
    }


class TestPlotROCCurves:
    """Test plot_roc_curves function."""

    def test_basic_plot(self, sample_roc_results):
        """Test basic ROC curve plotting."""
        fig = plot_roc_curves(sample_roc_results)

        assert fig is not None
        assert isinstance(fig, plt.Figure)

        # Check that axes exist
        assert len(fig.axes) >= 1

        plt.close(fig)

    def test_plot_with_all_options(self, sample_roc_results):
        """Test plotting with all options enabled."""
        fig = plot_roc_curves(
            sample_roc_results,
            title="Test ROC Curves",
            show_optimal=True,
            show_diagonal=True,
            show_legend=True,
            show_auc_in_legend=True,
            color_by_auc=True,
        )

        assert fig is not None
        plt.close(fig)

    def test_plot_without_legend(self, sample_roc_results):
        """Test plotting without legend."""
        fig = plot_roc_curves(
            sample_roc_results,
            show_legend=False,
        )

        assert fig is not None
        plt.close(fig)

    def test_plot_single_score(self, sample_roc_results):
        """Test plotting single score."""
        single_result = {"good_score": sample_roc_results["good_score"]}
        fig = plot_roc_curves(single_result)

        assert fig is not None
        plt.close(fig)

    def test_save_to_file(self, sample_roc_results):
        """Test saving plot to file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "roc_curves"
            fig = plot_roc_curves(sample_roc_results, output_path=output_path)

            # Check that file was created
            assert (Path(tmpdir) / "roc_curves.png").exists()

            plt.close(fig)

    def test_empty_results_raises_error(self):
        """Test that empty results raise ValueError."""
        with pytest.raises(ValueError, match="No ROC results"):
            plot_roc_curves({})


class TestPlotROCCurveSingle:
    """Test plot_roc_curve_single function."""

    def test_basic_single_plot(self, sample_roc_results):
        """Test basic single ROC curve plotting."""
        fig = plot_roc_curve_single(sample_roc_results["good_score"])

        assert fig is not None
        assert isinstance(fig, plt.Figure)

        plt.close(fig)

    def test_with_all_options(self, sample_roc_results):
        """Test single plot with all options."""
        fig = plot_roc_curve_single(
            sample_roc_results["good_score"],
            title="Custom Title",
            show_optimal=True,
            show_diagonal=True,
            show_stats=True,
        )

        assert fig is not None
        plt.close(fig)

    def test_without_stats(self, sample_roc_results):
        """Test single plot without stats box."""
        fig = plot_roc_curve_single(
            sample_roc_results["good_score"],
            show_stats=False,
        )

        assert fig is not None
        plt.close(fig)

    def test_save_to_file(self, sample_roc_results):
        """Test saving single plot to file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "single_roc"
            fig = plot_roc_curve_single(
                sample_roc_results["good_score"],
                output_path=output_path,
            )

            assert (Path(tmpdir) / "single_roc.png").exists()
            plt.close(fig)


class TestPlotAUCComparison:
    """Test plot_auc_comparison function."""

    def test_horizontal_bar_chart(self, sample_roc_results):
        """Test horizontal bar chart."""
        fig = plot_auc_comparison(sample_roc_results, horizontal=True)

        assert fig is not None
        plt.close(fig)

    def test_vertical_bar_chart(self, sample_roc_results):
        """Test vertical bar chart."""
        fig = plot_auc_comparison(sample_roc_results, horizontal=False)

        assert fig is not None
        plt.close(fig)

    def test_with_values(self, sample_roc_results):
        """Test with value annotations."""
        fig = plot_auc_comparison(sample_roc_results, show_values=True)

        assert fig is not None
        plt.close(fig)

    def test_color_by_auc(self, sample_roc_results):
        """Test coloring by AUC quality."""
        fig = plot_auc_comparison(sample_roc_results, color_by_auc=True)

        assert fig is not None
        plt.close(fig)

    def test_empty_results_raises_error(self):
        """Test that empty results raise ValueError."""
        with pytest.raises(ValueError, match="No ROC results"):
            plot_auc_comparison({})


class TestPlotMissingnessSummary:
    """Test plot_missingness_summary function."""

    def test_basic_missingness_plot(self, sample_missingness):
        """Test basic missingness summary plot."""
        fig = plot_missingness_summary(sample_missingness)

        assert fig is not None
        plt.close(fig)

    def test_with_threshold_line(self, sample_missingness):
        """Test with threshold line."""
        fig = plot_missingness_summary(
            sample_missingness,
            max_missingness_threshold=0.3,
        )

        assert fig is not None
        plt.close(fig)

    def test_save_to_file(self, sample_missingness):
        """Test saving to file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "missingness"
            fig = plot_missingness_summary(
                sample_missingness,
                output_path=output_path,
            )

            assert (Path(tmpdir) / "missingness.png").exists()
            plt.close(fig)

    def test_empty_data_raises_error(self):
        """Test that empty data raises ValueError."""
        with pytest.raises(ValueError, match="No missingness data"):
            plot_missingness_summary({})


class TestPlotPSROCSummaryDashboard:
    """Test plot_psroc_summary_dashboard function."""

    def test_full_dashboard(self, sample_roc_results, sample_missingness):
        """Test full summary dashboard."""
        fig = plot_psroc_summary_dashboard(
            sample_roc_results,
            sample_missingness,
        )

        assert fig is not None
        # Dashboard should have multiple subplots
        assert len(fig.axes) >= 4
        plt.close(fig)

    def test_dashboard_with_threshold(self, sample_roc_results, sample_missingness):
        """Test dashboard with missingness threshold."""
        fig = plot_psroc_summary_dashboard(
            sample_roc_results,
            sample_missingness,
            max_missingness_threshold=0.3,
        )

        assert fig is not None
        plt.close(fig)

    def test_dashboard_empty_results(self, sample_missingness):
        """Test dashboard with empty ROC results."""
        fig = plot_psroc_summary_dashboard({}, sample_missingness)

        assert fig is not None
        plt.close(fig)

    def test_dashboard_empty_missingness(self, sample_roc_results):
        """Test dashboard with empty missingness."""
        fig = plot_psroc_summary_dashboard(sample_roc_results, {})

        assert fig is not None
        plt.close(fig)

    def test_save_dashboard(self, sample_roc_results, sample_missingness):
        """Test saving dashboard to file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "dashboard"
            fig = plot_psroc_summary_dashboard(
                sample_roc_results,
                sample_missingness,
                output_path=output_path,
            )

            assert (Path(tmpdir) / "dashboard.png").exists()
            plt.close(fig)


@pytest.fixture
def sample_collection_metrics():
    """Create sample collection metrics (multi-panel) for testing."""
    labels = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1])
    scores_a = {
        "CADD_phred": np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]),
        "REVEL_score": np.array([0.2, 0.3, 0.5, 0.6, 0.7, 0.4, 0.5, 0.8, 0.9, 0.95]),
    }
    scores_b = {
        "CADD_phred": np.array([0.3, 0.2, 0.4, 0.5, 0.45, 0.55, 0.65, 0.75, 0.85, 0.9]),
        "REVEL_score": np.array([0.1, 0.15, 0.2, 0.3, 0.35, 0.7, 0.8, 0.85, 0.9, 0.95]),
    }
    scores_c = {
        "CADD_phred": np.array([0.5, 0.4, 0.6, 0.5, 0.55, 0.45, 0.5, 0.55, 0.6, 0.5]),
        "REVEL_score": np.array([0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.75, 0.8]),
    }
    return {
        "Hereditary Cancer": compute_roc_metrics(labels, scores_a),
        "Cardiomyopathy": compute_roc_metrics(labels, scores_b),
        "Epilepsy": compute_roc_metrics(labels, scores_c),
    }


class TestPlotCollectionHeatmap:
    """Test plot_collection_heatmap function."""

    def test_basic_heatmap(self, sample_collection_metrics):
        """Test basic heatmap creation."""
        fig = plot_collection_heatmap(sample_collection_metrics)

        assert fig is not None
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_heatmap_without_ci(self, sample_collection_metrics):
        """Test heatmap without CI annotations."""
        fig = plot_collection_heatmap(sample_collection_metrics, show_ci=False)

        assert fig is not None
        plt.close(fig)

    def test_heatmap_without_values(self, sample_collection_metrics):
        """Test heatmap without value annotations."""
        fig = plot_collection_heatmap(sample_collection_metrics, show_values=False)

        assert fig is not None
        plt.close(fig)

    def test_heatmap_sort_by_name(self, sample_collection_metrics):
        """Test heatmap with scores sorted alphabetically."""
        fig = plot_collection_heatmap(sample_collection_metrics, sort_scores_by="name")

        assert fig is not None
        plt.close(fig)

    def test_heatmap_custom_figsize(self, sample_collection_metrics):
        """Test heatmap with custom figure size."""
        fig = plot_collection_heatmap(sample_collection_metrics, figsize=(14, 6))

        assert fig is not None
        plt.close(fig)

    def test_save_to_file(self, sample_collection_metrics):
        """Test saving heatmap to file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "heatmap"
            fig = plot_collection_heatmap(
                sample_collection_metrics, output_path=output_path
            )

            assert (Path(tmpdir) / "heatmap.png").exists()
            plt.close(fig)

    def test_empty_collection_raises_error(self):
        """Test that empty collection raises ValueError."""
        with pytest.raises(ValueError, match="No collection metrics"):
            plot_collection_heatmap({})

    def test_no_scores_raises_error(self):
        """Test that collection with no scores raises ValueError."""
        with pytest.raises(ValueError, match="No scores found"):
            plot_collection_heatmap({"group_a": {}, "group_b": {}})

    def test_missing_score_in_group(self):
        """Test heatmap where a score is missing from one group."""
        labels = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1])
        scores_full = {
            "CADD_phred": np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]),
            "REVEL_score": np.array(
                [0.2, 0.3, 0.5, 0.6, 0.7, 0.4, 0.5, 0.8, 0.9, 0.95]
            ),
        }
        scores_partial = {
            "CADD_phred": np.array(
                [0.3, 0.2, 0.4, 0.5, 0.45, 0.55, 0.65, 0.75, 0.85, 0.9]
            ),
        }

        collection = {
            "Panel_A": compute_roc_metrics(labels, scores_full),
            "Panel_B": compute_roc_metrics(labels, scores_partial),
        }

        fig = plot_collection_heatmap(collection)
        assert fig is not None
        plt.close(fig)

    def test_with_bootstrap_ci(self):
        """Test heatmap with bootstrap CI values populated."""
        labels = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1])
        scores = {
            "score_a": np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]),
        }

        collection = {
            "Group_1": compute_roc_metrics(labels, scores, n_bootstrap=100),
            "Group_2": compute_roc_metrics(labels, scores, n_bootstrap=100),
        }

        # Verify CI values exist
        for metrics in collection.values():
            for roc in metrics.values():
                assert roc.auc_ci_lower is not None
                assert roc.auc_ci_upper is not None

        fig = plot_collection_heatmap(collection, show_ci=True)
        assert fig is not None
        plt.close(fig)


class TestPlotIntegration:
    """Integration tests for plotting functions."""

    def test_end_to_end_workflow(self):
        """Test complete workflow from data to plots."""
        # Create test data
        np.random.seed(42)
        n = 100
        labels = np.array([0] * (n // 2) + [1] * (n // 2))

        # Create scores with different qualities
        scores = {
            "excellent": np.concatenate(
                [
                    np.random.normal(0.3, 0.1, n // 2),
                    np.random.normal(0.8, 0.1, n // 2),
                ]
            ),
            "good": np.concatenate(
                [
                    np.random.normal(0.4, 0.15, n // 2),
                    np.random.normal(0.7, 0.15, n // 2),
                ]
            ),
            "with_missing": np.concatenate(
                [
                    np.random.normal(0.35, 0.12, n // 2),
                    np.random.normal(0.75, 0.12, n // 2),
                ]
            ),
        }

        # Add some missing values
        scores["with_missing"][np.random.choice(n, 10, replace=False)] = np.nan

        # Compute ROC metrics
        results = compute_roc_metrics(labels, scores, max_missingness=0.5)
        missingness = compute_all_missingness(scores, max_missingness=0.3)

        # Test all plotting functions
        with tempfile.TemporaryDirectory() as tmpdir:
            # ROC curves
            fig1 = plot_roc_curves(
                results,
                output_path=Path(tmpdir) / "roc_curves",
            )
            assert (Path(tmpdir) / "roc_curves.png").exists()
            plt.close(fig1)

            # Single curve
            fig2 = plot_roc_curve_single(
                results["excellent"],
                output_path=Path(tmpdir) / "single",
            )
            assert (Path(tmpdir) / "single.png").exists()
            plt.close(fig2)

            # AUC comparison
            fig3 = plot_auc_comparison(
                results,
                output_path=Path(tmpdir) / "auc",
            )
            assert (Path(tmpdir) / "auc.png").exists()
            plt.close(fig3)

            # Missingness summary
            fig4 = plot_missingness_summary(
                missingness,
                output_path=Path(tmpdir) / "missingness",
                max_missingness_threshold=0.3,
            )
            assert (Path(tmpdir) / "missingness.png").exists()
            plt.close(fig4)

            # Dashboard
            fig5 = plot_psroc_summary_dashboard(
                results,
                missingness,
                output_path=Path(tmpdir) / "dashboard",
                max_missingness_threshold=0.3,
            )
            assert (Path(tmpdir) / "dashboard.png").exists()
            plt.close(fig5)


if __name__ == "__main__":
    print("Running PSROC plotting tests...")

    # Create test data
    labels = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1])
    scores = {
        "good_score": np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]),
        "medium_score": np.array([0.2, 0.3, 0.5, 0.6, 0.7, 0.4, 0.5, 0.8, 0.9, 0.95]),
    }
    results = compute_roc_metrics(labels, scores)
    missingness = compute_all_missingness(scores)

    print("\n1. Testing plot_roc_curves")
    fig = plot_roc_curves(results)
    print("  ✓ Created ROC curves plot")
    plt.close(fig)

    print("\n2. Testing plot_roc_curve_single")
    fig = plot_roc_curve_single(results["good_score"])
    print("  ✓ Created single ROC curve plot")
    plt.close(fig)

    print("\n3. Testing plot_auc_comparison")
    fig = plot_auc_comparison(results)
    print("  ✓ Created AUC comparison plot")
    plt.close(fig)

    print("\n4. Testing plot_missingness_summary")
    fig = plot_missingness_summary(missingness)
    print("  ✓ Created missingness summary plot")
    plt.close(fig)

    print("\n5. Testing plot_psroc_summary_dashboard")
    fig = plot_psroc_summary_dashboard(results, missingness)
    print("  ✓ Created summary dashboard")
    plt.close(fig)

    print("\n✅ All plotting tests passed!")
