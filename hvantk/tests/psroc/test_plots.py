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

from hvantk.psroc.roc import (
    ROCResult,
    ScoreMissingness,
    compute_roc_metrics,
    compute_all_missingness,
)
from hvantk.psroc.plots import (
    plot_roc_curves,
    plot_roc_curve_single,
    plot_auc_comparison,
    plot_missingness_summary,
    plot_psroc_summary_dashboard,
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
            "excellent": np.concatenate([
                np.random.normal(0.3, 0.1, n // 2),
                np.random.normal(0.8, 0.1, n // 2),
            ]),
            "good": np.concatenate([
                np.random.normal(0.4, 0.15, n // 2),
                np.random.normal(0.7, 0.15, n // 2),
            ]),
            "with_missing": np.concatenate([
                np.random.normal(0.35, 0.12, n // 2),
                np.random.normal(0.75, 0.12, n // 2),
            ]),
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
