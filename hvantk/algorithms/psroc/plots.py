"""
ROC curve visualization functions for PSROC pipeline.

This module provides static (matplotlib) plotting functions for ROC curves
and related metrics visualization. Follows the styling conventions established
in hvantk.algorithms.visualization.base.

Main Functions:
    - plot_roc_curves: Multi-score ROC curve overlay plot
    - plot_roc_curve_single: Single score ROC curve with detailed annotations
    - plot_auc_comparison: Bar chart comparing AUC across scores
    - plot_missingness_summary: Missingness rates visualization
    - plot_collection_heatmap: Cross-panel AUC heatmap for run_collection results
"""

import logging
from typing import Dict, Optional, Union, List, Tuple
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from hvantk.algorithms.psroc.roc import ROCResult, ScoreMissingness

# Try to import visualization utilities
try:
    from hvantk.algorithms.visualization.base import set_default_style, save_figure, get_colors
except ImportError:
    # Fallback implementations if base module not available
    def set_default_style(style: str = "default", **kwargs) -> None:
        plt.style.use("default")

    def save_figure(fig, filename, dpi=300, formats=None, **kwargs) -> None:
        if formats is None:
            formats = ["png"]
        for fmt in formats:
            fig.savefig(f"{filename}.{fmt}", dpi=dpi, bbox_inches="tight")

    def get_colors(n_colors: int, palette: str = "deep", **kwargs) -> List[str]:
        cmap = plt.get_cmap("tab10")
        return [cmap(i % 10) for i in range(n_colors)]


logger = logging.getLogger(__name__)

# PSROC-specific color scheme
PSROC_COLORS = {
    "diagonal": "#888888",  # Gray for random classifier line
    "optimal": "#FFD700",  # Gold for optimal threshold markers
    "good_auc": "#2E8B57",  # Sea Green for AUC >= 0.8
    "medium_auc": "#FF8C00",  # Dark Orange for 0.6 <= AUC < 0.8
    "poor_auc": "#DC143C",  # Crimson for AUC < 0.6
    "included": "#4682B4",  # Steel Blue for included scores
    "excluded": "#DC143C",  # Crimson for excluded scores
}


def _get_auc_color(auc: float) -> str:
    """Get color based on AUC value."""
    if auc >= 0.8:
        return PSROC_COLORS["good_auc"]
    elif auc >= 0.6:
        return PSROC_COLORS["medium_auc"]
    else:
        return PSROC_COLORS["poor_auc"]


def plot_roc_curves(
    results: Dict[str, ROCResult],
    output_path: Optional[Union[str, Path]] = None,
    title: str = "PSROC: Variant Pathogenicity Prediction",
    figsize: Tuple[float, float] = (10, 8),
    show_optimal: bool = True,
    show_diagonal: bool = True,
    show_legend: bool = True,
    show_auc_in_legend: bool = True,
    style: str = "default",
    color_by_auc: bool = False,
    **kwargs,
) -> plt.Figure:
    """
    Plot ROC curves for multiple prediction scores.

    Creates an overlay plot with one ROC curve per score, optionally marking
    optimal operating points and showing AUC values in the legend.

    Parameters
    ----------
    results : Dict[str, ROCResult]
        Dictionary mapping score name to ROCResult objects.
    output_path : str or Path, optional
        Path to save the figure. If None, figure is not saved.
    title : str
        Plot title.
    figsize : tuple
        Figure size (width, height) in inches.
    show_optimal : bool
        Whether to mark optimal threshold points on each curve.
    show_diagonal : bool
        Whether to show the diagonal reference line (random classifier).
    show_legend : bool
        Whether to show the legend.
    show_auc_in_legend : bool
        Whether to include AUC values in legend labels.
    style : str
        Matplotlib style ('default', 'publication', etc.).
    color_by_auc : bool
        If True, color curves by AUC quality; otherwise use distinct colors.
    **kwargs
        Additional keyword arguments passed to save_figure().

    Returns
    -------
    matplotlib.figure.Figure
        The created figure object.

    Examples
    --------
    >>> from hvantk.algorithms.psroc import compute_roc_metrics
    >>> from hvantk.algorithms.psroc.plots import plot_roc_curves
    >>> results = compute_roc_metrics(labels, scores)
    >>> fig = plot_roc_curves(results, output_path="roc_curves.png")
    """
    set_default_style(style)

    if len(results) == 0:
        raise ValueError("No ROC results to plot")

    fig, ax = plt.subplots(figsize=figsize)

    # Get colors for each score
    n_scores = len(results)
    if color_by_auc:
        colors = {name: _get_auc_color(r.auc) for name, r in results.items()}
    else:
        color_list = get_colors(n_scores)
        colors = {name: color_list[i] for i, name in enumerate(results.keys())}

    # Sort results by AUC (descending) for consistent ordering
    sorted_results = sorted(results.items(), key=lambda x: x[1].auc, reverse=True)

    # Plot each ROC curve
    for score_name, roc_result in sorted_results:
        color = colors[score_name]

        # Create label
        if show_auc_in_legend:
            if roc_result.auc_ci_lower is not None:
                label = (
                    f"{score_name} (AUC = {roc_result.auc:.3f} "
                    f"[{roc_result.auc_ci_lower:.3f}\u2013"
                    f"{roc_result.auc_ci_upper:.3f}])"
                )
            else:
                label = f"{score_name} (AUC = {roc_result.auc:.3f})"
        else:
            label = score_name

        # Plot ROC curve
        ax.plot(
            roc_result.fpr,
            roc_result.tpr,
            color=color,
            linewidth=2,
            label=label,
        )

        # Mark optimal threshold point
        if show_optimal:
            # Find the point closest to the optimal threshold
            opt_idx = np.argmin(
                np.abs(roc_result.thresholds - roc_result.optimal_threshold)
            )
            ax.scatter(
                roc_result.fpr[opt_idx],
                roc_result.tpr[opt_idx],
                color=color,
                s=100,
                marker="o",
                edgecolors="black",
                linewidths=1,
                zorder=5,
            )

    # Plot diagonal reference line
    if show_diagonal:
        ax.plot(
            [0, 1],
            [0, 1],
            color=PSROC_COLORS["diagonal"],
            linestyle="--",
            linewidth=1,
            label="Random (AUC = 0.500)",
        )

    # Customize plot
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.05])
    ax.set_xlabel("False Positive Rate (1 - Specificity)", fontsize=12)
    ax.set_ylabel("True Positive Rate (Sensitivity)", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")

    # Add grid
    ax.grid(True, alpha=0.3)

    # Add legend
    if show_legend:
        ax.legend(loc="lower right", fontsize=10)

    plt.tight_layout()

    # Save figure if path provided
    if output_path:
        save_figure(fig, str(output_path), **kwargs)
        plt.close(fig)

    return fig


def plot_roc_curve_single(
    roc_result: ROCResult,
    output_path: Optional[Union[str, Path]] = None,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (10, 8),
    show_optimal: bool = True,
    show_diagonal: bool = True,
    show_stats: bool = True,
    style: str = "default",
    **kwargs,
) -> plt.Figure:
    """
    Plot a single ROC curve with detailed annotations.

    Provides a detailed view of a single score's ROC curve including
    optimal threshold annotation, statistics box, and AUC shading.

    Parameters
    ----------
    roc_result : ROCResult
        ROC analysis result for a single score.
    output_path : str or Path, optional
        Path to save the figure.
    title : str, optional
        Plot title. If None, uses score name.
    figsize : tuple
        Figure size (width, height) in inches.
    show_optimal : bool
        Whether to mark and annotate the optimal threshold.
    show_diagonal : bool
        Whether to show the diagonal reference line.
    show_stats : bool
        Whether to show a statistics box.
    style : str
        Matplotlib style.
    **kwargs
        Additional keyword arguments passed to save_figure().

    Returns
    -------
    matplotlib.figure.Figure
        The created figure object.
    """
    set_default_style(style)

    if title is None:
        title = f"ROC Curve: {roc_result.score_name}"

    fig, ax = plt.subplots(figsize=figsize)

    # Get color based on AUC
    curve_color = _get_auc_color(roc_result.auc)

    # Fill area under curve
    ax.fill_between(
        roc_result.fpr,
        roc_result.tpr,
        alpha=0.2,
        color=curve_color,
    )

    # Plot ROC curve
    ax.plot(
        roc_result.fpr,
        roc_result.tpr,
        color=curve_color,
        linewidth=2.5,
        label=f"AUC = {roc_result.auc:.3f}",
    )

    # Plot diagonal reference line
    if show_diagonal:
        ax.plot(
            [0, 1],
            [0, 1],
            color=PSROC_COLORS["diagonal"],
            linestyle="--",
            linewidth=1,
            label="Random Classifier",
        )

    # Mark optimal threshold point
    if show_optimal:
        opt_idx = np.argmin(
            np.abs(roc_result.thresholds - roc_result.optimal_threshold)
        )
        opt_fpr = roc_result.fpr[opt_idx]
        opt_tpr = roc_result.tpr[opt_idx]

        ax.scatter(
            opt_fpr,
            opt_tpr,
            color=PSROC_COLORS["optimal"],
            s=150,
            marker="*",
            edgecolors="black",
            linewidths=1,
            zorder=5,
            label=f"Optimal (t={roc_result.optimal_threshold:.3f})",
        )

        # Add annotation for optimal point
        ax.annotate(
            f"Sens: {roc_result.sensitivity_at_optimal:.2f}\n"
            f"Spec: {roc_result.specificity_at_optimal:.2f}",
            xy=(opt_fpr, opt_tpr),
            xytext=(opt_fpr + 0.15, opt_tpr - 0.15),
            fontsize=9,
            arrowprops=dict(arrowstyle="->", color="black", lw=1),
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
        )

    # Add statistics box
    if show_stats:
        auc_line = f"AUC: {roc_result.auc:.3f}"
        if roc_result.auc_ci_lower is not None:
            auc_line += (
                f" [{roc_result.auc_ci_lower:.3f}"
                f"\u2013{roc_result.auc_ci_upper:.3f}]"
            )
        stats_text = (
            f"Score: {roc_result.score_name}\n"
            f"{auc_line}\n"
            f"Optimal Threshold: {roc_result.optimal_threshold:.3f}\n"
            f"Sensitivity: {roc_result.sensitivity_at_optimal:.3f}\n"
            f"Specificity: {roc_result.specificity_at_optimal:.3f}\n"
            f"Variants: {roc_result.n_variants_used}\n"
            f"Missingness: {roc_result.missingness.missingness_rate:.1%}"
        )
        ax.text(
            0.98,
            0.02,
            stats_text,
            transform=ax.transAxes,
            verticalalignment="bottom",
            horizontalalignment="right",
            fontsize=9,
            fontfamily="monospace",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.9),
        )

    # Customize plot
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.05])
    ax.set_xlabel("False Positive Rate (1 - Specificity)", fontsize=12)
    ax.set_ylabel("True Positive Rate (Sensitivity)", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower right", fontsize=10)

    plt.tight_layout()

    if output_path:
        save_figure(fig, str(output_path), **kwargs)
        plt.close(fig)

    return fig


def plot_auc_comparison(
    results: Dict[str, ROCResult],
    output_path: Optional[Union[str, Path]] = None,
    title: str = "AUC Comparison Across Prediction Scores",
    figsize: Tuple[float, float] = (12, 6),
    horizontal: bool = True,
    show_values: bool = True,
    color_by_auc: bool = True,
    style: str = "default",
    **kwargs,
) -> plt.Figure:
    """
    Create a bar chart comparing AUC values across prediction scores.

    Parameters
    ----------
    results : Dict[str, ROCResult]
        Dictionary mapping score name to ROCResult objects.
    output_path : str or Path, optional
        Path to save the figure.
    title : str
        Plot title.
    figsize : tuple
        Figure size (width, height) in inches.
    horizontal : bool
        If True, creates horizontal bars; otherwise vertical.
    show_values : bool
        Whether to annotate bars with AUC values.
    color_by_auc : bool
        If True, color bars by AUC quality.
    style : str
        Matplotlib style.
    **kwargs
        Additional keyword arguments passed to save_figure().

    Returns
    -------
    matplotlib.figure.Figure
        The created figure object.
    """
    set_default_style(style)

    if len(results) == 0:
        raise ValueError("No ROC results to plot")

    # Sort by AUC descending
    sorted_items = sorted(results.items(), key=lambda x: x[1].auc, reverse=True)
    names = [item[0] for item in sorted_items]
    aucs = [item[1].auc for item in sorted_items]

    # Get colors
    if color_by_auc:
        colors = [_get_auc_color(auc) for auc in aucs]
    else:
        colors = get_colors(len(aucs))

    fig, ax = plt.subplots(figsize=figsize)

    if horizontal:
        y_pos = np.arange(len(names))

        for i, (name, roc_result) in enumerate(sorted_items):
            ci_lower = roc_result.auc_ci_lower
            ci_upper = roc_result.auc_ci_upper
            xerr = (
                [[roc_result.auc - ci_lower], [ci_upper - roc_result.auc]]
                if ci_lower is not None
                else None
            )
            ax.errorbar(
                roc_result.auc,
                y_pos[i],
                xerr=xerr,
                fmt="o",
                color=colors[i],
                markersize=8,
                capsize=4,
                capthick=1.5,
                elinewidth=2,
                zorder=5,
            )

        ax.set_yticks(y_pos)
        ax.set_yticklabels(names)
        ax.set_xlabel("AUC (Area Under ROC Curve)", fontsize=12)
        ax.set_xlim([0, 1.0])

        # Add value annotations
        if show_values:
            for i, (name, roc_result) in enumerate(sorted_items):
                ci_lower = roc_result.auc_ci_lower
                ci_upper = roc_result.auc_ci_upper
                if ci_lower is not None:
                    text = (
                        f"{roc_result.auc:.3f} " f"[{ci_lower:.3f}\u2013{ci_upper:.3f}]"
                    )
                    x_text = ci_upper + 0.02
                else:
                    text = f"{roc_result.auc:.3f}"
                    x_text = roc_result.auc + 0.02
                ax.text(x_text, y_pos[i], text, va="center", fontsize=10)

        # Add reference lines
        ax.axvline(
            0.5,
            color=PSROC_COLORS["diagonal"],
            linestyle="--",
            alpha=0.5,
            label="Random",
        )
        ax.axvline(
            0.8,
            color=PSROC_COLORS["good_auc"],
            linestyle="--",
            alpha=0.5,
            label="Good (0.8)",
        )

    else:
        x_pos = np.arange(len(names))

        for i, (name, roc_result) in enumerate(sorted_items):
            ci_lower = roc_result.auc_ci_lower
            ci_upper = roc_result.auc_ci_upper
            yerr = (
                [[roc_result.auc - ci_lower], [ci_upper - roc_result.auc]]
                if ci_lower is not None
                else None
            )
            ax.errorbar(
                x_pos[i],
                roc_result.auc,
                yerr=yerr,
                fmt="o",
                color=colors[i],
                markersize=8,
                capsize=4,
                capthick=1.5,
                elinewidth=2,
                zorder=5,
            )

        ax.set_xticks(x_pos)
        ax.set_xticklabels(names, rotation=45, ha="right")
        ax.set_ylabel("AUC (Area Under ROC Curve)", fontsize=12)
        ax.set_ylim([0, 1.0])

        # Add value annotations
        if show_values:
            for i, (name, roc_result) in enumerate(sorted_items):
                ci_lower = roc_result.auc_ci_lower
                ci_upper = roc_result.auc_ci_upper
                if ci_lower is not None:
                    text = (
                        f"{roc_result.auc:.3f}\n"
                        f"[{ci_lower:.3f}\u2013{ci_upper:.3f}]"
                    )
                    y_text = ci_upper + 0.02
                else:
                    text = f"{roc_result.auc:.3f}"
                    y_text = roc_result.auc + 0.02
                ax.text(x_pos[i], y_text, text, ha="center", fontsize=10)

        # Add reference lines
        ax.axhline(
            0.5,
            color=PSROC_COLORS["diagonal"],
            linestyle="--",
            alpha=0.5,
            label="Random",
        )
        ax.axhline(
            0.8,
            color=PSROC_COLORS["good_auc"],
            linestyle="--",
            alpha=0.5,
            label="Good (0.8)",
        )

    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.legend(loc="lower right" if horizontal else "upper right", fontsize=9)
    ax.grid(True, alpha=0.3, axis="x" if horizontal else "y")

    plt.tight_layout()

    if output_path:
        save_figure(fig, str(output_path), **kwargs)
        plt.close(fig)

    return fig


def plot_missingness_summary(
    missingness: Dict[str, ScoreMissingness],
    output_path: Optional[Union[str, Path]] = None,
    title: str = "Score Missingness Summary",
    figsize: Tuple[float, float] = (12, 6),
    max_missingness_threshold: Optional[float] = None,
    style: str = "default",
    **kwargs,
) -> plt.Figure:
    """
    Create a bar chart showing missingness rates for each score.

    Parameters
    ----------
    missingness : Dict[str, ScoreMissingness]
        Dictionary mapping score name to ScoreMissingness objects.
    output_path : str or Path, optional
        Path to save the figure.
    title : str
        Plot title.
    figsize : tuple
        Figure size (width, height) in inches.
    max_missingness_threshold : float, optional
        If provided, draws a threshold line.
    style : str
        Matplotlib style.
    **kwargs
        Additional keyword arguments passed to save_figure().

    Returns
    -------
    matplotlib.figure.Figure
        The created figure object.
    """
    set_default_style(style)

    if len(missingness) == 0:
        raise ValueError("No missingness data to plot")

    # Sort by missingness rate
    sorted_items = sorted(
        missingness.items(), key=lambda x: x[1].missingness_rate, reverse=True
    )
    names = [item[0] for item in sorted_items]
    rates = [item[1].missingness_rate * 100 for item in sorted_items]  # Convert to %
    included = [item[1].included_in_analysis for item in sorted_items]

    # Color by inclusion status
    colors = [
        PSROC_COLORS["included"] if inc else PSROC_COLORS["excluded"]
        for inc in included
    ]

    fig, ax = plt.subplots(figsize=figsize)

    y_pos = np.arange(len(names))
    bars = ax.barh(y_pos, rates, color=colors, edgecolor="black", linewidth=0.5)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(names)
    ax.set_xlabel("Missingness Rate (%)", fontsize=12)
    max_rate = max(rates) if rates else 0
    ax.set_xlim([0, max(max_rate * 1.1, 1)])

    # Add value annotations
    for bar, rate, inc in zip(bars, rates, included):
        status = "included" if inc else "EXCLUDED"
        ax.text(
            bar.get_width() + 1,
            bar.get_y() + bar.get_height() / 2,
            f"{rate:.1f}% ({status})",
            va="center",
            fontsize=9,
        )

    # Add threshold line
    if max_missingness_threshold is not None:
        threshold_pct = max_missingness_threshold * 100
        ax.axvline(
            threshold_pct,
            color=PSROC_COLORS["excluded"],
            linestyle="--",
            linewidth=2,
            label=f"Threshold ({threshold_pct:.0f}%)",
        )

    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.grid(True, alpha=0.3, axis="x")

    # Add legend for colors (combine with threshold if present)
    from matplotlib.patches import Patch

    legend_elements = [
        Patch(facecolor=PSROC_COLORS["included"], edgecolor="black", label="Included"),
        Patch(facecolor=PSROC_COLORS["excluded"], edgecolor="black", label="Excluded"),
    ]

    # Include threshold line in legend if present
    if max_missingness_threshold is not None:
        from matplotlib.lines import Line2D

        threshold_pct = max_missingness_threshold * 100
        legend_elements.append(
            Line2D(
                [0],
                [0],
                color=PSROC_COLORS["excluded"],
                linestyle="--",
                linewidth=2,
                label=f"Threshold ({threshold_pct:.0f}%)",
            )
        )

    ax.legend(handles=legend_elements, loc="lower right", fontsize=10)

    try:
        plt.tight_layout()
    except ValueError:
        pass  # save_figure uses bbox_inches="tight" as fallback

    if output_path:
        save_figure(fig, str(output_path), **kwargs)
        plt.close(fig)

    return fig


def plot_psroc_summary_dashboard(
    results: Dict[str, ROCResult],
    missingness: Dict[str, ScoreMissingness],
    output_path: Optional[Union[str, Path]] = None,
    title: str = "PSROC Analysis Summary",
    figsize: Tuple[float, float] = (16, 12),
    max_missingness_threshold: Optional[float] = None,
    n_genes: int = 0,
    n_pathogenic: int = 0,
    n_benign: int = 0,
    style: str = "default",
    **kwargs,
) -> plt.Figure:
    """
    Create a comprehensive dashboard combining ROC curves, AUC comparison,
    and missingness summary.

    Parameters
    ----------
    results : Dict[str, ROCResult]
        Dictionary mapping score name to ROCResult objects.
    missingness : Dict[str, ScoreMissingness]
        Dictionary mapping score name to ScoreMissingness objects.
    output_path : str or Path, optional
        Path to save the figure.
    title : str
        Dashboard title.
    figsize : tuple
        Figure size (width, height) in inches.
    max_missingness_threshold : float, optional
        If provided, shows threshold line on missingness plot.
    n_genes : int
        Number of genes in the panel.
    n_pathogenic : int
        Number of pathogenic variants used.
    n_benign : int
        Number of benign variants used.
    style : str
        Matplotlib style.
    **kwargs
        Additional keyword arguments passed to save_figure().

    Returns
    -------
    matplotlib.figure.Figure
        The created figure object.
    """
    set_default_style(style)

    fig = plt.figure(figsize=figsize)
    gs = GridSpec(2, 2, figure=fig, hspace=0.3, wspace=0.3)

    # Top left: ROC curves
    ax1 = fig.add_subplot(gs[0, 0])
    if results:
        n_scores = len(results)
        color_list = get_colors(n_scores)
        colors = {name: color_list[i] for i, name in enumerate(results.keys())}

        sorted_results = sorted(results.items(), key=lambda x: x[1].auc, reverse=True)
        for score_name, roc_result in sorted_results:
            if roc_result.auc_ci_lower is not None:
                legend_label = (
                    f"{score_name} ({roc_result.auc:.3f} "
                    f"[{roc_result.auc_ci_lower:.3f}\u2013"
                    f"{roc_result.auc_ci_upper:.3f}])"
                )
            else:
                legend_label = f"{score_name} ({roc_result.auc:.3f})"
            ax1.plot(
                roc_result.fpr,
                roc_result.tpr,
                color=colors[score_name],
                linewidth=2,
                label=legend_label,
            )

        ax1.plot(
            [0, 1], [0, 1], color=PSROC_COLORS["diagonal"], linestyle="--", linewidth=1
        )
        ax1.set_xlim([0, 1])
        ax1.set_ylim([0, 1.05])
        ax1.set_xlabel("False Positive Rate")
        ax1.set_ylabel("True Positive Rate")
        ax1.set_title("ROC Curves")
        ax1.legend(loc="lower right", fontsize=8)
        ax1.grid(True, alpha=0.3)
    else:
        ax1.text(0.5, 0.5, "No ROC Results\nAvailable", ha="center", va="center")
        ax1.set_title("ROC Curves")

    # Top right: AUC comparison (dot + CI)
    ax2 = fig.add_subplot(gs[0, 1])
    if results:
        sorted_items = sorted(results.items(), key=lambda x: x[1].auc, reverse=True)
        names = [item[0] for item in sorted_items]
        colors_bar = [_get_auc_color(item[1].auc) for item in sorted_items]

        y_pos = np.arange(len(names))
        for i, (name, roc) in enumerate(sorted_items):
            xerr = (
                [[roc.auc - roc.auc_ci_lower], [roc.auc_ci_upper - roc.auc]]
                if roc.auc_ci_lower is not None
                else None
            )
            ax2.errorbar(
                roc.auc,
                y_pos[i],
                xerr=xerr,
                fmt="o",
                color=colors_bar[i],
                markersize=6,
                capsize=3,
                capthick=1,
                elinewidth=1.5,
                zorder=5,
            )
        ax2.set_yticks(y_pos)
        ax2.set_yticklabels(names, fontsize=9)
        ax2.set_xlabel("AUC")
        ax2.set_xlim([0, 1])
        ax2.axvline(0.5, color=PSROC_COLORS["diagonal"], linestyle="--", alpha=0.5)
        ax2.axvline(0.8, color=PSROC_COLORS["good_auc"], linestyle="--", alpha=0.5)
        ax2.set_title("AUC Comparison")
        ax2.grid(True, alpha=0.3, axis="x")
    else:
        ax2.text(0.5, 0.5, "No AUC Data\nAvailable", ha="center", va="center")
        ax2.set_title("AUC Comparison")

    # Bottom left: Missingness summary
    ax3 = fig.add_subplot(gs[1, 0])
    if missingness:
        sorted_miss = sorted(
            missingness.items(), key=lambda x: x[1].missingness_rate, reverse=True
        )
        names = [item[0] for item in sorted_miss]
        rates = [item[1].missingness_rate * 100 for item in sorted_miss]
        included = [item[1].included_in_analysis for item in sorted_miss]
        colors_miss = [
            PSROC_COLORS["included"] if inc else PSROC_COLORS["excluded"]
            for inc in included
        ]

        y_pos = np.arange(len(names))
        ax3.barh(y_pos, rates, color=colors_miss, edgecolor="black", linewidth=0.5)
        ax3.set_yticks(y_pos)
        ax3.set_yticklabels(names, fontsize=9)
        ax3.set_xlabel("Missingness (%)")

        if max_missingness_threshold is not None:
            ax3.axvline(
                max_missingness_threshold * 100,
                color=PSROC_COLORS["excluded"],
                linestyle="--",
                linewidth=2,
            )

        ax3.set_title("Score Missingness")
        ax3.grid(True, alpha=0.3, axis="x")
    else:
        ax3.text(0.5, 0.5, "No Missingness Data\nAvailable", ha="center", va="center")
        ax3.set_title("Score Missingness")

    # Bottom right: Summary statistics
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.axis("off")

    summary_text = "PSROC Analysis Summary\n" + "=" * 40 + "\n\n"

    # Panel and variant info
    if n_genes > 0:
        summary_text += f"GENE PANEL: {n_genes} genes\n"
    if n_pathogenic > 0 or n_benign > 0:
        total_labeled = n_pathogenic + n_benign
        bp_ratio = f"{n_benign / n_pathogenic:.1f}:1" if n_pathogenic > 0 else "N/A"
        summary_text += f"VARIANTS: {total_labeled} (P={n_pathogenic}, B={n_benign})\n"
        summary_text += f"B/P RATIO: {bp_ratio}\n"
    summary_text += "\n"

    if results:
        best_score = max(results.items(), key=lambda x: x[1].auc)
        worst_score = min(results.items(), key=lambda x: x[1].auc)

        summary_text += f"SCORES ANALYZED: {len(results)}\n\n"
        summary_text += "Best Score:\n"
        summary_text += f"  {best_score[0]}\n"
        best_auc_line = f"  AUC: {best_score[1].auc:.3f}"
        if best_score[1].auc_ci_lower is not None:
            best_auc_line += (
                f" [{best_score[1].auc_ci_lower:.3f}"
                f"\u2013{best_score[1].auc_ci_upper:.3f}]"
            )
        summary_text += best_auc_line + "\n"
        summary_text += f"  Optimal Threshold: {best_score[1].optimal_threshold:.3f}\n"
        summary_text += f"  Sensitivity: {best_score[1].sensitivity_at_optimal:.3f}\n"
        summary_text += f"  Specificity: {best_score[1].specificity_at_optimal:.3f}\n\n"

        if len(results) > 1:
            summary_text += "Worst Score:\n"
            summary_text += f"  {worst_score[0]}\n"
            summary_text += f"  AUC: {worst_score[1].auc:.3f}\n\n"

    if missingness:
        n_included = sum(1 for m in missingness.values() if m.included_in_analysis)
        n_excluded = len(missingness) - n_included
        summary_text += "MISSINGNESS:\n"
        summary_text += f"  Scores Included: {n_included}\n"
        summary_text += f"  Scores Excluded: {n_excluded}\n"
        if max_missingness_threshold is not None:
            summary_text += f"  Threshold: {max_missingness_threshold:.0%}\n"

    ax4.text(
        0.05,
        0.95,
        summary_text,
        transform=ax4.transAxes,
        verticalalignment="top",
        fontfamily="monospace",
        fontsize=10,
        bbox=dict(boxstyle="round", facecolor="lightgray", alpha=0.9),
    )

    fig.suptitle(title, fontsize=16, fontweight="bold")

    if output_path:
        save_figure(fig, str(output_path), **kwargs)
        plt.close(fig)

    return fig


def plot_collection_heatmap(
    collection_metrics: Dict[str, Dict[str, ROCResult]],
    output_path: Optional[Union[str, Path]] = None,
    title: str = "AUC Across Gene Set Panels",
    figsize: Optional[Tuple[float, float]] = None,
    show_ci: bool = True,
    show_values: bool = True,
    sort_scores_by: str = "mean_auc",
    cmap: str = "RdYlGn",
    vmin: float = 0.5,
    vmax: float = 1.0,
    style: str = "default",
    **kwargs,
) -> plt.Figure:
    """
    Create a heatmap of AUC values across gene set panels and scores.

    Provides a single summary view of ``run_collection()`` results with
    rows = prediction scores, columns = gene set groups, and cells colored
    by AUC.  Optionally annotates each cell with the AUC value and 95% CI.

    Parameters
    ----------
    collection_metrics : Dict[str, Dict[str, ROCResult]]
        Outer key = group name, inner key = score name, value = ROCResult.
        Typically built as ``{name: result.metrics for name, result in
        run_collection_results.items()}``.
    output_path : str or Path, optional
        Path to save the figure (without extension).
    title : str
        Plot title.
    figsize : tuple, optional
        Figure size (width, height) in inches.  If None, auto-sized based
        on the number of groups and scores.
    show_ci : bool
        Whether to show 95% CI below the AUC value in each cell.
    show_values : bool
        Whether to annotate cells with AUC values at all.
    sort_scores_by : str
        How to order rows.  ``"mean_auc"`` (default) sorts by mean AUC
        across panels (descending).  ``"name"`` sorts alphabetically.
    cmap : str
        Matplotlib colormap name.
    vmin : float
        Minimum value for color scale.
    vmax : float
        Maximum value for color scale.
    style : str
        Matplotlib style.
    **kwargs
        Additional keyword arguments passed to ``save_figure()``.

    Returns
    -------
    matplotlib.figure.Figure
        The created figure object.

    Raises
    ------
    ValueError
        If ``collection_metrics`` is empty.
    """
    set_default_style(style)

    if not collection_metrics:
        raise ValueError("No collection metrics to plot")

    # Collect all group names and score names
    group_names = sorted(collection_metrics.keys())
    score_names_set: set = set()
    for metrics in collection_metrics.values():
        score_names_set.update(metrics.keys())

    if not score_names_set:
        raise ValueError("No scores found in collection metrics")

    # Sort scores
    score_names = sorted(score_names_set)
    if sort_scores_by == "mean_auc":

        def _mean_auc(score: str) -> float:
            aucs = [
                collection_metrics[g][score].auc
                for g in group_names
                if score in collection_metrics[g]
            ]
            return np.mean(aucs) if aucs else 0.0

        score_names = sorted(score_names, key=_mean_auc, reverse=True)

    n_scores = len(score_names)
    n_groups = len(group_names)

    # Build AUC matrix (scores x groups), NaN where missing
    auc_matrix = np.full((n_scores, n_groups), np.nan)
    for j, group in enumerate(group_names):
        for i, score in enumerate(score_names):
            if score in collection_metrics[group]:
                auc_matrix[i, j] = collection_metrics[group][score].auc

    # Auto-size figure if not provided
    if figsize is None:
        width = max(8, 2.0 + n_groups * 1.8)
        height = max(4, 1.5 + n_scores * 0.7)
        figsize = (width, height)

    fig, ax = plt.subplots(figsize=figsize)

    # Plot heatmap
    masked = np.ma.array(auc_matrix, mask=np.isnan(auc_matrix))
    colormap = plt.get_cmap(cmap)
    colormap.set_bad(color="#e0e0e0")

    im = ax.imshow(
        masked,
        cmap=colormap,
        vmin=vmin,
        vmax=vmax,
        aspect="auto",
        interpolation="nearest",
    )

    # Add colorbar
    cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label("AUC", fontsize=11)

    # Set tick labels
    ax.set_xticks(np.arange(n_groups))
    ax.set_yticks(np.arange(n_scores))
    ax.set_xticklabels(group_names, rotation=45, ha="right", fontsize=9)
    ax.set_yticklabels(score_names, fontsize=10)

    # Annotate cells
    if show_values:
        for i, score in enumerate(score_names):
            for j, group in enumerate(group_names):
                if score not in collection_metrics[group]:
                    ax.text(
                        j,
                        i,
                        "n/a",
                        ha="center",
                        va="center",
                        fontsize=8,
                        color="#999999",
                        style="italic",
                    )
                    continue

                roc = collection_metrics[group][score]
                auc_val = roc.auc

                # Choose text color for readability against cell color
                text_color = "white" if auc_val < (vmin + vmax) / 2 else "black"

                if show_ci and roc.auc_ci_lower is not None:
                    cell_text = (
                        f"{auc_val:.3f}\n"
                        f"[{roc.auc_ci_lower:.3f}\u2013{roc.auc_ci_upper:.3f}]"
                    )
                    ax.text(
                        j,
                        i,
                        cell_text,
                        ha="center",
                        va="center",
                        fontsize=8,
                        color=text_color,
                        fontweight="bold",
                        linespacing=1.4,
                    )
                else:
                    ax.text(
                        j,
                        i,
                        f"{auc_val:.3f}",
                        ha="center",
                        va="center",
                        fontsize=9,
                        color=text_color,
                        fontweight="bold",
                    )

    # Grid lines between cells
    ax.set_xticks(np.arange(n_groups + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(n_scores + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="minor", bottom=False, left=False)

    ax.set_title(title, fontsize=14, fontweight="bold", pad=12)
    ax.set_xlabel("Gene Set Panel", fontsize=11)

    plt.tight_layout()

    if output_path:
        save_figure(fig, str(output_path), **kwargs)
        plt.close(fig)

    return fig
