"""
Quality Control (QC) plotting functions for hvantk.

This module provides comprehensive visualization functions for genomic variant QC data,
including sample-level and variant-level quality metrics visualization using matplotlib.

The module provides static plots suitable for:
- Publication-quality figures
- Batch processing and automation
- PDF/PNG export for reports
- Consistent styling across analyses

Main Functions:
    Sample QC:
    - plot_sample_call_rate_distribution: Sample call rate histograms with thresholds
    - plot_sample_titv_distribution: Ti/Tv ratio distributions with expected ranges
    - plot_sample_depth_distribution: Depth statistics distributions
    - plot_sample_qc_overview: Multi-panel sample QC dashboard

    Variant QC:
    - plot_variant_call_rate_distribution: Variant call rate histograms with thresholds
    - plot_allele_frequency_spectrum: Allele frequency distributions with MAF lines
    - plot_hwe_pvalues: Hardy-Weinberg equilibrium p-value plots
    - plot_variant_qc_overview: Multi-panel variant QC dashboard

    Combined:
    - plot_qc_summary_dashboard: Comprehensive 12-panel QC dashboard

Example:
    >>> from hvantk.hgc import compute_full_qc
    >>> from hvantk.visualization.qc_plots import plot_qc_summary_dashboard
    >>> qc_results = compute_full_qc(mt)
    >>> fig = plot_qc_summary_dashboard(qc_results)
    >>> plt.show()
"""

import logging
from typing import Optional, Union, Dict, Tuple
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Try to import optional dependencies
try:
    import seaborn as sns

    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False

try:
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots
    import plotly.offline as pyo

    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False

from .base import set_default_style, save_figure

logger = logging.getLogger(__name__)

# Default color schemes for QC plots
QC_COLORS = {
    "pass": "#2E8B57",  # Sea Green
    "warn": "#FF8C00",  # Dark Orange
    "fail": "#DC143C",  # Crimson
    "neutral": "#4682B4",  # Steel Blue
    "highlight": "#FFD700",  # Gold
    "background": "#F5F5F5",  # White Smoke
}

QC_THRESHOLDS = {
    "sample_call_rate": {"good": 0.95, "acceptable": 0.85},
    "variant_call_rate": {"good": 0.95, "acceptable": 0.85},
    "ti_tv_ratio": {
        "good_min": 1.8,
        "good_max": 2.2,
        "acceptable_min": 1.5,
        "acceptable_max": 2.5,
    },
    "hwe_pvalue": {"fail": 1e-6, "warn": 1e-4},
}


def _prepare_sample_qc_data(sample_df: pd.DataFrame) -> pd.DataFrame:
    """
    Prepare sample QC data for plotting by flattening nested columns.

    Parameters
    ----------
    sample_df : pd.DataFrame
        Sample QC metrics DataFrame with potentially nested columns

    Returns
    -------
    pd.DataFrame
        Flattened DataFrame with accessible column names
    """
    df = sample_df.copy()

    # Flatten nested column names, handling duplicates intelligently
    new_columns = []
    seen_names = {}

    for col in df.columns:
        if isinstance(col, str) and "." in col:
            # Extract the metric name after the last dot
            base_name = col.split(".")[-1]

            # Handle duplicates by keeping context
            if base_name in seen_names:
                # For duplicates, use more context (e.g., 'dp_mean', 'gq_mean')
                parts = col.split(".")
                if len(parts) >= 3:
                    context_name = f"{parts[-2]}_{base_name}"
                else:
                    context_name = f"{seen_names[base_name]}_{base_name}"

                # Update the previous occurrence
                prev_idx = seen_names[base_name]
                new_columns[prev_idx] = (
                    f"{df.columns[prev_idx].split('.')[-2]}_{base_name}"
                )

                new_columns.append(context_name)
            else:
                seen_names[base_name] = len(new_columns)
                new_columns.append(base_name)
        else:
            new_columns.append(col)

    df.columns = new_columns
    return df


def _prepare_variant_qc_data(variant_df: pd.DataFrame) -> pd.DataFrame:
    """
    Prepare variant QC data for plotting by flattening nested columns.

    Parameters
    ----------
    variant_df : pd.DataFrame
        Variant QC metrics DataFrame with potentially nested columns

    Returns
    -------
    pd.DataFrame
        Flattened DataFrame with accessible column names
    """
    df = variant_df.copy()

    # Flatten nested column names
    new_columns = []
    for col in df.columns:
        if isinstance(col, str) and "." in col:
            new_name = col.split(".")[-1]
            new_columns.append(new_name)
        else:
            new_columns.append(col)

    df.columns = new_columns

    # Handle array columns (AC, AF) by extracting alternate allele values
    if "AC" in df.columns:
        # Extract alternate allele count (index 1)
        try:
            df["AC_alt"] = df["AC"].apply(
                lambda x: (
                    x[1] if isinstance(x, (list, np.ndarray)) and len(x) > 1 else x
                )
            )
        except (TypeError, IndexError):
            logger.warning("Could not extract alternate allele count from AC column")

    if "AF" in df.columns:
        # Extract alternate allele frequency (index 1)
        try:
            df["AF_alt"] = df["AF"].apply(
                lambda x: (
                    x[1] if isinstance(x, (list, np.ndarray)) and len(x) > 1 else x
                )
            )
        except (TypeError, IndexError):
            logger.warning(
                "Could not extract alternate allele frequency from AF column"
            )

    return df


def _add_threshold_lines(
    ax, metric: str, thresholds: Optional[Dict] = None, orientation: str = "vertical"
) -> None:
    """
    Add threshold lines to plots.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Plot axes
    metric : str
        QC metric name
    thresholds : dict, optional
        Custom thresholds dictionary
    orientation : str
        'vertical' or 'horizontal' threshold lines
    """
    if thresholds is None:
        thresholds = QC_THRESHOLDS.get(metric, {})

    if not thresholds:
        return

    line_func = ax.axvline if orientation == "vertical" else ax.axhline

    if "good" in thresholds:
        line_func(
            thresholds["good"],
            color=QC_COLORS["pass"],
            linestyle="--",
            alpha=0.7,
            label="Good threshold",
        )

    if "acceptable" in thresholds:
        line_func(
            thresholds["acceptable"],
            color=QC_COLORS["warn"],
            linestyle="--",
            alpha=0.7,
            label="Acceptable threshold",
        )

    if "fail" in thresholds:
        line_func(
            thresholds["fail"],
            color=QC_COLORS["fail"],
            linestyle="--",
            alpha=0.7,
            label="Fail threshold",
        )

    # Handle range thresholds (e.g., Ti/Tv ratio)
    if "good_min" in thresholds and "good_max" in thresholds:
        line_func(
            thresholds["good_min"], color=QC_COLORS["pass"], linestyle="--", alpha=0.7
        )
        line_func(
            thresholds["good_max"],
            color=QC_COLORS["pass"],
            linestyle="--",
            alpha=0.7,
            label="Good range",
        )


def plot_sample_call_rate_distribution(
    sample_df: pd.DataFrame,
    thresholds: Optional[Dict] = None,
    bins: int = 50,
    kde: bool = True,
    outliers: bool = True,
    style: str = "default",
    figsize: Tuple[float, float] = (10, 6),
    save_path: Optional[Union[str, Path]] = None,
    **kwargs,
) -> plt.Figure:
    """
    Plot distribution of sample call rates.

    Parameters
    ----------
    sample_df : pd.DataFrame
        Sample QC metrics DataFrame
    thresholds : dict, optional
        Custom call rate thresholds
    bins : int
        Number of histogram bins
    kde : bool
        Whether to overlay kernel density estimate
    outliers : bool
        Whether to highlight outliers
    style : str
        Plot style ('default', 'publication', etc.)
    figsize : tuple
        Figure size (width, height)
    save_path : str or Path, optional
        Path to save the figure
    **kwargs
        Additional plotting arguments

    Returns
    -------
    matplotlib.figure.Figure
        The created figure
    """
    set_default_style(style)

    # Prepare data
    df = _prepare_sample_qc_data(sample_df)

    if "call_rate" not in df.columns:
        raise ValueError("Sample DataFrame must contain 'call_rate' column")

    call_rates = df["call_rate"].dropna()

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Plot histogram
    n, bins_array, patches = ax.hist(
        call_rates,
        bins=bins,
        alpha=0.7,
        color=QC_COLORS["neutral"],
        edgecolor="black",
        linewidth=0.5,
    )

    # Add KDE if requested and seaborn is available
    if kde and HAS_SEABORN:
        try:
            sns.kdeplot(
                data=call_rates, ax=ax, color="red", linewidth=2, warn_singular=False
            )
        except Exception as e:
            logger.warning(f"Could not add KDE overlay: {e}")

    # Add threshold lines
    _add_threshold_lines(ax, "sample_call_rate", thresholds)

    # Highlight outliers
    if outliers:
        default_thresholds = QC_THRESHOLDS["sample_call_rate"]
        threshold = (
            thresholds.get("acceptable", default_thresholds.get("acceptable", 0.85))
            if thresholds
            else default_thresholds.get("acceptable", 0.85)
        )

        outlier_mask = call_rates < threshold
        if outlier_mask.any():
            outlier_samples = call_rates[outlier_mask]
            ax.hist(
                outlier_samples,
                bins=bins_array,
                alpha=0.8,
                color=QC_COLORS["fail"],
                label=f"Low call rate (n={len(outlier_samples)})",
            )

    # Customize plot
    ax.set_xlabel("Sample Call Rate", fontsize=12)
    ax.set_ylabel("Number of Samples", fontsize=12)
    ax.set_title(
        f"Sample Call Rate Distribution (n={len(call_rates)})",
        fontsize=14,
        fontweight="bold",
    )

    # Add statistics
    stats_text = f"Mean: {call_rates.mean():.3f}\nStd: {call_rates.std():.3f}\nMin: {call_rates.min():.3f}\nMax: {call_rates.max():.3f}"
    ax.text(
        0.02,
        0.98,
        stats_text,
        transform=ax.transAxes,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
    )

    # Add legend if there are threshold lines or outliers
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(loc="upper right")

    plt.tight_layout()

    if save_path:
        save_figure(fig, save_path, **kwargs)

    return fig


def plot_sample_titv_distribution(
    sample_df: pd.DataFrame,
    thresholds: Optional[Dict] = None,
    bins: int = 50,
    kde: bool = True,
    style: str = "default",
    figsize: Tuple[float, float] = (10, 6),
    save_path: Optional[Union[str, Path]] = None,
    **kwargs,
) -> plt.Figure:
    """
    Plot distribution of sample Ti/Tv ratios.

    Parameters
    ----------
    sample_df : pd.DataFrame
        Sample QC metrics DataFrame
    thresholds : dict, optional
        Custom Ti/Tv ratio thresholds
    bins : int
        Number of histogram bins
    kde : bool
        Whether to overlay kernel density estimate
    style : str
        Plot style
    figsize : tuple
        Figure size (width, height)
    save_path : str or Path, optional
        Path to save the figure
    **kwargs
        Additional plotting arguments

    Returns
    -------
    matplotlib.figure.Figure
        The created figure
    """
    set_default_style(style)

    # Prepare data
    df = _prepare_sample_qc_data(sample_df)

    if "r_ti_tv" not in df.columns:
        raise ValueError("Sample DataFrame must contain 'r_ti_tv' column")

    ti_tv_ratios = df["r_ti_tv"].dropna()

    # Filter extreme outliers for better visualization
    q1, q99 = np.percentile(ti_tv_ratios, [1, 99])
    filtered_ratios = ti_tv_ratios[(ti_tv_ratios >= q1) & (ti_tv_ratios <= q99)]

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Plot histogram
    ax.hist(
        filtered_ratios,
        bins=bins,
        alpha=0.7,
        color=QC_COLORS["neutral"],
        edgecolor="black",
        linewidth=0.5,
    )

    # Add KDE if requested
    if kde and HAS_SEABORN:
        try:
            sns.kdeplot(
                data=filtered_ratios,
                ax=ax,
                color="red",
                linewidth=2,
                warn_singular=False,
            )
        except Exception as e:
            logger.warning(f"Could not add KDE overlay: {e}")

    # Add threshold lines for expected Ti/Tv range
    _add_threshold_lines(ax, "ti_tv_ratio", thresholds)

    # Customize plot
    ax.set_xlabel("Ti/Tv Ratio", fontsize=12)
    ax.set_ylabel("Number of Samples", fontsize=12)
    ax.set_title(
        f"Sample Ti/Tv Ratio Distribution (n={len(ti_tv_ratios)})",
        fontsize=14,
        fontweight="bold",
    )

    # Add statistics
    stats_text = f"Mean: {ti_tv_ratios.mean():.3f}\nStd: {ti_tv_ratios.std():.3f}\nMedian: {ti_tv_ratios.median():.3f}"
    ax.text(
        0.98,
        0.98,
        stats_text,
        transform=ax.transAxes,
        verticalalignment="top",
        horizontalalignment="right",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
    )

    # Add legend if there are threshold lines
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(loc="upper left")

    plt.tight_layout()

    if save_path:
        save_figure(fig, save_path, **kwargs)

    return fig


def plot_sample_depth_distribution(
    sample_df: pd.DataFrame,
    bins: int = 50,
    kde: bool = True,
    log_scale: bool = False,
    style: str = "default",
    figsize: Tuple[float, float] = (10, 6),
    save_path: Optional[Union[str, Path]] = None,
    **kwargs,
) -> plt.Figure:
    """
    Plot distribution of sample depth statistics.

    Parameters
    ----------
    sample_df : pd.DataFrame
        Sample QC metrics DataFrame
    bins : int
        Number of histogram bins
    kde : bool
        Whether to overlay kernel density estimate
    log_scale : bool
        Whether to use log scale for y-axis
    style : str
        Plot style
    figsize : tuple
        Figure size (width, height)
    save_path : str or Path, optional
        Path to save the figure
    **kwargs
        Additional plotting arguments

    Returns
    -------
    matplotlib.figure.Figure
        The created figure
    """
    set_default_style(style)

    # Prepare data
    df = _prepare_sample_qc_data(sample_df)

    # Look for depth mean column
    depth_col = None
    for col in ["dp_stats_mean", "mean", "sample_qc_mean_dp", "dp_mean"]:
        if col in df.columns:
            depth_col = col
            break

    if depth_col is None:
        raise ValueError(
            "Sample DataFrame must contain depth statistics (mean depth column)"
        )

    depths = df[depth_col].dropna()

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Plot histogram
    ax.hist(
        depths,
        bins=bins,
        alpha=0.7,
        color=QC_COLORS["neutral"],
        edgecolor="black",
        linewidth=0.5,
    )

    # Add KDE if requested
    if kde and HAS_SEABORN:
        try:
            sns.kdeplot(
                data=depths, ax=ax, color="red", linewidth=2, warn_singular=False
            )
        except Exception as e:
            logger.warning(f"Could not add KDE overlay: {e}")

    # Set log scale if requested
    if log_scale:
        ax.set_yscale("log")

    # Customize plot
    ax.set_xlabel("Mean Depth", fontsize=12)
    ax.set_ylabel("Number of Samples", fontsize=12)
    ax.set_title(
        f"Sample Depth Distribution (n={len(depths)})", fontsize=14, fontweight="bold"
    )

    # Add statistics
    stats_text = f"Mean: {depths.mean():.2f}\nStd: {depths.std():.2f}\nMedian: {depths.median():.2f}"
    ax.text(
        0.98,
        0.98,
        stats_text,
        transform=ax.transAxes,
        verticalalignment="top",
        horizontalalignment="right",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
    )

    plt.tight_layout()

    if save_path:
        save_figure(fig, save_path, **kwargs)

    return fig


def plot_sample_qc_overview(
    sample_df: pd.DataFrame,
    thresholds: Optional[Dict] = None,
    style: str = "default",
    figsize: Tuple[float, float] = (15, 10),
    save_path: Optional[Union[str, Path]] = None,
    **kwargs,
) -> plt.Figure:
    """
    Create a comprehensive overview of sample QC metrics.

    Parameters
    ----------
    sample_df : pd.DataFrame
        Sample QC metrics DataFrame
    thresholds : dict, optional
        Custom thresholds for different metrics
    style : str
        Plot style
    figsize : tuple
        Figure size (width, height)
    save_path : str or Path, optional
        Path to save the figure
    **kwargs
        Additional plotting arguments

    Returns
    -------
    matplotlib.figure.Figure
        The created figure
    """
    set_default_style(style)

    # Prepare data
    df = _prepare_sample_qc_data(sample_df)

    # Create figure with subplots
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.3)

    # Plot 1: Call rate distribution
    ax1 = fig.add_subplot(gs[0, 0])
    if "call_rate" in df.columns:
        call_rates = df["call_rate"].dropna()
        ax1.hist(
            call_rates,
            bins=30,
            alpha=0.7,
            color=QC_COLORS["neutral"],
            edgecolor="black",
        )
        _add_threshold_lines(ax1, "sample_call_rate", thresholds)
        ax1.set_xlabel("Call Rate")
        ax1.set_ylabel("Count")
        ax1.set_title("Sample Call Rate")
    else:
        ax1.text(
            0.5,
            0.5,
            "Call Rate\nNot Available",
            ha="center",
            va="center",
            transform=ax1.transAxes,
        )

    # Plot 2: Ti/Tv ratio distribution
    ax2 = fig.add_subplot(gs[0, 1])
    if "r_ti_tv" in df.columns:
        ti_tv = df["r_ti_tv"].dropna()
        q1, q99 = np.percentile(ti_tv, [1, 99])
        filtered_ti_tv = ti_tv[(ti_tv >= q1) & (ti_tv <= q99)]
        ax2.hist(
            filtered_ti_tv,
            bins=30,
            alpha=0.7,
            color=QC_COLORS["neutral"],
            edgecolor="black",
        )
        _add_threshold_lines(ax2, "ti_tv_ratio", thresholds)
        ax2.set_xlabel("Ti/Tv Ratio")
        ax2.set_ylabel("Count")
        ax2.set_title("Ti/Tv Ratio")
    else:
        ax2.text(
            0.5,
            0.5,
            "Ti/Tv Ratio\nNot Available",
            ha="center",
            va="center",
            transform=ax2.transAxes,
        )

    # Plot 3: Mean depth distribution
    ax3 = fig.add_subplot(gs[0, 2])
    depth_col = None
    for col in ["dp_stats_mean", "mean", "sample_qc_mean_dp", "dp_mean"]:
        if col in df.columns:
            depth_col = col
            break

    if depth_col:
        depths = df[depth_col].dropna()
        if len(depths) > 0:
            ax3.hist(
                depths,
                bins=30,
                alpha=0.7,
                color=QC_COLORS["neutral"],
                edgecolor="black",
            )
            ax3.set_xlabel("Mean Depth")
            ax3.set_ylabel("Count")
            ax3.set_title("Mean Depth")
        else:
            ax3.text(
                0.5,
                0.5,
                "Mean Depth\nNo Data",
                ha="center",
                va="center",
                transform=ax3.transAxes,
            )
    else:
        ax3.text(
            0.5,
            0.5,
            "Mean Depth\nNot Available",
            ha="center",
            va="center",
            transform=ax3.transAxes,
        )

    # Plot 4: Genotype counts (stacked bar if available)
    ax4 = fig.add_subplot(gs[1, 0])
    genotype_cols = ["n_hom_ref", "n_het", "n_hom_var"]
    available_cols = [col for col in genotype_cols if col in df.columns]

    if available_cols:
        genotype_data = df[available_cols].sum()
        colors = [QC_COLORS["pass"], QC_COLORS["warn"], QC_COLORS["fail"]]
        ax4.pie(
            genotype_data.values,
            labels=genotype_data.index,
            colors=colors[: len(genotype_data)],
            autopct="%1.1f%%",
        )
        ax4.set_title("Genotype Distribution")
    else:
        ax4.text(
            0.5,
            0.5,
            "Genotype Counts\nNot Available",
            ha="center",
            va="center",
            transform=ax4.transAxes,
        )

    # Plot 5: Singleton rate
    ax5 = fig.add_subplot(gs[1, 1])
    if "n_singleton" in df.columns and "n_called" in df.columns:
        singleton_rate = df["n_singleton"] / df["n_called"]
        singleton_rate = singleton_rate.dropna()
        ax5.hist(
            singleton_rate,
            bins=30,
            alpha=0.7,
            color=QC_COLORS["neutral"],
            edgecolor="black",
        )
        ax5.set_xlabel("Singleton Rate")
        ax5.set_ylabel("Count")
        ax5.set_title("Singleton Rate")
    else:
        ax5.text(
            0.5,
            0.5,
            "Singleton Rate\nNot Available",
            ha="center",
            va="center",
            transform=ax5.transAxes,
        )

    # Plot 6: Het/Hom ratio
    ax6 = fig.add_subplot(gs[1, 2])
    if "r_het_hom_var" in df.columns:
        het_hom = df["r_het_hom_var"].dropna()
        ax6.hist(
            het_hom, bins=30, alpha=0.7, color=QC_COLORS["neutral"], edgecolor="black"
        )
        ax6.set_xlabel("Het/Hom Ratio")
        ax6.set_ylabel("Count")
        ax6.set_title("Het/Hom Ratio")
    else:
        ax6.text(
            0.5,
            0.5,
            "Het/Hom Ratio\nNot Available",
            ha="center",
            va="center",
            transform=ax6.transAxes,
        )

    # Add overall title
    fig.suptitle(
        f"Sample QC Overview (n={len(df)} samples)", fontsize=16, fontweight="bold"
    )

    if save_path:
        save_figure(fig, save_path, **kwargs)

    return fig


def plot_variant_call_rate_distribution(
    variant_df: pd.DataFrame,
    thresholds: Optional[Dict] = None,
    bins: int = 50,
    kde: bool = True,
    style: str = "default",
    figsize: Tuple[float, float] = (10, 6),
    save_path: Optional[Union[str, Path]] = None,
    **kwargs,
) -> plt.Figure:
    """
    Plot distribution of variant call rates.

    Parameters
    ----------
    variant_df : pd.DataFrame
        Variant QC metrics DataFrame
    thresholds : dict, optional
        Custom call rate thresholds
    bins : int
        Number of histogram bins
    kde : bool
        Whether to overlay kernel density estimate
    style : str
        Plot style
    figsize : tuple
        Figure size (width, height)
    save_path : str or Path, optional
        Path to save the figure
    **kwargs
        Additional plotting arguments

    Returns
    -------
    matplotlib.figure.Figure
        The created figure
    """
    set_default_style(style)

    # Prepare data
    df = _prepare_variant_qc_data(variant_df)

    if "call_rate" not in df.columns:
        raise ValueError("Variant DataFrame must contain 'call_rate' column")

    call_rates = df["call_rate"].dropna()

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Plot histogram
    ax.hist(
        call_rates,
        bins=bins,
        alpha=0.7,
        color=QC_COLORS["neutral"],
        edgecolor="black",
        linewidth=0.5,
    )

    # Add KDE if requested
    if kde and HAS_SEABORN:
        try:
            sns.kdeplot(
                data=call_rates, ax=ax, color="red", linewidth=2, warn_singular=False
            )
        except Exception as e:
            logger.warning(f"Could not add KDE overlay: {e}")

    # Add threshold lines
    _add_threshold_lines(ax, "variant_call_rate", thresholds)

    # Customize plot
    ax.set_xlabel("Variant Call Rate", fontsize=12)
    ax.set_ylabel("Number of Variants", fontsize=12)
    ax.set_title(
        f"Variant Call Rate Distribution (n={len(call_rates)})",
        fontsize=14,
        fontweight="bold",
    )

    # Add statistics
    stats_text = f"Mean: {call_rates.mean():.3f}\nStd: {call_rates.std():.3f}\nMin: {call_rates.min():.3f}"
    ax.text(
        0.02,
        0.98,
        stats_text,
        transform=ax.transAxes,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
    )

    # Add legend if there are threshold lines
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(loc="upper right")

    plt.tight_layout()

    if save_path:
        save_figure(fig, save_path, **kwargs)

    return fig


def plot_allele_frequency_spectrum(
    variant_df: pd.DataFrame,
    bins: int = 50,
    log_scale: bool = True,
    style: str = "default",
    figsize: Tuple[float, float] = (10, 6),
    save_path: Optional[Union[str, Path]] = None,
    **kwargs,
) -> plt.Figure:
    """
    Plot allele frequency spectrum.

    Parameters
    ----------
    variant_df : pd.DataFrame
        Variant QC metrics DataFrame
    bins : int
        Number of histogram bins
    log_scale : bool
        Whether to use log scale for y-axis
    style : str
        Plot style
    figsize : tuple
        Figure size (width, height)
    save_path : str or Path, optional
        Path to save the figure
    **kwargs
        Additional plotting arguments

    Returns
    -------
    matplotlib.figure.Figure
        The created figure
    """
    set_default_style(style)

    # Prepare data
    df = _prepare_variant_qc_data(variant_df)

    # Look for allele frequency column
    af_col = None
    for col in ["AF_alt", "AF"]:
        if col in df.columns:
            af_col = col
            break

    if af_col is None:
        raise ValueError(
            "Variant DataFrame must contain allele frequency column (AF or AF_alt)"
        )

    # Get allele frequencies, handling both single values and arrays
    if af_col == "AF":
        # Try to extract alternate allele frequency from array
        try:
            afs = df[af_col].apply(
                lambda x: (
                    x[1] if isinstance(x, (list, np.ndarray)) and len(x) > 1 else x
                )
            )
        except:
            afs = df[af_col]
    else:
        afs = df[af_col]

    afs = afs.dropna()

    # Filter to remove monomorphic sites (AF = 0 or 1)
    afs_filtered = afs[(afs > 0) & (afs < 1)]

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Plot histogram
    ax.hist(
        afs_filtered,
        bins=bins,
        alpha=0.7,
        color=QC_COLORS["neutral"],
        edgecolor="black",
        linewidth=0.5,
    )

    # Set log scale if requested
    if log_scale:
        ax.set_yscale("log")

    # Add vertical lines for common frequency categories
    ax.axvline(0.01, color=QC_COLORS["warn"], linestyle="--", alpha=0.7, label="1% MAF")
    ax.axvline(0.05, color=QC_COLORS["pass"], linestyle="--", alpha=0.7, label="5% MAF")

    # Customize plot
    ax.set_xlabel("Allele Frequency", fontsize=12)
    ax.set_ylabel("Number of Variants", fontsize=12)
    ax.set_title(
        f"Allele Frequency Spectrum (n={len(afs_filtered)} polymorphic)",
        fontsize=14,
        fontweight="bold",
    )

    # Add statistics
    stats_text = f"Mean: {afs_filtered.mean():.4f}\nMedian: {afs_filtered.median():.4f}\n<1%: {(afs_filtered < 0.01).sum()}\n<5%: {(afs_filtered < 0.05).sum()}"
    ax.text(
        0.98,
        0.98,
        stats_text,
        transform=ax.transAxes,
        verticalalignment="top",
        horizontalalignment="right",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
    )

    ax.legend(loc="upper right")
    plt.tight_layout()

    if save_path:
        save_figure(fig, save_path, **kwargs)

    return fig


def plot_hwe_pvalues(
    variant_df: pd.DataFrame,
    thresholds: Optional[Dict] = None,
    log_transform: bool = True,
    bins: int = 50,
    style: str = "default",
    figsize: Tuple[float, float] = (10, 6),
    save_path: Optional[Union[str, Path]] = None,
    **kwargs,
) -> plt.Figure:
    """
    Plot Hardy-Weinberg equilibrium p-values.

    Parameters
    ----------
    variant_df : pd.DataFrame
        Variant QC metrics DataFrame
    thresholds : dict, optional
        Custom HWE p-value thresholds
    log_transform : bool
        Whether to plot -log10(p-values)
    bins : int
        Number of histogram bins
    style : str
        Plot style
    figsize : tuple
        Figure size (width, height)
    save_path : str or Path, optional
        Path to save the figure
    **kwargs
        Additional plotting arguments

    Returns
    -------
    matplotlib.figure.Figure
        The created figure
    """
    set_default_style(style)

    # Prepare data
    df = _prepare_variant_qc_data(variant_df)

    if "p_value_hwe" not in df.columns:
        raise ValueError("Variant DataFrame must contain 'p_value_hwe' column")

    pvalues = df["p_value_hwe"].dropna()

    # Remove p-values of 0 (which cause issues with log transform)
    pvalues = pvalues[pvalues > 0]

    if log_transform:
        plot_values = -np.log10(pvalues)
        xlabel = "-log10(HWE p-value)"
        title_suffix = "(-log10 scale)"
    else:
        plot_values = pvalues
        xlabel = "HWE p-value"
        title_suffix = ""

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Plot histogram
    ax.hist(
        plot_values,
        bins=bins,
        alpha=0.7,
        color=QC_COLORS["neutral"],
        edgecolor="black",
        linewidth=0.5,
    )

    # Add threshold lines
    if thresholds is None:
        thresholds = QC_THRESHOLDS.get("hwe_pvalue", {})

    if thresholds:
        for name, threshold in thresholds.items():
            if log_transform:
                line_val = -np.log10(threshold)
            else:
                line_val = threshold

            color = QC_COLORS["fail"] if name == "fail" else QC_COLORS["warn"]
            ax.axvline(
                line_val,
                color=color,
                linestyle="--",
                alpha=0.7,
                label=f"{name.title()}: {threshold:.0e}",
            )

    # Customize plot
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel("Number of Variants", fontsize=12)
    ax.set_title(
        f"Hardy-Weinberg Equilibrium p-values {title_suffix} (n={len(pvalues)})",
        fontsize=14,
        fontweight="bold",
    )

    # Add statistics
    if log_transform:
        stats_text = (
            f"Mean: {plot_values.mean():.2f}\nMedian: {plot_values.median():.2f}"
        )
    else:
        stats_text = (
            f"Mean: {plot_values.mean():.2e}\nMedian: {plot_values.median():.2e}"
        )

    ax.text(
        0.98,
        0.98,
        stats_text,
        transform=ax.transAxes,
        verticalalignment="top",
        horizontalalignment="right",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
    )

    # Add legend if there are threshold lines
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(loc="upper left" if log_transform else "upper right")

    plt.tight_layout()

    if save_path:
        save_figure(fig, save_path, **kwargs)

    return fig


def plot_variant_qc_overview(
    variant_df: pd.DataFrame,
    thresholds: Optional[Dict] = None,
    style: str = "default",
    figsize: Tuple[float, float] = (15, 10),
    save_path: Optional[Union[str, Path]] = None,
    **kwargs,
) -> plt.Figure:
    """
    Create a comprehensive overview of variant QC metrics.

    Parameters
    ----------
    variant_df : pd.DataFrame
        Variant QC metrics DataFrame
    thresholds : dict, optional
        Custom thresholds for different metrics
    style : str
        Plot style
    figsize : tuple
        Figure size (width, height)
    save_path : str or Path, optional
        Path to save the figure
    **kwargs
        Additional plotting arguments

    Returns
    -------
    matplotlib.figure.Figure
        The created figure
    """
    set_default_style(style)

    # Prepare data
    df = _prepare_variant_qc_data(variant_df)

    # Create figure with subplots
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.3)

    # Plot 1: Call rate distribution
    ax1 = fig.add_subplot(gs[0, 0])
    if "call_rate" in df.columns:
        call_rates = df["call_rate"].dropna()
        ax1.hist(
            call_rates,
            bins=30,
            alpha=0.7,
            color=QC_COLORS["neutral"],
            edgecolor="black",
        )
        _add_threshold_lines(ax1, "variant_call_rate", thresholds)
        ax1.set_xlabel("Call Rate")
        ax1.set_ylabel("Count")
        ax1.set_title("Variant Call Rate")
    else:
        ax1.text(
            0.5,
            0.5,
            "Call Rate\nNot Available",
            ha="center",
            va="center",
            transform=ax1.transAxes,
        )

    # Plot 2: Allele frequency spectrum
    ax2 = fig.add_subplot(gs[0, 1])
    af_col = "AF_alt" if "AF_alt" in df.columns else "AF"
    if af_col in df.columns:
        if af_col == "AF":
            try:
                afs = df[af_col].apply(
                    lambda x: (
                        x[1] if isinstance(x, (list, np.ndarray)) and len(x) > 1 else x
                    )
                )
            except:
                afs = df[af_col]
        else:
            afs = df[af_col]

        afs_filtered = afs.dropna()
        afs_filtered = afs_filtered[(afs_filtered > 0) & (afs_filtered < 1)]
        if len(afs_filtered) > 0:
            ax2.hist(
                afs_filtered,
                bins=30,
                alpha=0.7,
                color=QC_COLORS["neutral"],
                edgecolor="black",
            )
        ax2.set_xlabel("Allele Frequency")
        ax2.set_ylabel("Count")
        ax2.set_title("AF Spectrum")
    else:
        ax2.text(
            0.5,
            0.5,
            "Allele Frequency\nNot Available",
            ha="center",
            va="center",
            transform=ax2.transAxes,
        )

    # Plot 3: HWE p-values
    ax3 = fig.add_subplot(gs[0, 2])
    if "p_value_hwe" in df.columns:
        pvals = df["p_value_hwe"].dropna()
        pvals = pvals[pvals > 0]
        if len(pvals) > 0:
            log_pvals = -np.log10(pvals)
            ax3.hist(
                log_pvals,
                bins=30,
                alpha=0.7,
                color=QC_COLORS["neutral"],
                edgecolor="black",
            )
            _add_threshold_lines(ax3, "hwe_pvalue", thresholds)
        ax3.set_xlabel("-log10(HWE p-val)")
        ax3.set_ylabel("Count")
        ax3.set_title("HWE p-values")
    else:
        ax3.text(
            0.5,
            0.5,
            "HWE p-values\nNot Available",
            ha="center",
            va="center",
            transform=ax3.transAxes,
        )

    # Plot 4: Allele count distribution
    ax4 = fig.add_subplot(gs[1, 0])
    ac_col = "AC_alt" if "AC_alt" in df.columns else "AC"
    if ac_col in df.columns:
        if ac_col == "AC":
            try:
                acs = df[ac_col].apply(
                    lambda x: (
                        x[1] if isinstance(x, (list, np.ndarray)) and len(x) > 1 else x
                    )
                )
            except:
                acs = df[ac_col]
        else:
            acs = df[ac_col]

        acs_filtered = acs.dropna()
        acs_filtered = acs_filtered[acs_filtered > 0]
        if len(acs_filtered) > 0:
            ax4.hist(
                acs_filtered,
                bins=30,
                alpha=0.7,
                color=QC_COLORS["neutral"],
                edgecolor="black",
            )
        ax4.set_xlabel("Allele Count")
        ax4.set_ylabel("Count")
        ax4.set_title("Allele Count")
        ax4.set_yscale("log")
    else:
        ax4.text(
            0.5,
            0.5,
            "Allele Count\nNot Available",
            ha="center",
            va="center",
            transform=ax4.transAxes,
        )

    # Plot 5: Genotype distribution
    ax5 = fig.add_subplot(gs[1, 1])
    genotype_cols = ["n_het", "n_non_ref"]
    available_cols = [col for col in genotype_cols if col in df.columns]

    if available_cols:
        genotype_data = df[available_cols].sum()
        ax5.pie(genotype_data.values, labels=genotype_data.index, autopct="%1.1f%%")
        ax5.set_title("Genotype Distribution")
    else:
        ax5.text(
            0.5,
            0.5,
            "Genotype Counts\nNot Available",
            ha="center",
            va="center",
            transform=ax5.transAxes,
        )

    # Plot 6: Missing data pattern
    ax6 = fig.add_subplot(gs[1, 2])
    if "n_not_called" in df.columns and "n_called" in df.columns:
        missing_rate = df["n_not_called"] / (df["n_called"] + df["n_not_called"])
        missing_rate = missing_rate.dropna()
        if len(missing_rate) > 0:
            ax6.hist(
                missing_rate,
                bins=30,
                alpha=0.7,
                color=QC_COLORS["neutral"],
                edgecolor="black",
            )
        ax6.set_xlabel("Missing Rate")
        ax6.set_ylabel("Count")
        ax6.set_title("Missing Data")
    else:
        ax6.text(
            0.5,
            0.5,
            "Missing Data\nNot Available",
            ha="center",
            va="center",
            transform=ax6.transAxes,
        )

    # Add overall title
    fig.suptitle(
        f"Variant QC Overview (n={len(df)} variants)", fontsize=16, fontweight="bold"
    )

    if save_path:
        save_figure(fig, save_path, **kwargs)

    return fig


def plot_qc_summary_dashboard(
    qc_results,
    thresholds: Optional[Dict] = None,
    style: str = "default",
    figsize: Tuple[float, float] = (20, 12),
    save_path: Optional[Union[str, Path]] = None,
    **kwargs,
) -> plt.Figure:
    """
    Create a comprehensive QC summary dashboard combining sample and variant metrics.

    Parameters
    ----------
    qc_results : QCMetrics or dict
        QC results object or dictionary with 'sample' and 'variant' DataFrames
    thresholds : dict, optional
        Custom thresholds for different metrics
    style : str
        Plot style
    figsize : tuple
        Figure size (width, height)
    save_path : str or Path, optional
        Path to save the figure
    **kwargs
        Additional plotting arguments

    Returns
    -------
    matplotlib.figure.Figure
        The created figure
    """
    set_default_style(style)

    # Extract DataFrames from qc_results
    if hasattr(qc_results, "get_sample_metrics_df"):
        sample_df = qc_results.get_sample_metrics_df()
        variant_df = qc_results.get_variant_metrics_df()
    elif isinstance(qc_results, dict):
        sample_df = qc_results.get("sample")
        variant_df = qc_results.get("variant")
    else:
        raise ValueError(
            "qc_results must be QCMetrics object or dictionary with 'sample' and 'variant' DataFrames"
        )

    # Create figure with subplots
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(3, 4, figure=fig, hspace=0.4, wspace=0.3)

    # Sample QC plots (top row)
    if sample_df is not None and not sample_df.empty:
        sample_data = _prepare_sample_qc_data(sample_df)

        # Sample call rate
        ax1 = fig.add_subplot(gs[0, 0])
        if "call_rate" in sample_data.columns:
            call_rates = sample_data["call_rate"].dropna()
            ax1.hist(
                call_rates,
                bins=20,
                alpha=0.7,
                color=QC_COLORS["neutral"],
                edgecolor="black",
            )
            _add_threshold_lines(ax1, "sample_call_rate", thresholds)
            ax1.set_xlabel("Sample Call Rate")
            ax1.set_ylabel("Count")
            ax1.set_title("Sample Call Rate")

        # Sample Ti/Tv
        ax2 = fig.add_subplot(gs[0, 1])
        if "r_ti_tv" in sample_data.columns:
            ti_tv = sample_data["r_ti_tv"].dropna()
            q1, q99 = np.percentile(ti_tv, [1, 99])
            filtered_ti_tv = ti_tv[(ti_tv >= q1) & (ti_tv <= q99)]
            ax2.hist(
                filtered_ti_tv,
                bins=20,
                alpha=0.7,
                color=QC_COLORS["neutral"],
                edgecolor="black",
            )
            ax2.set_xlabel("Ti/Tv Ratio")
            ax2.set_ylabel("Count")
            ax2.set_title("Sample Ti/Tv")

        # Sample depth
        ax3 = fig.add_subplot(gs[0, 2])
        depth_col = None
        for col in ["mean", "dp_stats_mean", "sample_qc_mean_dp"]:
            if col in sample_data.columns:
                depth_col = col
                break

        if depth_col:
            depths = sample_data[depth_col].dropna()
            ax3.hist(
                depths,
                bins=20,
                alpha=0.7,
                color=QC_COLORS["neutral"],
                edgecolor="black",
            )
            ax3.set_xlabel("Mean Depth")
            ax3.set_ylabel("Count")
            ax3.set_title("Sample Depth")

        # Sample genotype summary
        ax4 = fig.add_subplot(gs[0, 3])
        genotype_cols = ["n_hom_ref", "n_het", "n_hom_var"]
        available_cols = [col for col in genotype_cols if col in sample_data.columns]

        if available_cols:
            genotype_data = sample_data[available_cols].sum()
            colors = [QC_COLORS["pass"], QC_COLORS["warn"], QC_COLORS["fail"]]
            ax4.pie(
                genotype_data.values,
                labels=genotype_data.index,
                colors=colors[: len(genotype_data)],
                autopct="%1.0f%%",
            )
            ax4.set_title("Sample Genotypes")

    # Variant QC plots (middle row)
    if variant_df is not None and not variant_df.empty:
        variant_data = _prepare_variant_qc_data(variant_df)

        # Variant call rate
        ax5 = fig.add_subplot(gs[1, 0])
        if "call_rate" in variant_data.columns:
            call_rates = variant_data["call_rate"].dropna()
            ax5.hist(
                call_rates,
                bins=20,
                alpha=0.7,
                color=QC_COLORS["neutral"],
                edgecolor="black",
            )
            _add_threshold_lines(ax5, "variant_call_rate", thresholds)
            ax5.set_xlabel("Variant Call Rate")
            ax5.set_ylabel("Count")
            ax5.set_title("Variant Call Rate")

        # Allele frequency spectrum
        ax6 = fig.add_subplot(gs[1, 1])
        af_col = "AF_alt" if "AF_alt" in variant_data.columns else "AF"
        if af_col in variant_data.columns:
            if af_col == "AF":
                try:
                    afs = variant_data[af_col].apply(
                        lambda x: (
                            x[1]
                            if isinstance(x, (list, np.ndarray)) and len(x) > 1
                            else x
                        )
                    )
                except:
                    afs = variant_data[af_col]
            else:
                afs = variant_data[af_col]

            afs_filtered = afs.dropna()
            afs_filtered = afs_filtered[(afs_filtered > 0) & (afs_filtered < 1)]
            if len(afs_filtered) > 0:
                ax6.hist(
                    afs_filtered,
                    bins=20,
                    alpha=0.7,
                    color=QC_COLORS["neutral"],
                    edgecolor="black",
                )
            ax6.set_xlabel("Allele Frequency")
            ax6.set_ylabel("Count")
            ax6.set_title("AF Spectrum")

        # HWE p-values
        ax7 = fig.add_subplot(gs[1, 2])
        if "p_value_hwe" in variant_data.columns:
            pvals = variant_data["p_value_hwe"].dropna()
            pvals = pvals[pvals > 0]
            if len(pvals) > 0:
                log_pvals = -np.log10(pvals)
                ax7.hist(
                    log_pvals,
                    bins=20,
                    alpha=0.7,
                    color=QC_COLORS["neutral"],
                    edgecolor="black",
                )
            ax7.set_xlabel("-log10(HWE p-val)")
            ax7.set_ylabel("Count")
            ax7.set_title("HWE p-values")

        # Allele count distribution
        ax8 = fig.add_subplot(gs[1, 3])
        ac_col = "AC_alt" if "AC_alt" in variant_data.columns else "AC"
        if ac_col in variant_data.columns:
            if ac_col == "AC":
                try:
                    acs = variant_data[ac_col].apply(
                        lambda x: (
                            x[1]
                            if isinstance(x, (list, np.ndarray)) and len(x) > 1
                            else x
                        )
                    )
                except:
                    acs = variant_data[ac_col]
            else:
                acs = variant_data[ac_col]

            acs_filtered = acs.dropna()
            acs_filtered = acs_filtered[acs_filtered > 0]
            if len(acs_filtered) > 0:
                ax8.hist(
                    acs_filtered,
                    bins=20,
                    alpha=0.7,
                    color=QC_COLORS["neutral"],
                    edgecolor="black",
                )
                ax8.set_yscale("log")
            ax8.set_xlabel("Allele Count")
            ax8.set_ylabel("Count (log)")
            ax8.set_title("Allele Count")

    # Summary statistics (bottom row)
    ax9 = fig.add_subplot(gs[2, :2])
    ax9.axis("off")

    # Create summary text
    summary_text = "QC Summary Statistics\n" + "=" * 50 + "\n\n"

    if sample_df is not None and not sample_df.empty:
        sample_data = _prepare_sample_qc_data(sample_df)
        summary_text += f"SAMPLES (n={len(sample_data)}):\n"

        if "call_rate" in sample_data.columns:
            cr = sample_data["call_rate"].dropna()
            summary_text += f"  Call Rate: {cr.mean():.3f} ± {cr.std():.3f}\n"

        if "r_ti_tv" in sample_data.columns:
            titv = sample_data["r_ti_tv"].dropna()
            summary_text += f"  Ti/Tv Ratio: {titv.mean():.3f} ± {titv.std():.3f}\n"

        depth_col = None
        for col in ["dp_stats_mean", "mean", "dp_mean"]:
            if col in sample_data.columns:
                depth_col = col
                break

        if depth_col:
            depth = sample_data[depth_col].dropna()
            summary_text += f"  Mean Depth: {depth.mean():.1f} ± {depth.std():.1f}\n"

    summary_text += "\n"

    if variant_df is not None and not variant_df.empty:
        variant_data = _prepare_variant_qc_data(variant_df)
        summary_text += f"VARIANTS (n={len(variant_data)}):\n"

        if "call_rate" in variant_data.columns:
            cr = variant_data["call_rate"].dropna()
            summary_text += f"  Call Rate: {cr.mean():.3f} ± {cr.std():.3f}\n"

        af_col = "AF_alt" if "AF_alt" in variant_data.columns else "AF"
        if af_col in variant_data.columns:
            if af_col == "AF":
                try:
                    afs = variant_data[af_col].apply(
                        lambda x: (
                            x[1]
                            if isinstance(x, (list, np.ndarray)) and len(x) > 1
                            else x
                        )
                    )
                except:
                    afs = variant_data[af_col]
            else:
                afs = variant_data[af_col]

            afs = afs.dropna()
            afs = afs[(afs > 0) & (afs < 1)]
            if len(afs) > 0:
                summary_text += f"  Mean AF: {afs.mean():.4f}\n"
                summary_text += f"  Rare variants (<1%): {(afs < 0.01).sum()}\n"

        if "p_value_hwe" in variant_data.columns:
            pvals = variant_data["p_value_hwe"].dropna()
            pvals = pvals[pvals > 0]
            if len(pvals) > 0:
                failing_hwe = (pvals < 1e-6).sum()
                summary_text += f"  HWE failures (p<1e-6): {failing_hwe}\n"

    ax9.text(
        0.05,
        0.95,
        summary_text,
        transform=ax9.transAxes,
        verticalalignment="top",
        fontfamily="monospace",
        fontsize=10,
        bbox=dict(boxstyle="round", facecolor="lightgray", alpha=0.8),
    )

    # Add quality indicators (bottom right)
    ax10 = fig.add_subplot(gs[2, 2:])
    ax10.axis("off")

    # Create quality assessment
    quality_text = "Quality Assessment\n" + "=" * 30 + "\n\n"

    # Sample quality indicators
    if sample_df is not None and not sample_df.empty:
        sample_data = _prepare_sample_qc_data(sample_df)
        quality_text += "SAMPLE QUALITY:\n"

        if "call_rate" in sample_data.columns:
            cr = sample_data["call_rate"].dropna()
            good_samples = (cr >= 0.95).sum()
            total_samples = len(cr)
            pct = good_samples / total_samples * 100
            status = "✓" if pct >= 90 else "⚠" if pct >= 75 else "✗"
            quality_text += f"  Call Rate ≥95%: {status} {pct:.1f}%\n"

        if "r_ti_tv" in sample_data.columns:
            titv = sample_data["r_ti_tv"].dropna()
            good_titv = ((titv >= 1.8) & (titv <= 2.2)).sum()
            total = len(titv)
            pct = good_titv / total * 100 if total > 0 else 0
            status = "✓" if pct >= 80 else "⚠" if pct >= 60 else "✗"
            quality_text += f"  Ti/Tv in range: {status} {pct:.1f}%\n"

    quality_text += "\n"

    # Variant quality indicators
    if variant_df is not None and not variant_df.empty:
        variant_data = _prepare_variant_qc_data(variant_df)
        quality_text += "VARIANT QUALITY:\n"

        if "call_rate" in variant_data.columns:
            cr = variant_data["call_rate"].dropna()
            good_variants = (cr >= 0.95).sum()
            total_variants = len(cr)
            pct = good_variants / total_variants * 100
            status = "✓" if pct >= 80 else "⚠" if pct >= 60 else "✗"
            quality_text += f"  Call Rate ≥95%: {status} {pct:.1f}%\n"

        if "p_value_hwe" in variant_data.columns:
            pvals = variant_data["p_value_hwe"].dropna()
            pvals = pvals[pvals > 0]
            if len(pvals) > 0:
                good_hwe = (pvals >= 1e-6).sum()
                total = len(pvals)
                pct = good_hwe / total * 100
                status = "✓" if pct >= 95 else "⚠" if pct >= 90 else "✗"
                quality_text += f"  HWE p≥1e-6: {status} {pct:.1f}%\n"

    ax10.text(
        0.05,
        0.95,
        quality_text,
        transform=ax10.transAxes,
        verticalalignment="top",
        fontfamily="monospace",
        fontsize=10,
        bbox=dict(boxstyle="round", facecolor="lightblue", alpha=0.8),
    )

    # Add overall title
    fig.suptitle("Quality Control Summary Dashboard", fontsize=18, fontweight="bold")

    if save_path:
        save_figure(fig, save_path, **kwargs)

    return fig
