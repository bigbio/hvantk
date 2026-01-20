"""
Interactive QC Plots using Plotly

This module provides interactive QC visualizations using Plotly for genomic variant data.
Interactive plots allow zooming, panning, hover tooltips, and data exploration.

The module provides:
- Interactive histograms with hover tooltips
- Multi-panel dashboards with synchronized controls
- Scatter plots for correlation analysis
- Zoomable and pannable plots
- HTML export for sharing

Requirements:
    plotly: pip install plotly

Example:
    >>> from hvantk.hgc import compute_full_qc
    >>> from hvantk.visualization.interactive_qc import plot_interactive_qc_dashboard
    >>> qc_results = compute_full_qc(mt)
    >>> fig = plot_interactive_qc_dashboard(qc_results)
    >>> fig.show()
"""

import logging
from typing import Optional, Dict, Union, Tuple
from pathlib import Path

import numpy as np
import pandas as pd

# Try to import plotly components
try:
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots

    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False

from .qc_plots import (
    _prepare_sample_qc_data,
    _prepare_variant_qc_data,
    QC_COLORS,
    QC_THRESHOLDS,
)

logger = logging.getLogger(__name__)

# Plotly color scheme (matching our QC colors)
PLOTLY_COLORS = {
    "pass": "#2E8B57",  # Sea Green
    "warn": "#FF8C00",  # Dark Orange
    "fail": "#DC143C",  # Crimson
    "neutral": "#4682B4",  # Steel Blue
    "highlight": "#FFD700",  # Gold
    "background": "#F5F5F5",  # White Smoke
}


def check_plotly_available():
    """Check if plotly is available and raise informative error if not."""
    if not HAS_PLOTLY:
        raise ImportError(
            "Plotly is required for interactive plotting. "
            "Install with: pip install plotly"
        )


def plot_interactive_sample_call_rates(
    sample_df: pd.DataFrame,
    thresholds: Optional[Dict] = None,
    title: str = "Sample Call Rate Distribution",
    **kwargs,
) -> go.Figure:
    """
    Create interactive sample call rate distribution plot.

    Parameters
    ----------
    sample_df : pd.DataFrame
        Sample QC metrics DataFrame
    thresholds : dict, optional
        Custom thresholds for highlighting
    title : str
        Plot title
    **kwargs
        Additional plotly arguments

    Returns
    -------
    plotly.graph_objects.Figure
        Interactive plotly figure
    """
    check_plotly_available()

    # Prepare data
    df = _prepare_sample_qc_data(sample_df)

    if "call_rate" not in df.columns:
        raise ValueError("Sample DataFrame must contain 'call_rate' column")

    call_rates = df["call_rate"].dropna()

    # Create histogram
    fig = go.Figure()

    # Add main histogram
    fig.add_trace(
        go.Histogram(
            x=call_rates,
            nbinsx=50,
            name="Call Rates",
            marker_color=PLOTLY_COLORS["neutral"],
            opacity=0.7,
            hovertemplate="<b>Call Rate Range:</b> %{x}<br>"
            + "<b>Sample Count:</b> %{y}<br>"
            + "<extra></extra>",
        )
    )

    # Add threshold lines
    if thresholds is None:
        thresholds = QC_THRESHOLDS.get("sample_call_rate", {})

    if "good" in thresholds:
        fig.add_vline(
            x=thresholds["good"],
            line_dash="dash",
            line_color=PLOTLY_COLORS["pass"],
            annotation_text=f"Good (≥{thresholds['good']:.2f})",
            annotation_position="top",
        )

    if "acceptable" in thresholds:
        fig.add_vline(
            x=thresholds["acceptable"],
            line_dash="dash",
            line_color=PLOTLY_COLORS["warn"],
            annotation_text=f"Acceptable (≥{thresholds['acceptable']:.2f})",
            annotation_position="top",
        )

    # Customize layout
    fig.update_layout(
        title={
            "text": f"{title}<br><sup>n={len(call_rates)} samples | Mean: {call_rates.mean():.3f} ± {call_rates.std():.3f}</sup>",
            "x": 0.5,
            "xanchor": "center",
        },
        xaxis_title="Sample Call Rate",
        yaxis_title="Number of Samples",
        template="plotly_white",
        hovermode="x unified",
        **kwargs,
    )

    return fig


def plot_interactive_sample_titv(
    sample_df: pd.DataFrame,
    thresholds: Optional[Dict] = None,
    title: str = "Sample Ti/Tv Ratio Distribution",
    **kwargs,
) -> go.Figure:
    """
    Create interactive sample Ti/Tv ratio distribution plot.

    Parameters
    ----------
    sample_df : pd.DataFrame
        Sample QC metrics DataFrame
    thresholds : dict, optional
        Custom thresholds for highlighting
    title : str
        Plot title
    **kwargs
        Additional plotly arguments

    Returns
    -------
    plotly.graph_objects.Figure
        Interactive plotly figure
    """
    check_plotly_available()

    # Prepare data
    df = _prepare_sample_qc_data(sample_df)

    if "r_ti_tv" not in df.columns:
        raise ValueError("Sample DataFrame must contain 'r_ti_tv' column")

    ti_tv_ratios = df["r_ti_tv"].dropna()

    # Filter extreme outliers for better visualization
    q1, q99 = np.percentile(ti_tv_ratios, [1, 99])
    filtered_ratios = ti_tv_ratios[(ti_tv_ratios >= q1) & (ti_tv_ratios <= q99)]

    # Create histogram
    fig = go.Figure()

    # Add main histogram
    fig.add_trace(
        go.Histogram(
            x=filtered_ratios,
            nbinsx=50,
            name="Ti/Tv Ratios",
            marker_color=PLOTLY_COLORS["neutral"],
            opacity=0.7,
            hovertemplate="<b>Ti/Tv Range:</b> %{x}<br>"
            + "<b>Sample Count:</b> %{y}<br>"
            + "<extra></extra>",
        )
    )

    # Add expected range for Ti/Tv
    if thresholds is None:
        thresholds = QC_THRESHOLDS.get("ti_tv_ratio", {})

    if "good_min" in thresholds and "good_max" in thresholds:
        fig.add_vrect(
            x0=thresholds["good_min"],
            x1=thresholds["good_max"],
            fillcolor=PLOTLY_COLORS["pass"],
            opacity=0.2,
            layer="below",
            line_width=0,
            annotation_text="Expected Range",
            annotation_position="top left",
        )

    # Customize layout
    fig.update_layout(
        title={
            "text": f"{title}<br><sup>n={len(ti_tv_ratios)} samples | Mean: {ti_tv_ratios.mean():.3f} ± {ti_tv_ratios.std():.3f}</sup>",
            "x": 0.5,
            "xanchor": "center",
        },
        xaxis_title="Ti/Tv Ratio",
        yaxis_title="Number of Samples",
        template="plotly_white",
        hovermode="x unified",
        **kwargs,
    )

    return fig


def plot_interactive_variant_call_rates(
    variant_df: pd.DataFrame,
    thresholds: Optional[Dict] = None,
    title: str = "Variant Call Rate Distribution",
    **kwargs,
) -> go.Figure:
    """
    Create interactive variant call rate distribution plot.

    Parameters
    ----------
    variant_df : pd.DataFrame
        Variant QC metrics DataFrame
    thresholds : dict, optional
        Custom thresholds for highlighting
    title : str
        Plot title
    **kwargs
        Additional plotly arguments

    Returns
    -------
    plotly.graph_objects.Figure
        Interactive plotly figure
    """
    check_plotly_available()

    # Prepare data
    df = _prepare_variant_qc_data(variant_df)

    if "call_rate" not in df.columns:
        raise ValueError("Variant DataFrame must contain 'call_rate' column")

    call_rates = df["call_rate"].dropna()

    # Create histogram
    fig = go.Figure()

    # Add main histogram
    fig.add_trace(
        go.Histogram(
            x=call_rates,
            nbinsx=50,
            name="Call Rates",
            marker_color=PLOTLY_COLORS["neutral"],
            opacity=0.7,
            hovertemplate="<b>Call Rate Range:</b> %{x}<br>"
            + "<b>Variant Count:</b> %{y}<br>"
            + "<extra></extra>",
        )
    )

    # Add threshold lines
    if thresholds is None:
        thresholds = QC_THRESHOLDS.get("variant_call_rate", {})

    if "good" in thresholds:
        fig.add_vline(
            x=thresholds["good"],
            line_dash="dash",
            line_color=PLOTLY_COLORS["pass"],
            annotation_text=f"Good (≥{thresholds['good']:.2f})",
            annotation_position="top",
        )

    if "acceptable" in thresholds:
        fig.add_vline(
            x=thresholds["acceptable"],
            line_dash="dash",
            line_color=PLOTLY_COLORS["warn"],
            annotation_text=f"Acceptable (≥{thresholds['acceptable']:.2f})",
            annotation_position="top",
        )

    # Customize layout
    fig.update_layout(
        title={
            "text": f"{title}<br><sup>n={len(call_rates)} variants | Mean: {call_rates.mean():.3f} ± {call_rates.std():.3f}</sup>",
            "x": 0.5,
            "xanchor": "center",
        },
        xaxis_title="Variant Call Rate",
        yaxis_title="Number of Variants",
        template="plotly_white",
        hovermode="x unified",
        **kwargs,
    )

    return fig


def plot_interactive_allele_frequencies(
    variant_df: pd.DataFrame,
    title: str = "Allele Frequency Spectrum",
    log_scale: bool = True,
    **kwargs,
) -> go.Figure:
    """
    Create interactive allele frequency spectrum plot.

    Parameters
    ----------
    variant_df : pd.DataFrame
        Variant QC metrics DataFrame
    title : str
        Plot title
    log_scale : bool
        Whether to use log scale for y-axis
    **kwargs
        Additional plotly arguments

    Returns
    -------
    plotly.graph_objects.Figure
        Interactive plotly figure
    """
    check_plotly_available()

    # Prepare data
    df = _prepare_variant_qc_data(variant_df)

    # Look for allele frequency column
    af_col = "AF_alt" if "AF_alt" in df.columns else "AF"
    if af_col not in df.columns:
        raise ValueError("Variant DataFrame must contain allele frequency column")

    # Get allele frequencies
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

    afs = afs.dropna()
    afs_filtered = afs[(afs > 0) & (afs < 1)]

    # Create histogram
    fig = go.Figure()

    # Add main histogram
    fig.add_trace(
        go.Histogram(
            x=afs_filtered,
            nbinsx=50,
            name="Allele Frequencies",
            marker_color=PLOTLY_COLORS["neutral"],
            opacity=0.7,
            hovertemplate="<b>AF Range:</b> %{x}<br>"
            + "<b>Variant Count:</b> %{y}<br>"
            + "<extra></extra>",
        )
    )

    # Add frequency category lines
    fig.add_vline(
        x=0.01,
        line_dash="dash",
        line_color=PLOTLY_COLORS["warn"],
        annotation_text="1% MAF",
        annotation_position="top",
    )

    fig.add_vline(
        x=0.05,
        line_dash="dash",
        line_color=PLOTLY_COLORS["pass"],
        annotation_text="5% MAF",
        annotation_position="top",
    )

    # Customize layout
    yaxis_type = "log" if log_scale else "linear"

    fig.update_layout(
        title={
            "text": f"{title}<br><sup>n={len(afs_filtered)} polymorphic variants | Mean AF: {afs_filtered.mean():.4f}</sup>",
            "x": 0.5,
            "xanchor": "center",
        },
        xaxis_title="Allele Frequency",
        yaxis_title="Number of Variants",
        yaxis_type=yaxis_type,
        template="plotly_white",
        hovermode="x unified",
        **kwargs,
    )

    return fig


def plot_interactive_hwe_pvalues(
    variant_df: pd.DataFrame,
    thresholds: Optional[Dict] = None,
    title: str = "Hardy-Weinberg Equilibrium p-values",
    log_transform: bool = True,
    **kwargs,
) -> go.Figure:
    """
    Create interactive Hardy-Weinberg equilibrium p-values plot.

    Parameters
    ----------
    variant_df : pd.DataFrame
        Variant QC metrics DataFrame
    thresholds : dict, optional
        Custom HWE p-value thresholds
    title : str
        Plot title
    log_transform : bool
        Whether to plot -log10(p-values)
    **kwargs
        Additional plotly arguments

    Returns
    -------
    plotly.graph_objects.Figure
        Interactive plotly figure
    """
    check_plotly_available()

    # Prepare data
    df = _prepare_variant_qc_data(variant_df)

    if "p_value_hwe" not in df.columns:
        raise ValueError("Variant DataFrame must contain 'p_value_hwe' column")

    pvalues = df["p_value_hwe"].dropna()
    pvalues = pvalues[pvalues > 0]  # Remove p=0 values

    if log_transform:
        plot_values = -np.log10(pvalues)
        xlabel = "-log10(HWE p-value)"
        title_suffix = " (-log10 scale)"
    else:
        plot_values = pvalues
        xlabel = "HWE p-value"
        title_suffix = ""

    # Create histogram
    fig = go.Figure()

    # Add main histogram
    fig.add_trace(
        go.Histogram(
            x=plot_values,
            nbinsx=50,
            name="HWE p-values",
            marker_color=PLOTLY_COLORS["neutral"],
            opacity=0.7,
            hovertemplate=f"<b>{xlabel} Range:</b> %{{x}}<br>"
            + "<b>Variant Count:</b> %{y}<br>"
            + "<extra></extra>",
        )
    )

    # Add threshold lines
    if thresholds is None:
        thresholds = QC_THRESHOLDS.get("hwe_pvalue", {})

    for name, threshold in thresholds.items():
        if log_transform:
            line_val = -np.log10(threshold)
        else:
            line_val = threshold

        color = PLOTLY_COLORS["fail"] if name == "fail" else PLOTLY_COLORS["warn"]
        fig.add_vline(
            x=line_val,
            line_dash="dash",
            line_color=color,
            annotation_text=f"{name.title()}: {threshold:.0e}",
            annotation_position="top",
        )

    # Customize layout
    fig.update_layout(
        title={
            "text": f"{title}{title_suffix}<br><sup>n={len(pvalues)} variants</sup>",
            "x": 0.5,
            "xanchor": "center",
        },
        xaxis_title=xlabel,
        yaxis_title="Number of Variants",
        template="plotly_white",
        hovermode="x unified",
        **kwargs,
    )

    return fig


def plot_interactive_sample_scatter(
    sample_df: pd.DataFrame,
    x_metric: str = "call_rate",
    y_metric: str = "r_ti_tv",
    color_metric: Optional[str] = None,
    title: Optional[str] = None,
    **kwargs,
) -> go.Figure:
    """
    Create interactive scatter plot of sample QC metrics.

    Parameters
    ----------
    sample_df : pd.DataFrame
        Sample QC metrics DataFrame
    x_metric : str
        Column name for x-axis metric
    y_metric : str
        Column name for y-axis metric
    color_metric : str, optional
        Column name for color mapping
    title : str, optional
        Plot title
    **kwargs
        Additional plotly arguments

    Returns
    -------
    plotly.graph_objects.Figure
        Interactive plotly figure
    """
    check_plotly_available()

    # Prepare data
    df = _prepare_sample_qc_data(sample_df)

    # Check required columns
    required_cols = [x_metric, y_metric]
    if color_metric:
        required_cols.append(color_metric)

    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing columns: {missing_cols}")

    # Create scatter plot
    if color_metric:
        fig = px.scatter(
            df,
            x=x_metric,
            y=y_metric,
            color=color_metric,
            color_continuous_scale="viridis",
            hover_data=[
                col
                for col in df.columns
                if col not in [x_metric, y_metric, color_metric]
            ][:5],
            **kwargs,
        )
    else:
        fig = px.scatter(
            df,
            x=x_metric,
            y=y_metric,
            hover_data=[col for col in df.columns if col not in [x_metric, y_metric]][
                :5
            ],
            **kwargs,
        )
        fig.update_traces(marker_color=PLOTLY_COLORS["neutral"])

    # Customize layout
    if title is None:
        title = f"Sample {y_metric.replace('_', ' ').title()} vs {x_metric.replace('_', ' ').title()}"

    fig.update_layout(
        title={
            "text": f"{title}<br><sup>n={len(df)} samples</sup>",
            "x": 0.5,
            "xanchor": "center",
        },
        xaxis_title=x_metric.replace("_", " ").title(),
        yaxis_title=y_metric.replace("_", " ").title(),
        template="plotly_white",
    )

    return fig


def plot_interactive_qc_dashboard(
    qc_results, title: str = "Interactive QC Dashboard", **kwargs
) -> go.Figure:
    """
    Create comprehensive interactive QC dashboard.

    Parameters
    ----------
    qc_results : QCMetrics
        QC results object with computed metrics
    title : str
        Dashboard title
    **kwargs
        Additional plotly arguments

    Returns
    -------
    plotly.graph_objects.Figure
        Interactive plotly dashboard
    """
    check_plotly_available()

    # Get data
    sample_df = qc_results.get_sample_metrics_df() if qc_results.has_sample_qc else None
    variant_df = (
        qc_results.get_variant_metrics_df() if qc_results.has_variant_qc else None
    )

    # Create subplots
    subplot_titles = []
    if sample_df is not None:
        subplot_titles.extend(["Sample Call Rates", "Sample Ti/Tv Ratios"])
    if variant_df is not None:
        subplot_titles.extend(["Variant Call Rates", "Allele Frequencies"])

    n_plots = len(subplot_titles)
    cols = 2
    rows = (n_plots + cols - 1) // cols

    fig = make_subplots(
        rows=rows,
        cols=cols,
        subplot_titles=subplot_titles,
        horizontal_spacing=0.1,
        vertical_spacing=0.12,
    )

    plot_idx = 1

    # Sample plots
    if sample_df is not None:
        sample_data = _prepare_sample_qc_data(sample_df)

        # Sample call rates
        if "call_rate" in sample_data.columns:
            call_rates = sample_data["call_rate"].dropna()

            fig.add_trace(
                go.Histogram(
                    x=call_rates,
                    nbinsx=30,
                    name="Sample Call Rates",
                    marker_color=PLOTLY_COLORS["neutral"],
                    opacity=0.7,
                    showlegend=False,
                    hovertemplate="Call Rate: %{x}<br>Count: %{y}<extra></extra>",
                ),
                row=(plot_idx - 1) // cols + 1,
                col=(plot_idx - 1) % cols + 1,
            )
            plot_idx += 1

        # Sample Ti/Tv ratios
        if "r_ti_tv" in sample_data.columns:
            ti_tv = sample_data["r_ti_tv"].dropna()
            q1, q99 = np.percentile(ti_tv, [1, 99])
            filtered_ti_tv = ti_tv[(ti_tv >= q1) & (ti_tv <= q99)]

            fig.add_trace(
                go.Histogram(
                    x=filtered_ti_tv,
                    nbinsx=30,
                    name="Ti/Tv Ratios",
                    marker_color=PLOTLY_COLORS["pass"],
                    opacity=0.7,
                    showlegend=False,
                    hovertemplate="Ti/Tv Ratio: %{x}<br>Count: %{y}<extra></extra>",
                ),
                row=(plot_idx - 1) // cols + 1,
                col=(plot_idx - 1) % cols + 1,
            )
            plot_idx += 1

    # Variant plots
    if variant_df is not None:
        variant_data = _prepare_variant_qc_data(variant_df)

        # Variant call rates
        if "call_rate" in variant_data.columns:
            call_rates = variant_data["call_rate"].dropna()

            fig.add_trace(
                go.Histogram(
                    x=call_rates,
                    nbinsx=30,
                    name="Variant Call Rates",
                    marker_color=PLOTLY_COLORS["warn"],
                    opacity=0.7,
                    showlegend=False,
                    hovertemplate="Call Rate: %{x}<br>Count: %{y}<extra></extra>",
                ),
                row=(plot_idx - 1) // cols + 1,
                col=(plot_idx - 1) % cols + 1,
            )
            plot_idx += 1

        # Allele frequencies
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
                fig.add_trace(
                    go.Histogram(
                        x=afs_filtered,
                        nbinsx=30,
                        name="Allele Frequencies",
                        marker_color=PLOTLY_COLORS["highlight"],
                        opacity=0.7,
                        showlegend=False,
                        hovertemplate="Allele Freq: %{x}<br>Count: %{y}<extra></extra>",
                    ),
                    row=(plot_idx - 1) // cols + 1,
                    col=(plot_idx - 1) % cols + 1,
                )

                # Set log scale for AF plot
                fig.update_yaxes(
                    type="log",
                    row=(plot_idx - 1) // cols + 1,
                    col=(plot_idx - 1) % cols + 1,
                )

    # Update layout
    fig.update_layout(
        title={"text": title, "x": 0.5, "xanchor": "center", "font": {"size": 20}},
        template="plotly_white",
        height=300 * rows,
        **kwargs,
    )

    return fig


def save_interactive_plot(
    fig: go.Figure, output_path: Union[str, Path], format: str = "html", **kwargs
):
    """
    Save interactive plotly figure to file.

    Parameters
    ----------
    fig : plotly.graph_objects.Figure
        Plotly figure to save
    output_path : str or Path
        Output file path
    format : str
        Output format ('html', 'png', 'pdf', 'svg')
    **kwargs
        Additional arguments for plotly save functions
    """
    output_path = Path(output_path)

    if format == "html":
        fig.write_html(str(output_path), **kwargs)
    elif format in ["png", "pdf", "svg"]:
        fig.write_image(str(output_path), format=format, **kwargs)
    else:
        raise ValueError(f"Unsupported format: {format}")

    logger.info(f"Interactive plot saved: {output_path}")
