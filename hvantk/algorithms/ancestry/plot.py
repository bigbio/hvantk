"""Visualization functions for ancestry inference.

This module provides plotting utilities for ancestry inference results,
including PCA scatter plots, variance explained charts, ancestry distributions,
and confusion matrices.

The main plotting function `plot_pca_scatter` supports:
- Dual encoding: color by population, shape by source (query/reference)
- Filtering to specific populations for zoom views
- Query samples shown as distinct markers ("Undefined" or predicted ancestry)
- Customizable color palettes and marker styles
"""

import logging
from typing import Dict, List, Optional, Tuple, Union, Any

import numpy as np
import pandas as pd

from hvantk.algorithms.ancestry.constants import (
    SUPERPOP_COLORS,
    POPULATION_COLORS,
    POPULATION_NAMES,
    SOURCE_COL,
    KNOWN_ANCESTRY_COL,
    PREDICTED_ANCESTRY_COL,
    ANCESTRY_PROB_COL,
)

logger = logging.getLogger(__name__)

# Marker styles for different sources
DEFAULT_MARKERS = {
    "reference": "o",  # Filled circle
    "query": "x",  # X marker
}

# Default figure size
DEFAULT_FIGSIZE = (10, 8)


def _get_matplotlib():
    """Lazily import matplotlib to avoid import errors if not installed."""
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from matplotlib.lines import Line2D

        return plt, mpatches, Line2D
    except ImportError as e:
        raise ImportError(
            "matplotlib is required for plotting. "
            "Install with: pip install matplotlib"
        ) from e


def plot_pca_scatter(
    scores_df: pd.DataFrame,
    pc_x: int = 1,
    pc_y: int = 2,
    color_by: str = "ancestry",
    shape_by: str = "source",
    colors: Optional[Dict[str, str]] = None,
    markers: Optional[Dict[str, str]] = None,
    filter_populations: Optional[List[str]] = None,
    filter_source: Optional[str] = None,
    show_query_as_undefined: bool = False,
    query_label: str = "Query",
    reference_label: str = "Reference",
    undefined_label: str = "Undefined",
    undefined_color: str = "#8B4513",  # Brown (matches reference plot)
    title: Optional[str] = None,
    figsize: Tuple[float, float] = DEFAULT_FIGSIZE,
    alpha: float = 0.7,
    s_reference: int = 40,
    s_query: int = 60,
    legend_loc: str = "best",
    show_legend: bool = True,
    ax: Optional[Any] = None,
) -> Any:
    """Create PCA scatter plot with dual encoding (color + shape).

    This function creates publication-quality PCA scatter plots similar to
    standard population genetics figures. Reference samples are shown as
    filled circles colored by their known ancestry, while query samples
    are shown as X markers.

    Parameters
    ----------
    scores_df : pd.DataFrame
        DataFrame with PC scores. Must contain columns:
        - PC1, PC2, ... (PC scores)
        - _ancestry_source (str): "query" or "reference"
        - _known_ancestry (str): Known ancestry for reference samples
        - predicted_ancestry (str): Predicted ancestry (optional)
    pc_x : int
        Principal component for x-axis (1-indexed). Default: 1.
    pc_y : int
        Principal component for y-axis (1-indexed). Default: 2.
    color_by : str
        How to color points. Options:
        - "ancestry": Reference by known, query by predicted (or undefined)
        - "source": Color by query vs reference
        - "predicted": All by predicted ancestry
        Default: "ancestry".
    shape_by : str
        How to shape points. Options:
        - "source": Different shapes for query vs reference
        - None: Same shape for all
        Default: "source".
    colors : dict, optional
        Population -> color mapping. Defaults to SUPERPOP_COLORS.
    markers : dict, optional
        Source -> marker mapping. Defaults to {"reference": "o", "query": "x"}.
    filter_populations : list, optional
        Only show samples from these populations (for zoom views).
    filter_source : str, optional
        Only show "query" or "reference" samples.
    show_query_as_undefined : bool
        If True, show all query samples as "Undefined" regardless of prediction.
        Default: False.
    query_label : str
        Label for query samples in legend. Default: "Query".
    reference_label : str
        Label for reference samples in legend. Default: "Reference".
    undefined_label : str
        Label for undefined/unassigned samples. Default: "Undefined".
    undefined_color : str
        Color for undefined samples. Default: "#8B4513" (brown).
    title : str, optional
        Plot title. Auto-generated if None.
    figsize : tuple
        Figure size (width, height). Default: (10, 8).
    alpha : float
        Point transparency. Default: 0.7.
    s_reference : int
        Marker size for reference samples. Default: 40.
    s_query : int
        Marker size for query samples. Default: 60.
    legend_loc : str
        Legend location. Default: "best".
    show_legend : bool
        Whether to show legend. Default: True.
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. Creates new figure if None.

    Returns
    -------
    matplotlib.figure.Figure
        The matplotlib figure object.

    Example
    -------
    >>> # Basic plot
    >>> fig = plot_pca_scatter(result.get_predictions_df())
    >>> fig.savefig("pca_plot.png", dpi=300)

    >>> # Zoom to European cluster
    >>> fig = plot_pca_scatter(
    ...     result.get_predictions_df(),
    ...     filter_populations=["EUR"],
    ...     title="European Cluster"
    ... )
    """
    plt, mpatches, Line2D = _get_matplotlib()

    # Set up colors and markers
    if colors is None:
        colors = SUPERPOP_COLORS.copy()
    if markers is None:
        markers = DEFAULT_MARKERS.copy()

    # Add undefined color
    colors[undefined_label] = undefined_color
    colors["unassigned"] = undefined_color

    # Get PC column names
    pc_x_col = f"PC{pc_x}"
    pc_y_col = f"PC{pc_y}"

    if pc_x_col not in scores_df.columns or pc_y_col not in scores_df.columns:
        raise ValueError(
            f"PC columns {pc_x_col} and/or {pc_y_col} not found in DataFrame"
        )

    # Work with a copy
    df = scores_df.copy()

    # Determine color column based on color_by parameter
    if color_by == "ancestry":
        # Reference: known ancestry, Query: predicted (or undefined)
        if show_query_as_undefined:
            df["_plot_color"] = np.where(
                df[SOURCE_COL] == "reference", df[KNOWN_ANCESTRY_COL], undefined_label
            )
        else:
            # Use predicted ancestry for query, known for reference
            df["_plot_color"] = np.where(
                df[SOURCE_COL] == "reference",
                df[KNOWN_ANCESTRY_COL],
                df.get(PREDICTED_ANCESTRY_COL, undefined_label),
            )
    elif color_by == "source":
        df["_plot_color"] = df[SOURCE_COL]
        colors = {"query": undefined_color, "reference": "#1f77b4"}
    elif color_by == "predicted":
        df["_plot_color"] = df.get(PREDICTED_ANCESTRY_COL, undefined_label)
    else:
        raise ValueError(f"Invalid color_by: {color_by}")

    # Apply filters
    if filter_source is not None:
        df = df[df[SOURCE_COL] == filter_source]

    if filter_populations is not None:
        # Filter by either known or predicted ancestry
        mask = df[KNOWN_ANCESTRY_COL].isin(filter_populations) | df.get(
            PREDICTED_ANCESTRY_COL, pd.Series(dtype=str)
        ).isin(filter_populations)
        df = df[mask]

    if len(df) == 0:
        logger.warning("No samples to plot after filtering")
        fig, ax_plot = plt.subplots(figsize=figsize)
        ax_plot.text(0.5, 0.5, "No samples to display", ha="center", va="center")
        return fig

    # Create figure if needed
    if ax is None:
        fig, ax_plot = plt.subplots(figsize=figsize)
    else:
        ax_plot = ax
        fig = ax_plot.figure

    # Track populations and sources for legend
    plotted_populations = set()
    plotted_sources = set()

    # Plot reference samples first (background)
    if shape_by == "source":
        for source in ["reference", "query"]:
            source_df = df[df[SOURCE_COL] == source]
            if len(source_df) == 0:
                continue

            plotted_sources.add(source)
            marker = markers.get(source, "o")
            size = s_reference if source == "reference" else s_query

            for pop in source_df["_plot_color"].unique():
                pop_df = source_df[source_df["_plot_color"] == pop]
                color = colors.get(pop, "#7f7f7f")
                plotted_populations.add(pop)

                ax_plot.scatter(
                    pop_df[pc_x_col],
                    pop_df[pc_y_col],
                    c=color,
                    marker=marker,
                    s=size,
                    alpha=alpha,
                    label=None,  # We'll add custom legend
                    edgecolors="none" if marker == "o" else None,
                    linewidths=1.5 if marker == "x" else None,
                )
    else:
        # No shape distinction
        for pop in df["_plot_color"].unique():
            pop_df = df[df["_plot_color"] == pop]
            color = colors.get(pop, "#7f7f7f")
            plotted_populations.add(pop)

            ax_plot.scatter(
                pop_df[pc_x_col],
                pop_df[pc_y_col],
                c=color,
                s=s_reference,
                alpha=alpha,
                label=None,
            )

    # Add labels
    ax_plot.set_xlabel(pc_x_col, fontsize=12)
    ax_plot.set_ylabel(pc_y_col, fontsize=12)

    if title:
        ax_plot.set_title(title, fontsize=14)
    else:
        ax_plot.set_title(f"{pc_x_col} vs {pc_y_col}", fontsize=14)

    # Create legend
    if show_legend:
        legend_handles = []

        # Population legend (colors)
        pop_order = ["AFR", "AMR", "EAS", "EUR", "SAS", undefined_label, "unassigned"]
        for pop in pop_order:
            if pop in plotted_populations:
                color = colors.get(pop, "#7f7f7f")
                label = POPULATION_NAMES.get(pop, pop)
                if pop == undefined_label:
                    label = undefined_label
                elif pop == "unassigned":
                    label = "Unassigned"
                patch = mpatches.Patch(color=color, label=label)
                legend_handles.append(patch)

        # Add any other populations not in standard order
        for pop in sorted(plotted_populations):
            if pop not in pop_order:
                color = colors.get(pop, "#7f7f7f")
                label = POPULATION_NAMES.get(pop, pop)
                patch = mpatches.Patch(color=color, label=label)
                legend_handles.append(patch)

        # Source legend (shapes) - only if shape_by="source"
        if shape_by == "source" and len(plotted_sources) > 1:
            # Add separator
            legend_handles.append(mpatches.Patch(color="none", label=""))

            for source in ["reference", "query"]:
                if source in plotted_sources:
                    marker = markers.get(source, "o")
                    label = reference_label if source == "reference" else query_label
                    line = Line2D(
                        [0],
                        [0],
                        marker=marker,
                        color="gray",
                        linestyle="None",
                        markersize=8,
                        label=label,
                    )
                    legend_handles.append(line)

        ax_plot.legend(
            handles=legend_handles,
            loc=legend_loc,
            framealpha=0.9,
        )

    plt.tight_layout()
    return fig


def plot_pca_panel(
    scores_df: pd.DataFrame,
    filter_population: Optional[str] = None,
    **kwargs,
) -> Any:
    """Create a two-panel PCA figure (full view + zoomed cluster).

    Similar to the reference plot showing all populations in panel (a)
    and a zoomed view of one cluster in panel (b).

    Parameters
    ----------
    scores_df : pd.DataFrame
        DataFrame with PC scores and ancestry columns.
    filter_population : str, optional
        Population to show in zoomed panel. If None, shows EUR.
    **kwargs
        Additional arguments passed to plot_pca_scatter.

    Returns
    -------
    matplotlib.figure.Figure
        Two-panel figure.
    """
    plt, _, _ = _get_matplotlib()

    if filter_population is None:
        filter_population = "EUR"

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

    # Panel a: Full view
    plot_pca_scatter(
        scores_df,
        ax=ax1,
        title="a) All populations",
        **kwargs,
    )

    # Panel b: Zoomed view
    plot_pca_scatter(
        scores_df,
        ax=ax2,
        filter_populations=[filter_population],
        title=f"b) {POPULATION_NAMES.get(filter_population, filter_population)} cluster",
        **kwargs,
    )

    plt.tight_layout()
    return fig


def plot_variance_explained(
    eigenvalues: List[float],
    n_pcs: Optional[int] = None,
    cumulative: bool = True,
    figsize: Tuple[float, float] = (10, 6),
    ax: Optional[Any] = None,
) -> Any:
    """Plot variance explained by principal components.

    Parameters
    ----------
    eigenvalues : list
        Eigenvalues from PCA.
    n_pcs : int, optional
        Number of PCs to show. If None, shows all.
    cumulative : bool
        Whether to show cumulative variance line. Default: True.
    figsize : tuple
        Figure size. Default: (10, 6).
    ax : matplotlib.axes.Axes, optional
        Axes to plot on.

    Returns
    -------
    matplotlib.figure.Figure
        The figure.
    """
    plt, _, _ = _get_matplotlib()

    total = sum(eigenvalues)
    var_explained = [ev / total for ev in eigenvalues]

    if n_pcs is not None:
        var_explained = var_explained[:n_pcs]

    pcs = list(range(1, len(var_explained) + 1))

    if ax is None:
        fig, ax_plot = plt.subplots(figsize=figsize)
    else:
        ax_plot = ax
        fig = ax_plot.figure

    # Bar plot for individual variance
    bars = ax_plot.bar(
        pcs, var_explained, color="steelblue", alpha=0.7, label="Individual"
    )
    ax_plot.set_xlabel("Principal Component", fontsize=12)
    ax_plot.set_ylabel("Variance Explained", fontsize=12)
    ax_plot.set_title("Variance Explained by Principal Components", fontsize=14)

    # Cumulative line
    if cumulative:
        cumvar = np.cumsum(var_explained)
        ax2 = ax_plot.twinx()
        ax2.plot(pcs, cumvar, "r-o", markersize=5, label="Cumulative")
        ax2.set_ylabel("Cumulative Variance Explained", fontsize=12, color="red")
        ax2.tick_params(axis="y", labelcolor="red")
        ax2.set_ylim(0, 1.05)

        # Combined legend
        lines1, labels1 = ax_plot.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax_plot.legend(lines1 + lines2, labels1 + labels2, loc="center right")

    ax_plot.set_xticks(pcs)

    plt.tight_layout()
    return fig


def plot_ancestry_proportions(
    predictions_df: pd.DataFrame,
    ancestry_col: str = PREDICTED_ANCESTRY_COL,
    source_filter: Optional[str] = "query",
    colors: Optional[Dict[str, str]] = None,
    figsize: Tuple[float, float] = (10, 6),
    ax: Optional[Any] = None,
) -> Any:
    """Bar plot showing proportion of samples per ancestry.

    Parameters
    ----------
    predictions_df : pd.DataFrame
        DataFrame with predictions.
    ancestry_col : str
        Column containing ancestry labels. Default: "predicted_ancestry".
    source_filter : str, optional
        Only include samples from this source ("query", "reference", or None).
        Default: "query".
    colors : dict, optional
        Population -> color mapping.
    figsize : tuple
        Figure size.
    ax : matplotlib.axes.Axes, optional
        Axes to plot on.

    Returns
    -------
    matplotlib.figure.Figure
        The figure.
    """
    plt, _, _ = _get_matplotlib()

    if colors is None:
        colors = SUPERPOP_COLORS.copy()

    df = predictions_df.copy()
    if source_filter and SOURCE_COL in df.columns:
        df = df[df[SOURCE_COL] == source_filter]

    counts = df[ancestry_col].value_counts()
    total = len(df)
    proportions = counts / total

    if ax is None:
        fig, ax_plot = plt.subplots(figsize=figsize)
    else:
        ax_plot = ax
        fig = ax_plot.figure

    # Sort by standard population order
    pop_order = ["AFR", "AMR", "EAS", "EUR", "SAS", "unassigned"]
    sorted_pops = [p for p in pop_order if p in proportions.index]
    sorted_pops += [p for p in proportions.index if p not in pop_order]

    bar_colors = [colors.get(p, "#7f7f7f") for p in sorted_pops]
    bar_values = [proportions[p] for p in sorted_pops]

    bars = ax_plot.bar(sorted_pops, bar_values, color=bar_colors, alpha=0.8)

    # Add count labels on bars
    for bar, pop in zip(bars, sorted_pops):
        count = counts[pop]
        pct = proportions[pop] * 100
        ax_plot.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.01,
            f"{count}\n({pct:.1f}%)",
            ha="center",
            va="bottom",
            fontsize=10,
        )

    ax_plot.set_xlabel("Predicted Ancestry", fontsize=12)
    ax_plot.set_ylabel("Proportion", fontsize=12)
    ax_plot.set_title("Ancestry Distribution", fontsize=14)
    ax_plot.set_ylim(0, max(bar_values) * 1.2)

    plt.tight_layout()
    return fig


def plot_probability_distribution(
    predictions_df: pd.DataFrame,
    prob_col: str = ANCESTRY_PROB_COL,
    source_filter: Optional[str] = "query",
    bins: int = 50,
    figsize: Tuple[float, float] = (10, 6),
    ax: Optional[Any] = None,
) -> Any:
    """Histogram of prediction probabilities.

    Parameters
    ----------
    predictions_df : pd.DataFrame
        DataFrame with predictions.
    prob_col : str
        Column containing probabilities.
    source_filter : str, optional
        Only include samples from this source.
    bins : int
        Number of histogram bins.
    figsize : tuple
        Figure size.
    ax : matplotlib.axes.Axes, optional
        Axes to plot on.

    Returns
    -------
    matplotlib.figure.Figure
        The figure.
    """
    plt, _, _ = _get_matplotlib()

    df = predictions_df.copy()
    if source_filter and SOURCE_COL in df.columns:
        df = df[df[SOURCE_COL] == source_filter]

    probs = df[prob_col].dropna()

    if ax is None:
        fig, ax_plot = plt.subplots(figsize=figsize)
    else:
        ax_plot = ax
        fig = ax_plot.figure

    ax_plot.hist(probs, bins=bins, color="steelblue", alpha=0.7, edgecolor="black")
    ax_plot.axvline(
        x=0.75, color="red", linestyle="--", label="Default threshold (0.75)"
    )

    ax_plot.set_xlabel("Prediction Probability", fontsize=12)
    ax_plot.set_ylabel("Count", fontsize=12)
    ax_plot.set_title("Distribution of Ancestry Prediction Probabilities", fontsize=14)
    ax_plot.legend()

    plt.tight_layout()
    return fig


def plot_confusion_matrix(
    y_true: Union[np.ndarray, pd.Series],
    y_pred: Union[np.ndarray, pd.Series],
    labels: Optional[List[str]] = None,
    normalize: bool = True,
    figsize: Tuple[float, float] = (8, 6),
    cmap: str = "Blues",
    ax: Optional[Any] = None,
) -> Any:
    """Plot confusion matrix heatmap.

    Parameters
    ----------
    y_true : array-like
        True labels.
    y_pred : array-like
        Predicted labels.
    labels : list, optional
        Label order. Auto-detected if None.
    normalize : bool
        Whether to normalize by row (true label). Default: True.
    figsize : tuple
        Figure size.
    cmap : str
        Colormap name.
    ax : matplotlib.axes.Axes, optional
        Axes to plot on.

    Returns
    -------
    matplotlib.figure.Figure
        The figure.
    """
    plt, _, _ = _get_matplotlib()
    try:
        from sklearn.metrics import confusion_matrix as sk_confusion_matrix
    except ImportError:  # pragma: no cover
        raise RuntimeError(
            "Confusion matrix plotting requires scikit-learn. Install it with "
            "'poetry install --extras ancestry' (or --extras ml / --extras psroc)."
        )

    if labels is None:
        labels = sorted(set(y_true) | set(y_pred))

    cm = sk_confusion_matrix(y_true, y_pred, labels=labels)

    if normalize:
        cm = cm.astype("float") / cm.sum(axis=1)[:, np.newaxis]
        cm = np.nan_to_num(cm)

    if ax is None:
        fig, ax_plot = plt.subplots(figsize=figsize)
    else:
        ax_plot = ax
        fig = ax_plot.figure

    im = ax_plot.imshow(cm, interpolation="nearest", cmap=cmap)
    ax_plot.figure.colorbar(im, ax=ax_plot)

    # Show labels
    ax_plot.set(
        xticks=np.arange(len(labels)),
        yticks=np.arange(len(labels)),
        xticklabels=labels,
        yticklabels=labels,
        ylabel="True Ancestry",
        xlabel="Predicted Ancestry",
    )

    # Rotate x labels
    plt.setp(ax_plot.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    # Add text annotations
    fmt = ".2f" if normalize else "d"
    thresh = cm.max() / 2.0
    for i in range(len(labels)):
        for j in range(len(labels)):
            ax_plot.text(
                j,
                i,
                format(cm[i, j], fmt),
                ha="center",
                va="center",
                color="white" if cm[i, j] > thresh else "black",
            )

    title = "Confusion Matrix (Normalized)" if normalize else "Confusion Matrix"
    ax_plot.set_title(title, fontsize=14)

    plt.tight_layout()
    return fig


def encode_figure_to_base64(fig, format: str = "png", dpi: int = 150) -> str:
    """Encode matplotlib figure to a base64 data-URI string for HTML embedding.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure to encode.
    format : str
        Image format ('png', 'svg', etc.).
    dpi : int
        Resolution for raster formats.

    Returns
    -------
    str
        Data-URI string (``data:image/{format};base64,...``) suitable for an
        HTML ``<img src=...>`` attribute.
    """
    # Imported lazily so that importing this module does not require matplotlib
    # (kept optional via ``_get_matplotlib``); ``base`` imports it eagerly.
    from hvantk.algorithms.visualization.base import (
        encode_figure_to_base64 as _encode,
    )

    return _encode(fig, format=format, dpi=dpi, as_data_uri=True)


def close_figure(fig) -> None:
    """Close a matplotlib figure to free memory.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure to close.
    """
    plt, _, _ = _get_matplotlib()
    plt.close(fig)
