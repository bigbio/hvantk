"""
Plotting utilities for EnrichEx analyses.

This module provides the publication-quality visualizations described in
docs/planning/ENRICHEX_VISUALIZATION_DESIGN.md. The implementation uses
matplotlib (and seaborn when available) so that no new dependencies are
required beyond the existing hvantk[viz] extra.
"""

from __future__ import annotations

import base64
import io
import logging
from itertools import cycle
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Tuple

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

try:  # seaborn is optional but provides nicer defaults when present
    import seaborn as sns
except ImportError:  # pragma: no cover - exercised when seaborn is unavailable
    sns = None

logger = logging.getLogger(__name__)

_P_VALUE_COLUMNS = ("p_adjusted", "p_value")


def plot_enrichment_dotplot(
    results_df: pd.DataFrame,
    output_path: str,
    figsize: Tuple[int, int] = (10, 6),
    dpi: int = 300,
    top_n: Optional[int] = None,
    min_overlap: int = 1,
    sort_by: str = "p_adjusted",
    color_by: str = "odds_ratio",
    size_by: str = "n_overlap",
    size_range: Tuple[int, int] = (50, 500),
    cmap: str = "RdBu_r",
    alpha: float = 0.85,
    alpha_threshold: float = 0.05,
    show_threshold_line: bool = True,
    show_labels: bool = True,
    label_top_n: int = 5,
    title: Optional[str] = None,
    format: str = "png",
) -> plt.Figure:
    """
    Create a publication-quality enrichment dot plot.

    Parameters
    ----------
    results_df : pd.DataFrame
        Enrichment results from compute_overlap_enrichment_pandas().
        Required columns: gene_set_name, n_overlap, odds_ratio, p_value/p_adjusted
    output_path : str
        Output file path. The suffix will be replaced according to `format`.
    figsize : Tuple[int, int]
        Matplotlib figure size in inches.
    dpi : int
        Resolution when saving raster formats.
    top_n : Optional[int]
        Show the top N entries after sorting. If None, plots all rows.
    min_overlap : int
        Minimum overlap size to plot.
    sort_by : str
        Column used for ordering before selecting `top_n`.
    color_by : str
        Column used for color encoding. Use "significant" to color by p-value.
    size_by : str
        Column used to scale point size.
    size_range : Tuple[int, int]
        Minimum/maximum scatter sizes passed to matplotlib.
    cmap : str
        Matplotlib colormap for numeric color scales.
    alpha : float
        Dot opacity.
    alpha_threshold : float
        Significance threshold drawn as dashed line.
    show_threshold_line : bool
        Draw the horizontal threshold reference line when True.
    show_labels : bool
        Annotate the most significant gene sets.
    label_top_n : int
        Number of top entries (by p-value) to annotate.
    title : Optional[str]
        Plot title. Defaults to "Enrichment Dot Plot".
    format : str
        Export format passed to matplotlib (png, pdf, svg, ...).

    Returns
    -------
    matplotlib.figure.Figure
        The created figure (left open for further customization).
    """
    required_cols = {"gene_set_name", size_by}
    _check_dataframe(results_df, required_cols)
    p_col = _resolve_pvalue_column(results_df)

    df = results_df.copy()
    if "n_overlap" in df.columns:
        df = df[df["n_overlap"] >= max(min_overlap, 0)]
    if df.empty:
        raise ValueError("No enrichment results available to plot.")

    if sort_by not in df.columns:
        raise ValueError(f"Column '{sort_by}' not found in results.")
    ascending = sort_by in {"p_adjusted", "p_value"}
    df = df.sort_values(sort_by, ascending=ascending)
    if top_n is not None:
        df = df.head(top_n)
    if df.empty:
        raise ValueError("Filtering removed all rows to plot.")

    df["plot_order"] = range(len(df))

    if sns is not None:
        sns.set_style("whitegrid")

    y_values = -np.log10(df[p_col].clip(lower=1e-300))
    size_values = _scale_marker_sizes(df[size_by].astype(float).to_numpy(), size_range)
    color_values, color_config = _resolve_color_encoding(
        df, color_by=color_by, alpha_threshold=alpha_threshold, p_col=p_col
    )

    fig, ax = plt.subplots(figsize=figsize)
    scatter_kwargs = {
        "s": size_values,
        "alpha": alpha,
        "edgecolor": "black",
        "linewidth": 0.5,
    }
    if color_config["mode"] == "continuous":
        scatter = ax.scatter(
            df["plot_order"],
            y_values,
            c=color_values,
            cmap=cmap,
            **scatter_kwargs,
        )
        cbar_label = color_config.get("label") or color_by.replace("_", " ").title()
        # Position colorbar vertically on the right side (upper portion)
        cbar = fig.colorbar(
            scatter,
            ax=ax,
            label=cbar_label,
            orientation="vertical",
            pad=0.02,
            aspect=15,
            shrink=0.5,
            anchor=(0.0, 1.0),
        )
        cbar.ax.yaxis.set_label_position("right")
    else:
        scatter = ax.scatter(
            df["plot_order"],
            y_values,
            c=[color_config["mapping"][val] for val in color_values],
            **scatter_kwargs,
        )
        handles = [
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor=color_config["mapping"][val],
                markeredgecolor="black",
                markersize=8,
                label=str(val),
            )
            for val in color_config["order"]
        ]
        color_legend = ax.legend(
            handles=handles,
            title=color_by.replace("_", " ").title(),
            loc="upper left",
            frameon=False,
        )
        ax.add_artist(color_legend)

    ax.set_xticks(df["plot_order"])
    ax.set_xticklabels(df["gene_set_name"], rotation=45, ha="right")
    ax.set_ylabel("-log10(p-value)")
    ax.set_xlabel("Gene set")
    ax.set_title(title or "Enrichment Dot Plot")

    if show_threshold_line:
        threshold = -np.log10(max(alpha_threshold, 1e-300))
        ax.axhline(threshold, color="#666666", linestyle="--", linewidth=1)
        ax.text(
            0.99,
            threshold,
            f"FDR < {alpha_threshold:g}",
            ha="right",
            va="bottom",
            fontsize=10,
            color="#444444",
        )

    if show_labels and label_top_n > 0:
        label_df = df.nsmallest(label_top_n, p_col)
        for _, row in label_df.iterrows():
            ax.text(
                row["plot_order"],
                -np.log10(max(row[p_col], 1e-300)) + 0.1,
                row["gene_set_name"],
                ha="center",
                va="bottom",
                fontsize=9,
                rotation=45,
            )

    _add_size_legend(ax, size_by=size_by, values=df[size_by], size_range=size_range)
    fig.tight_layout()

    _save_figure(fig, output_path, format=format, dpi=dpi)
    return fig


def plot_burden_forest(
    results_df: pd.DataFrame,
    output_path: str,
    figsize: Tuple[int, int] = (9, 8),
    dpi: int = 300,
    phenotype_type: str = "binary",
    top_n: Optional[int] = None,
    sort_by: str = "odds_ratio",
    ascending: bool = False,
    color_by: str = "significant",
    colors: Optional[Dict[str, str]] = None,
    marker_size: int = 80,
    alpha: float = 0.9,
    log_scale: bool = True,
    xlim: Optional[Tuple[float, float]] = None,
    show_values: bool = True,
    show_ci: bool = True,
    show_pvalues: bool = True,
    pvalue_stars: bool = True,
    alpha_threshold: float = 0.05,
    title: Optional[str] = None,
    xlabel: Optional[str] = None,
    format: str = "png",
) -> plt.Figure:
    """
    Create a forest plot summarizing burden test results.

    Parameters mirror the implementation plan documented under section 3.2.
    """
    required_columns = {"gene_set_name", sort_by}
    effect_col = "odds_ratio" if phenotype_type == "binary" else "beta"
    required_columns.update({effect_col, "ci_lower", "ci_upper"})
    p_col = _resolve_pvalue_column(results_df)
    required_columns.add(p_col)
    _check_dataframe(results_df, required_columns)

    df = results_df.copy()
    if phenotype_type not in {"binary", "continuous"}:
        raise ValueError("phenotype_type must be 'binary' or 'continuous'.")
    df = df.sort_values(sort_by, ascending=ascending)
    if top_n is not None:
        df = df.head(top_n)
    if df.empty:
        raise ValueError("No burden results available to plot.")

    if sns is not None:
        sns.set_style("whitegrid")

    df.reset_index(drop=True, inplace=True)
    df["position"] = np.arange(len(df))[::-1]  # plot top row at top

    if color_by == "direction":
        direction_values = np.where(
            df[effect_col] >= (1.0 if phenotype_type == "binary" else 0.0),
            "risk",
            "protective",
        )
        mapping = _categorical_palette(["risk", "protective"])
        color_config = {
            "mode": "categorical",
            "mapping": mapping,
            "order": ["risk", "protective"],
        }
        color_values = direction_values
    else:
        color_values, color_config = _resolve_color_encoding(
            df,
            color_by=color_by,
            alpha_threshold=alpha_threshold,
            p_col=p_col,
        )
    color_config = _apply_custom_colors(color_config, colors, color_by)

    fig, ax = plt.subplots(figsize=figsize)
    if show_ci:
        ax.hlines(
            y=df["position"],
            xmin=df["ci_lower"],
            xmax=df["ci_upper"],
            color="#444444",
            linewidth=1.5,
        )

    scatter_kwargs = {
        "s": marker_size,
        "alpha": alpha,
        "edgecolor": "black",
        "linewidth": 0.5,
    }
    legend_handles: Optional[list] = None
    if color_config["mode"] == "continuous":
        scatter = ax.scatter(
            df[effect_col],
            df["position"],
            c=color_values,
            cmap="viridis",
            **scatter_kwargs,
        )
        fig.colorbar(
            scatter,
            ax=ax,
            label=color_config.get("label") or color_by.replace("_", " ").title(),
        )
    else:
        scatter = ax.scatter(
            df[effect_col],
            df["position"],
            c=[color_config["mapping"][val] for val in color_values],
            **scatter_kwargs,
        )
        legend_handles = [
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor=color_config["mapping"][val],
                markeredgecolor="black",
                label=_format_color_label(
                    color_by=color_by,
                    value=val,
                    alpha_threshold=alpha_threshold,
                ),
            )
            for val in color_config["order"]
        ]

    null_value = 1.0 if phenotype_type == "binary" else 0.0
    ax.axvline(null_value, color="#555555", linestyle="--", linewidth=1)

    ax.set_yticks(df["position"])
    ax.set_yticklabels(df["gene_set_name"])
    ax.set_xlabel(
        xlabel
        or (
            "Odds Ratio (log scale)"
            if phenotype_type == "binary" and log_scale
            else "Effect Size"
        )
    )
    ax.set_title(title or "Burden Forest Plot")
    if log_scale and phenotype_type == "binary":
        ax.set_xscale("log")

    if xlim:
        ax.set_xlim(xlim)

    # Add text annotations outside the plot area on the right side
    if show_values or show_pvalues:
        # Calculate x position outside the plot area
        xlim = ax.get_xlim()
        x_range = xlim[1] - xlim[0]
        text_x = xlim[1] + x_range * 0.05  # Position text 5% outside plot area

        for _, row in df.iterrows():
            text_chunks = []
            if show_values:
                text_chunks.append(
                    f"{effect_col.replace('_', ' ').title()}: {row[effect_col]:.2f} "
                    f"[{row['ci_lower']:.2f}-{row['ci_upper']:.2f}]"
                )
            if show_pvalues:
                stars = _format_pvalue(row[p_col], pvalue_stars)
                text_chunks.append(f"p={row[p_col]:.2e}{stars}")
            if text_chunks:
                ax.text(
                    text_x,
                    row["position"],
                    " | ".join(text_chunks),
                    ha="left",
                    va="center",
                    fontsize=8,
                )

    if legend_handles:
        legend = ax.legend(handles=legend_handles, loc="upper left", frameon=False)
        ax.add_artist(legend)

    fig.tight_layout()
    _save_figure(fig, output_path, format=format, dpi=dpi)
    return fig


def plot_enrichment_barplot(
    results_df: pd.DataFrame,
    output_path: str,
    value: str = "odds_ratio",
    top_n: int = 20,
    orientation: str = "horizontal",
    color_by: str = "significant",
    show_ci: bool = True,
    alpha_threshold: float = 0.05,
    figsize: Tuple[int, int] = (8, 6),
    dpi: int = 300,
    title: Optional[str] = None,
    format: str = "png",
) -> plt.Figure:
    """
    Create a compact bar plot for enrichment results.

    Parameters mirror the simplified design in docs/planning/ENRICHEX_VISUALIZATION_DESIGN.md.
    Use `alpha_threshold` to control the significance cutoff when `color_by="significant"`.
    """
    orientation = orientation.lower()
    if orientation not in {"horizontal", "vertical"}:
        raise ValueError("orientation must be 'horizontal' or 'vertical'.")
    _check_dataframe(results_df, {"gene_set_name"})
    df = results_df.copy()
    p_col = _resolve_pvalue_column(df)
    if p_col not in df.columns:
        raise ValueError("Results must contain p-values or adjusted p-values.")

    if value == "-log10_p":
        df[value] = -np.log10(df[p_col].clip(lower=1e-300))
    elif value not in df.columns:
        raise ValueError(f"Column '{value}' not found in results.")

    df = df.nsmallest(top_n, p_col)
    df = df.sort_values(value, ascending=orientation == "vertical")
    if df.empty:
        raise ValueError("No enrichment results available to plot.")

    color_values, color_config = _resolve_color_encoding(
        df,
        color_by=color_by,
        alpha_threshold=alpha_threshold,
        p_col=p_col,
    )
    if color_config["mode"] == "continuous":
        cmap = plt.get_cmap("viridis")
        values = color_values.astype(float)
        v_min = float(np.nanmin(values))
        v_max = float(np.nanmax(values))
        if np.isclose(v_min, v_max):
            colors = [cmap(0.5)] * len(values)
        else:
            normed = (values - v_min) / (v_max - v_min)
            colors = cmap(normed)
        legend_handles = None
    else:
        colors = [color_config["mapping"][val] for val in color_values]
        legend_handles = [
            Line2D(
                [0],
                [0],
                marker="s",
                color="w",
                markerfacecolor=color_config["mapping"][val],
                markeredgecolor="black",
                label=_format_color_label(
                    color_by=color_by,
                    value=val,
                    alpha_threshold=alpha_threshold,
                ),
            )
            for val in color_config["order"]
        ]

    fig, ax = plt.subplots(figsize=figsize)
    positions = np.arange(len(df))
    label = value.replace("_", " ").title()
    if orientation == "horizontal":
        ax.barh(positions, df[value], color=colors, alpha=0.9)
        ax.set_yticks(positions)
        ax.set_yticklabels(df["gene_set_name"])
        ax.set_xlabel(label)
    else:
        ax.bar(positions, df[value], color=colors, alpha=0.9)
        ax.set_xticks(positions)
        ax.set_xticklabels(df["gene_set_name"], rotation=45, ha="right")
        ax.set_ylabel(label)

    if show_ci and {"ci_lower", "ci_upper"}.issubset(df.columns):
        if value == "-log10_p":
            logger.warning(
                "Cannot plot confidence intervals when value == '-log10_p': CI columns are on the p-value scale and would be mismatched. Skipping error bars."
            )
        else:
            minus = df[value] - df["ci_lower"]
            plus = df["ci_upper"] - df[value]
            errors = np.vstack([minus, plus])
            if orientation == "horizontal":
                ax.errorbar(
                    df[value],
                    positions,
                    xerr=errors,
                    fmt="none",
                    ecolor="#444444",
                    capsize=3,
                )
            else:
                ax.errorbar(
                    positions,
                    df[value],
                    yerr=errors,
                    fmt="none",
                    ecolor="#444444",
                    capsize=3,
                )

    if legend_handles:
        ax.legend(handles=legend_handles, loc="best", frameon=False)

    ax.set_title(title or "Enrichment Bar Plot")
    ax.grid(axis="x" if orientation == "horizontal" else "y", linestyle="--", alpha=0.4)
    fig.tight_layout()

    _save_figure(fig, output_path, format=format, dpi=dpi)
    return fig


def encode_figure_to_base64(
    fig: plt.Figure,
    format: str = "png",
    dpi: int = 200,
) -> str:
    """
    Convert a matplotlib figure to a base64-encoded image string.
    """
    buffer = io.BytesIO()
    fig.savefig(buffer, format=format, dpi=dpi, bbox_inches="tight")
    buffer.seek(0)
    return base64.b64encode(buffer.read()).decode("utf-8")


def _check_dataframe(df: pd.DataFrame, columns: Iterable[str]) -> None:
    missing = [col for col in columns if col not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")


def _resolve_pvalue_column(df: pd.DataFrame) -> str:
    for col in _P_VALUE_COLUMNS:
        if col in df.columns:
            return col
    raise ValueError("Results DataFrame missing p-value columns.")


def _scale_marker_sizes(values: np.ndarray, size_range: Tuple[int, int]) -> np.ndarray:
    min_size, max_size = size_range
    if min_size < 0 or max_size <= 0 or max_size < min_size:
        raise ValueError("Invalid size_range values.")
    v_min = float(np.nanmin(values))
    v_max = float(np.nanmax(values))
    if np.isclose(v_min, v_max):
        return np.full_like(values, (min_size + max_size) / 2.0)
    scaled = (values - v_min) / (v_max - v_min)
    return scaled * (max_size - min_size) + min_size


def _resolve_color_encoding(
    df: pd.DataFrame,
    color_by: str,
    alpha_threshold: float,
    p_col: str,
) -> Tuple[np.ndarray, Dict[str, Any]]:
    if color_by == "significant":
        sig = _resolve_significance(df, p_col=p_col, threshold=alpha_threshold)
        mapping = {True: "#d62728", False: "#1f77b4"}
        return sig.to_numpy(), {
            "mode": "categorical",
            "mapping": mapping,
            "order": [True, False],
        }

    if color_by not in df.columns:
        raise ValueError(f"Column '{color_by}' not found for color encoding.")

    series = df[color_by]
    if pd.api.types.is_numeric_dtype(series):
        return series.to_numpy(), {"mode": "continuous", "label": color_by}

    categories = series.astype(str)
    unique_categories = sorted(categories.unique())
    mapping = _categorical_palette(unique_categories)
    return categories.to_numpy(), {
        "mode": "categorical",
        "mapping": mapping,
        "order": unique_categories,
    }


def _resolve_significance(
    df: pd.DataFrame,
    p_col: str,
    threshold: float = 0.05,
) -> pd.Series:
    if "significant" in df.columns:
        return df["significant"].astype(bool)
    return df[p_col] < threshold


def _format_pvalue(value: float, include_stars: bool) -> str:
    if not include_stars:
        return ""
    if value < 0.001:
        return " ***"
    if value < 0.01:
        return " **"
    if value < 0.05:
        return " *"
    return ""


def _add_size_legend(
    ax: plt.Axes,
    size_by: str,
    values: pd.Series,
    size_range: Tuple[int, int],
) -> None:
    if values.empty:
        return
    min_val = values.min()
    max_val = values.max()
    if min_val == max_val:
        legend_sizes = [min_val]
    else:
        legend_sizes = [min_val, (min_val + max_val) / 2.0, max_val]
    scaled_sizes = _scale_marker_sizes(np.array(legend_sizes), size_range)
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label=f"{val:.0f}",
            markerfacecolor="#bbbbbb",
            markeredgecolor="#666666",
            markersize=np.sqrt(size),
        )
        for val, size in zip(legend_sizes, scaled_sizes)
    ]
    # Position legend outside plot area on the right (lower portion)
    ax.legend(
        handles=handles,
        title=size_by.replace("_", " ").title(),
        loc="center left",
        bbox_to_anchor=(1.02, 0.2),
        frameon=False,
    )


def _save_figure(fig: plt.Figure, output_path: str, format: str, dpi: int) -> None:
    path = Path(output_path)
    if path.suffix.lower() != f".{format.lower()}":
        path = path.with_suffix(f".{format}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", format=format)
    logger.info("Saved figure to %s", path)


def _categorical_palette(categories: Iterable[str]) -> Dict[str, str]:
    unique = list(categories)
    if sns is not None:
        palette = sns.color_palette("Set2", len(unique))
    else:
        palette = plt.rcParams.get("axes.prop_cycle", None)
        palette = palette.by_key().get("color", []) if palette else []
        if not palette:
            palette = ["#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd"]
    mapping = {}
    color_iter = cycle(palette)
    for cat in unique:
        mapping[cat] = next(color_iter)
    return mapping


def _apply_custom_colors(
    color_config: Dict[str, Any],
    custom_colors: Optional[Dict[str, str]],
    color_by: str,
) -> Dict[str, Any]:
    if not custom_colors or color_config["mode"] != "categorical":
        return color_config

    mapping = dict(color_config["mapping"])
    for key, color in custom_colors.items():
        if color_by == "significant":
            if key == "significant":
                mapping[True] = color
            elif key == "not_significant":
                mapping[False] = color
        elif key in mapping:
            mapping[key] = color
    color_config["mapping"] = mapping
    return color_config


def _format_color_label(color_by: str, value: Any, alpha_threshold: float) -> str:
    if color_by == "significant":
        return (
            f"Significant (p < {alpha_threshold:g})"
            if bool(value)
            else "Not significant"
        )
    return str(value)
