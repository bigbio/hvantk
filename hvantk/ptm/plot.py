"""PTM-specific visualization functions.

Provides publication-quality plots for PTM-variant analysis results (Q1, Q3).
Follows the pattern of hvantk/enrichex/plot.py: matplotlib-based, optional
seaborn, each function saves to file and returns the Figure.

Plot functions:
    - plot_landscape_summary: P/LP and B/LB counts at PTM site vs proximal vs non-PTM
    - plot_overlap_by_category: P/LP overlap counts per PTM category
    - plot_distance_distribution: P/LP variant distance from nearest PTM site
    - plot_population_af: Mean allele frequency comparison across PTM strata
"""

from __future__ import annotations

import base64
import io
import logging
from pathlib import Path
from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np

try:
    import seaborn as sns
except ImportError:
    sns = None

from hvantk.ptm.analysis import PTMLandscapeResult, PTMPopulationResult

logger = logging.getLogger(__name__)

PTM_COLORS = {
    "pathogenic": "#D62728",
    "benign": "#1F77B4",
    "ptm_site": "#E377C2",
    "proximal": "#FF7F0E",
    "non_ptm": "#7F7F7F",
}


def plot_landscape_summary(
    result: PTMLandscapeResult,
    output_path: str,
    figsize: Tuple[int, int] = (8, 5),
    dpi: int = 300,
    title: Optional[str] = None,
    format: str = "png",
) -> plt.Figure:
    """Grouped bar chart of P/LP and B/LB counts at PTM site, proximal, and non-PTM.

    Parameters
    ----------
    result : PTMLandscapeResult
        Output of ptm_landscape.
    output_path : str
        Output file path.
    figsize : Tuple[int, int]
        Figure size in inches.
    dpi : int
        Resolution for raster formats.
    title : Optional[str]
        Plot title.
    format : str
        Export format (png, pdf, svg).

    Returns
    -------
    matplotlib.figure.Figure
    """
    if result.n_variants == 0:
        return _empty_figure(output_path, format=format, dpi=dpi,
                             title=title or "PTM-Variant Landscape")

    categories = ["PTM site", "Proximal", "Non-PTM"]
    path_counts = [
        result.n_ptm_site_pathogenic,
        result.n_ptm_proximal_pathogenic,
        result.n_pathogenic - result.n_ptm_site_pathogenic - result.n_ptm_proximal_pathogenic,
    ]
    benign_counts = [
        result.n_ptm_site_benign,
        result.n_ptm_proximal_benign,
        result.n_benign - result.n_ptm_site_benign - result.n_ptm_proximal_benign,
    ]

    if sns is not None:
        sns.set_style("whitegrid")

    x = np.arange(len(categories))
    width = 0.35

    fig, ax = plt.subplots(figsize=figsize)
    ax.bar(x - width / 2, path_counts, width, label="P/LP",
           color=PTM_COLORS["pathogenic"], edgecolor="black", linewidth=0.5)
    ax.bar(x + width / 2, benign_counts, width, label="B/LB",
           color=PTM_COLORS["benign"], edgecolor="black", linewidth=0.5)

    ax.set_xticks(x)
    ax.set_xticklabels(categories)
    ax.set_ylabel("Variant count")
    ax.set_title(title or "PTM-Variant Landscape Summary")
    ax.legend(frameon=False)

    # Annotate enrichment result
    ax.text(
        0.98, 0.95,
        f"Fisher OR={result.enrichment_odds_ratio:.2f}, p={result.enrichment_p_value:.2e}",
        transform=ax.transAxes, ha="right", va="top", fontsize=9,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.5),
    )

    ax.grid(axis="y", linestyle="--", alpha=0.4)
    fig.tight_layout()
    _save_figure(fig, output_path, format=format, dpi=dpi)
    return fig


def plot_overlap_by_category(
    result: PTMLandscapeResult,
    output_path: str,
    figsize: Tuple[int, int] = (8, 5),
    dpi: int = 300,
    title: Optional[str] = None,
    format: str = "png",
) -> plt.Figure:
    """Horizontal bar chart of P/LP variant counts per PTM category.

    Parameters
    ----------
    result : PTMLandscapeResult
        Output of ptm_landscape.
    output_path : str
        Output file path.

    Returns
    -------
    matplotlib.figure.Figure
    """
    if not result.overlap_by_category:
        return _empty_figure(output_path, format=format, dpi=dpi,
                             title=title or "P/LP Variants by PTM Category")

    sorted_cats = sorted(result.overlap_by_category.items(), key=lambda x: x[1])
    labels = [c[0] for c in sorted_cats]
    values = [c[1] for c in sorted_cats]

    if sns is not None:
        sns.set_style("whitegrid")

    fig, ax = plt.subplots(figsize=figsize)
    y = np.arange(len(labels))
    ax.barh(y, values, color=PTM_COLORS["pathogenic"], edgecolor="black", linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel("P/LP variant count")
    ax.set_title(title or "P/LP Variants by PTM Category")
    ax.grid(axis="x", linestyle="--", alpha=0.4)
    fig.tight_layout()
    _save_figure(fig, output_path, format=format, dpi=dpi)
    return fig


def plot_distance_distribution(
    result: PTMLandscapeResult,
    output_path: str,
    figsize: Tuple[int, int] = (8, 5),
    dpi: int = 300,
    title: Optional[str] = None,
    format: str = "png",
) -> plt.Figure:
    """Bar chart of P/LP variant distance (in residues) from nearest PTM site.

    Parameters
    ----------
    result : PTMLandscapeResult
        Output of ptm_landscape.
    output_path : str
        Output file path.

    Returns
    -------
    matplotlib.figure.Figure
    """
    if not result.distance_distribution:
        return _empty_figure(output_path, format=format, dpi=dpi,
                             title=title or "Distance to Nearest PTM Site")

    sorted_dist = sorted(result.distance_distribution.items())
    distances = [d[0] for d in sorted_dist]
    counts = [d[1] for d in sorted_dist]

    if sns is not None:
        sns.set_style("whitegrid")

    fig, ax = plt.subplots(figsize=figsize)

    colors = [PTM_COLORS["ptm_site"] if d == 0 else PTM_COLORS["proximal"]
              for d in distances]
    ax.bar(distances, counts, color=colors, edgecolor="black", linewidth=0.5)
    ax.set_xlabel("Distance to PTM site (residues)")
    ax.set_ylabel("P/LP variant count")
    ax.set_title(title or "P/LP Distance to Nearest PTM Site")
    ax.set_xticks(distances)
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    fig.tight_layout()
    _save_figure(fig, output_path, format=format, dpi=dpi)
    return fig


def plot_population_af(
    result: PTMPopulationResult,
    output_path: str,
    figsize: Tuple[int, int] = (7, 5),
    dpi: int = 300,
    title: Optional[str] = None,
    format: str = "png",
) -> plt.Figure:
    """Grouped bar chart comparing mean allele frequency across PTM strata.

    Parameters
    ----------
    result : PTMPopulationResult
        Output of ptm_population.
    output_path : str
        Output file path.

    Returns
    -------
    matplotlib.figure.Figure
    """
    if result.n_variants == 0:
        return _empty_figure(output_path, format=format, dpi=dpi,
                             title=title or "Allele Frequency at PTM Sites")

    categories = ["PTM site", "Proximal", "Non-PTM"]
    af_values = [
        result.mean_af_ptm_site,
        result.mean_af_ptm_proximal,
        result.mean_af_non_ptm,
    ]
    n_values = [result.n_ptm_site, result.n_ptm_proximal, result.n_non_ptm]
    colors = [PTM_COLORS["ptm_site"], PTM_COLORS["proximal"], PTM_COLORS["non_ptm"]]

    if sns is not None:
        sns.set_style("whitegrid")

    fig, ax = plt.subplots(figsize=figsize)
    x = np.arange(len(categories))
    bars = ax.bar(x, af_values, color=colors, edgecolor="black", linewidth=0.5)

    # Annotate counts on bars
    for bar, n in zip(bars, n_values):
        ax.text(
            bar.get_x() + bar.get_width() / 2, bar.get_height(),
            f"n={n:,}", ha="center", va="bottom", fontsize=9,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(categories)
    ax.set_ylabel("Mean allele frequency")
    ax.set_title(title or "Mean Allele Frequency by PTM Proximity")
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    fig.tight_layout()
    _save_figure(fig, output_path, format=format, dpi=dpi)
    return fig


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------


def encode_figure_to_base64(
    fig: plt.Figure,
    format: str = "png",
    dpi: int = 200,
) -> str:
    """Convert a matplotlib figure to a base64-encoded image string."""
    buffer = io.BytesIO()
    fig.savefig(buffer, format=format, dpi=dpi, bbox_inches="tight")
    buffer.seek(0)
    return base64.b64encode(buffer.read()).decode("utf-8")


def _save_figure(fig: plt.Figure, output_path: str, format: str, dpi: int) -> None:
    path = Path(output_path)
    if path.suffix.lower() != f".{format.lower()}":
        path = path.with_suffix(f".{format}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", format=format)
    logger.info("Saved figure to %s", path)


def _empty_figure(
    output_path: Optional[str] = None,
    format: str = "png",
    dpi: int = 300,
    title: str = "No Data",
    message: str = "No data available",
    figsize: Tuple[int, int] = (8, 4),
) -> plt.Figure:
    """Create a placeholder figure when there is no data to plot."""
    fig, ax = plt.subplots(figsize=figsize)
    ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=14,
            color="#888888", transform=ax.transAxes)
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.tight_layout()
    if output_path is not None:
        _save_figure(fig, output_path, format=format, dpi=dpi)
    return fig
