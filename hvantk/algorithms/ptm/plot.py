"""PTM-specific visualization functions.

Provides publication-quality plots for PTM-variant analysis results (Q1, Q3).
Follows the pattern of hvantk/enrichex/plot.py: matplotlib-based, optional
seaborn, each function saves to file and returns the Figure.

Plot functions:
    - plot_landscape_summary: P/LP and B/LB counts at PTM site vs proximal vs non-PTM
    - plot_overlap_by_category: P/LP overlap counts per PTM category
    - plot_distance_distribution: P/LP variant distance from nearest PTM site
    - plot_population_af: Mean allele frequency comparison across PTM strata
    - plot_source_overlap: UniProt vs PeptideAtlas site overlap and observation counts
"""

from __future__ import annotations

import logging
from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np

from hvantk.algorithms.visualization.base import save_figure_to_path

try:
    import seaborn as sns
except ImportError:
    sns = None

from hvantk.algorithms.ptm.analysis import PTMLandscapeResult, PTMPopulationResult

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
    figsize: Tuple[int, int] = (13, 5),
    dpi: int = 300,
    title: Optional[str] = None,
    format: str = "png",
) -> plt.Figure:
    """Two-panel landscape plot: variant counts (left) and % P/LP (right).

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
        return _empty_figure(
            output_path, format=format, dpi=dpi, title=title or "PTM-Variant Landscape"
        )

    categories = ["PTM site", "Proximal", "Non-PTM"]
    path_counts = [
        result.n_ptm_site_pathogenic,
        result.n_ptm_proximal_pathogenic,
        result.n_pathogenic
        - result.n_ptm_site_pathogenic
        - result.n_ptm_proximal_pathogenic,
    ]
    benign_counts = [
        result.n_ptm_site_benign,
        result.n_ptm_proximal_benign,
        result.n_benign - result.n_ptm_site_benign - result.n_ptm_proximal_benign,
    ]
    totals = [p + b for p, b in zip(path_counts, benign_counts)]
    pct_plp = [100 * p / t if t > 0 else 0 for p, t in zip(path_counts, totals)]

    if sns is not None:
        sns.set_style("whitegrid")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    x = np.arange(len(categories))
    width = 0.35

    # Left panel: grouped bar chart of counts
    bars_p = ax1.bar(
        x - width / 2,
        path_counts,
        width,
        label="P/LP",
        color=PTM_COLORS["pathogenic"],
        edgecolor="black",
        linewidth=0.5,
    )
    bars_b = ax1.bar(
        x + width / 2,
        benign_counts,
        width,
        label="B/LB",
        color=PTM_COLORS["benign"],
        edgecolor="black",
        linewidth=0.5,
    )

    for bar, n in zip(bars_p, path_counts):
        if n > 0:
            ax1.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height(),
                str(n),
                ha="center",
                va="bottom",
                fontsize=8,
            )
    for bar, n in zip(bars_b, benign_counts):
        if n > 0:
            ax1.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height(),
                str(n),
                ha="center",
                va="bottom",
                fontsize=8,
            )

    ax1.set_xticks(x)
    ax1.set_xticklabels(categories)
    ax1.set_ylabel("Variant count")
    ax1.set_title("Variant Counts")
    ax1.legend(frameon=False)
    ci_hi_str = (
        f"{result.enrichment_ci_high:.2f}" if result.enrichment_ci_high < 1e6 else "∞"
    )
    ax1.text(
        0.98,
        0.95,
        "PTM site + proximal vs non-PTM\n"
        f"OR={result.enrichment_odds_ratio:.2f} "
        f"({result.enrichment_ci_low:.2f}-{ci_hi_str})\n"
        f"p={result.enrichment_p_value:.2e}",
        transform=ax1.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.5),
    )
    ax1.grid(axis="y", linestyle="--", alpha=0.4)

    # Right panel: % P/LP per stratum
    bar_colors = [PTM_COLORS["ptm_site"], PTM_COLORS["proximal"], PTM_COLORS["non_ptm"]]
    bars = ax2.bar(x, pct_plp, color=bar_colors, edgecolor="black", linewidth=0.5)
    for bar, pct, tot in zip(bars, pct_plp, totals):
        ax2.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 1,
            f"{pct:.0f}%\n(n={tot})",
            ha="center",
            va="bottom",
            fontsize=9,
        )
    ax2.set_xticks(x)
    ax2.set_xticklabels(categories)
    ax2.set_ylabel("% Pathogenic (P/LP)")
    ax2.set_title("Pathogenic Proportion")
    y_max = max(pct_plp)
    upper = max(y_max * 1.35, y_max + 3, 1)
    ax2.set_ylim(0, min(upper, 105))
    ax2.axhline(
        y=100 * result.n_pathogenic / max(result.n_pathogenic + result.n_benign, 1),
        color="gray",
        linestyle="--",
        linewidth=1,
        label="Overall rate",
    )
    ax2.legend(frameon=False, fontsize=8)
    ax2.grid(axis="y", linestyle="--", alpha=0.4)

    fig.suptitle(
        title or "PTM-Variant Landscape Summary", fontsize=13, fontweight="bold"
    )
    fig.tight_layout()
    save_figure_to_path(fig, output_path, format=format, dpi=dpi)
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
        return _empty_figure(
            output_path,
            format=format,
            dpi=dpi,
            title=title or "P/LP Variants by PTM Category",
        )

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
    save_figure_to_path(fig, output_path, format=format, dpi=dpi)
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
        return _empty_figure(
            output_path,
            format=format,
            dpi=dpi,
            title=title or "Distance to Nearest PTM Site",
        )

    sorted_dist = sorted(result.distance_distribution.items())
    distances = [d[0] for d in sorted_dist]
    counts = [d[1] for d in sorted_dist]

    if sns is not None:
        sns.set_style("whitegrid")

    fig, ax = plt.subplots(figsize=figsize)

    colors = [
        PTM_COLORS["ptm_site"] if d == 0 else PTM_COLORS["proximal"] for d in distances
    ]
    ax.bar(distances, counts, color=colors, edgecolor="black", linewidth=0.5)
    ax.set_xlabel("Distance to PTM site (residues)")
    ax.set_ylabel("P/LP variant count")
    ax.set_title(title or "P/LP Distance to Nearest PTM Site")
    ax.set_xticks(distances)
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    fig.tight_layout()
    save_figure_to_path(fig, output_path, format=format, dpi=dpi)
    return fig


def plot_population_af(
    result: PTMPopulationResult,
    output_path: str,
    figsize: Tuple[int, int] = (13, 5),
    dpi: int = 300,
    title: Optional[str] = None,
    format: str = "png",
    ultra_rare_threshold: float = 1e-4,
) -> plt.Figure:
    """Two-panel population AF plot: strip plot (left) and % ultra-rare (right).

    Left panel shows individual variant AF values per stratum on a log10 y-axis
    with median lines. Right panel shows % of variants below the ultra-rare
    threshold per stratum.

    Parameters
    ----------
    result : PTMPopulationResult
        Output of ptm_population (must include per-variant AF arrays).
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
    ultra_rare_threshold : float
        AF threshold for "ultra-rare" classification (default: 1e-4).

    Returns
    -------
    matplotlib.figure.Figure
    """
    if result.n_variants == 0:
        return _empty_figure(
            output_path,
            format=format,
            dpi=dpi,
            title=title or "Allele Frequency at PTM Sites",
        )

    labels = ["PTM site", "Proximal", "Non-PTM"]
    af_arrays = [
        np.array(result.ptm_site_afs),
        np.array(result.proximal_afs),
        np.array(result.non_ptm_afs),
    ]
    colors = [PTM_COLORS["ptm_site"], PTM_COLORS["proximal"], PTM_COLORS["non_ptm"]]

    # Fall back to empty figure when AF data is not available (e.g., old JSON)
    if all(len(a) == 0 for a in af_arrays):
        return _empty_figure(
            output_path,
            format=format,
            dpi=dpi,
            title=title or "Allele Frequency at PTM Sites",
            message="Per-variant AF data not available (mean AF shown in overview)",
        )

    if sns is not None:
        sns.set_style("whitegrid")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    # --- Left panel: strip plot on log10 scale ---
    af_floor = 1e-7  # floor for log display of zero/near-zero values
    for i, (afs, color, label) in enumerate(zip(af_arrays, colors, labels)):
        if len(afs) == 0:
            continue
        afs_plot = np.maximum(afs, af_floor)
        jitter = np.random.default_rng(42).uniform(-0.15, 0.15, size=len(afs_plot))
        ax1.scatter(
            np.full(len(afs_plot), i) + jitter,
            afs_plot,
            color=color,
            alpha=0.5,
            s=18,
            edgecolors="none",
            label=label,
        )
        median_raw = float(np.median(afs))
        median_plot = max(median_raw, af_floor)
        ax1.hlines(median_plot, i - 0.3, i + 0.3, color="black", linewidth=2)
        ax1.text(
            i + 0.35, median_plot, f"med={median_raw:.1e}", va="center", fontsize=8
        )

    ax1.set_yscale("log")
    ax1.set_xticks(range(len(labels)))
    ax1.set_xticklabels(labels)
    ax1.set_ylabel("Allele frequency (log scale)")
    ax1.set_title("Per-Variant AF Distribution")
    stratum_counts = [result.n_ptm_site, result.n_ptm_proximal, result.n_non_ptm]
    n_labels = [f"n={n:,}" for n in stratum_counts]
    for i, nl in enumerate(n_labels):
        ax1.text(
            i, ax1.get_ylim()[0], nl, ha="center", va="top", fontsize=8, color="gray"
        )
    ax1.grid(axis="y", linestyle="--", alpha=0.3)

    # --- Right panel: % ultra-rare per stratum ---
    pct_ultra_rare = []
    for afs in af_arrays:
        if len(afs) > 0:
            pct_ultra_rare.append(100 * np.mean(afs < ultra_rare_threshold))
        else:
            pct_ultra_rare.append(0)

    x = np.arange(len(labels))
    bars = ax2.bar(x, pct_ultra_rare, color=colors, edgecolor="black", linewidth=0.5)
    for bar, pct, n_total in zip(bars, pct_ultra_rare, stratum_counts):
        ax2.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 1,
            f"{pct:.0f}%\n(n={n_total:,})",
            ha="center",
            va="bottom",
            fontsize=9,
        )
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels)
    ax2.set_ylabel(f"% variants with AF < {ultra_rare_threshold:.0e}")
    ax2.set_title("Constraint (% Ultra-Rare)")
    ax2.set_ylim(0, min(max(pct_ultra_rare) * 1.3 + 5, 105))
    ax2.grid(axis="y", linestyle="--", alpha=0.4)

    fig.suptitle(
        title or "Population Allele Frequency by PTM Proximity",
        fontsize=13,
        fontweight="bold",
    )
    fig.tight_layout()
    save_figure_to_path(fig, output_path, format=format, dpi=dpi)
    return fig


def plot_source_overlap(
    mapped_tsv: str,
    output_path: str,
    figsize: Tuple[int, int] = (14, 5),
    dpi: int = 300,
    title: Optional[str] = None,
    format: str = "png",
) -> plt.Figure:
    """Three-panel plot comparing UniProt-curated vs PeptideAtlas-observed phospho sites.

    Left panel: site counts by source overlap category (UniProt-only, both, PA-only).
    Center panel: box plot of log10(n_observations) for sites in both sources vs PA-only.
    Right panel: cumulative distribution of PeptideAtlas observation counts, split by
    whether the site is also in UniProt.

    Parameters
    ----------
    mapped_tsv : str
        Path to the combined 13-column mapped TSV (output of ptm build with both sources).
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
    import csv
    from collections import defaultdict

    # Read the combined TSV and group by (uniprot_id, residue_pos)
    site_sources: dict = defaultdict(lambda: {"sources": set(), "n_obs": 0})

    with open(mapped_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            key = (row["uniprot_id"], row["residue_pos"])
            site_sources[key]["sources"].add(row["source_db"])
            n_obs = int(row.get("n_observations", "0"))
            if n_obs > 0:
                site_sources[key]["n_obs"] = max(site_sources[key]["n_obs"], n_obs)

    if not site_sources:
        return _empty_figure(
            output_path,
            format=format,
            dpi=dpi,
            title=title or "Source Overlap",
            message="No sites found in mapped TSV",
        )

    # Categorize sites
    uniprot_only = []
    both_sources = []
    pa_only = []

    for _key, info in site_sources.items():
        has_up = "UniProt" in info["sources"]
        has_pa = "PeptideAtlas" in info["sources"]
        if has_up and has_pa:
            both_sources.append(info["n_obs"])
        elif has_up:
            uniprot_only.append(info["n_obs"])
        elif has_pa:
            pa_only.append(info["n_obs"])

    n_up = len(uniprot_only)
    n_both = len(both_sources)
    n_pa = len(pa_only)

    if sns is not None:
        sns.set_style("whitegrid")

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=figsize)

    # Color scheme
    c_up = "#2196F3"  # blue - UniProt
    c_both = "#9C27B0"  # purple - overlap
    c_pa = "#FF9800"  # orange - PeptideAtlas

    # --- Left panel: site counts by category ---
    categories = ["UniProt\nonly", "Both", "PeptideAtlas\nonly"]
    counts = [n_up, n_both, n_pa]
    colors = [c_up, c_both, c_pa]
    bars = ax1.bar(categories, counts, color=colors, edgecolor="black", linewidth=0.5)
    for bar, n in zip(bars, counts):
        ax1.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            f"{n:,}",
            ha="center",
            va="bottom",
            fontsize=9,
            fontweight="bold",
        )
    ax1.set_ylabel("Number of phospho sites")
    ax1.set_title("Site Overlap")
    ax1.grid(axis="y", linestyle="--", alpha=0.4)

    # --- Center panel: box plot of n_observations ---
    obs_both = np.array(both_sources, dtype=float)
    obs_pa = np.array(pa_only, dtype=float)

    box_data = []
    box_labels = []
    box_colors = []
    if len(obs_both) > 0:
        box_data.append(np.log10(np.maximum(obs_both, 1)))
        box_labels.append(f"Both\n(n={n_both:,})")
        box_colors.append(c_both)
    if len(obs_pa) > 0:
        box_data.append(np.log10(np.maximum(obs_pa, 1)))
        box_labels.append(f"PA only\n(n={n_pa:,})")
        box_colors.append(c_pa)

    if box_data:
        bp = ax2.boxplot(
            box_data,
            labels=box_labels,
            patch_artist=True,
            widths=0.5,
            showfliers=True,
            flierprops=dict(marker=".", markersize=3, alpha=0.3),
        )
        for patch, color in zip(bp["boxes"], box_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)
        for median in bp["medians"]:
            median.set_color("black")
            median.set_linewidth(2)
    else:
        ax2.text(
            0.5,
            0.5,
            "No PeptideAtlas\ndata",
            ha="center",
            va="center",
            transform=ax2.transAxes,
            fontsize=12,
            color="#888888",
        )

    ax2.set_ylabel("log$_{10}$(n_observations)")
    ax2.set_title("PeptideAtlas Evidence")
    ax2.grid(axis="y", linestyle="--", alpha=0.4)

    # --- Right panel: CDF of n_observations ---
    if len(obs_both) > 0:
        sorted_both = np.sort(obs_both)
        cdf_both = np.arange(1, len(sorted_both) + 1) / len(sorted_both)
        ax3.step(sorted_both, cdf_both, color=c_both, linewidth=2, label="Both sources")
    if len(obs_pa) > 0:
        sorted_pa = np.sort(obs_pa)
        cdf_pa = np.arange(1, len(sorted_pa) + 1) / len(sorted_pa)
        ax3.step(sorted_pa, cdf_pa, color=c_pa, linewidth=2, label="PA only")

    ax3.set_xscale("log")
    ax3.set_xlabel("n_observations")
    ax3.set_ylabel("Cumulative fraction")
    ax3.set_title("Observation Count CDF")
    ax3.legend(frameon=False, fontsize=9)
    ax3.grid(True, linestyle="--", alpha=0.4)

    fig.suptitle(
        title or "UniProt vs PeptideAtlas Phospho Site Comparison",
        fontsize=13,
        fontweight="bold",
    )
    fig.tight_layout()
    save_figure_to_path(fig, output_path, format=format, dpi=dpi)
    return fig


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------


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
    ax.text(
        0.5,
        0.5,
        message,
        ha="center",
        va="center",
        fontsize=14,
        color="#888888",
        transform=ax.transAxes,
    )
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.tight_layout()
    if output_path is not None:
        save_figure_to_path(fig, output_path, format=format, dpi=dpi)
    return fig
