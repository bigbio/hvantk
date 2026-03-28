"""
Visualization functions for QTL cascade analysis.

All functions return ``matplotlib.figure.Figure`` objects and optionally
save to disk.  Plot styling follows the enrichex/psroc conventions.
"""

import io
import base64
import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from hvantk.qtlcascade.constants import (
    CASCADE_CLASSES,
    CASCADE_CLASS_COLORS,
    CASCADE_CLASS_LABELS,
    DEFAULT_COLOC_H4_THRESHOLD,
)

logger = logging.getLogger(__name__)

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:  # pragma: no cover
    plt = None  # type: ignore[assignment]


def _require_matplotlib():
    if plt is None:
        raise ImportError(
            "matplotlib is required for QTL cascade plots. "
            "Install it with: pip install matplotlib"
        )


def _save_figure(fig, output_path: Optional[str], dpi: int = 300):
    """Save figure to disk if *output_path* is provided."""
    if output_path:
        p = Path(output_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(p), dpi=dpi, bbox_inches="tight")
        logger.info("Saved plot: %s", p)


def encode_figure_to_base64(fig) -> str:
    """Encode a matplotlib figure as a base64 PNG for HTML embedding."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=300, bbox_inches="tight")
    buf.seek(0)
    encoded = base64.b64encode(buf.read()).decode("utf-8")
    buf.close()
    return encoded


# ---------------------------------------------------------------------------
# Cascade class distribution
# ---------------------------------------------------------------------------


def plot_cascade_classes(
    class_counts: dict,
    output_path: Optional[str] = None,
    title: str = "QTL Cascade Classification",
    figsize: tuple = (8, 5),
) -> "plt.Figure":
    """Horizontal bar chart of cascade class counts.

    Parameters
    ----------
    class_counts : dict
        Mapping ``cascade_class`` → count.
    """
    _require_matplotlib()

    classes = [c for c in CASCADE_CLASSES if c in class_counts]
    counts = [class_counts[c] for c in classes]

    if not counts:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(
            0.5,
            0.5,
            "No data",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=14,
        )
        ax.set_title(title, fontsize=14, fontweight="bold")
        _save_figure(fig, output_path)
        return fig

    colors = [CASCADE_CLASS_COLORS[c] for c in classes]
    labels = [CASCADE_CLASS_LABELS[c] for c in classes]

    fig, ax = plt.subplots(figsize=figsize)
    y_pos = np.arange(len(classes))
    ax.barh(y_pos, counts, color=colors, edgecolor="black", linewidth=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=11)
    ax.set_xlabel("Number of variant-gene pairs", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.invert_yaxis()

    # Annotate counts
    for i, v in enumerate(counts):
        ax.text(v + max(counts) * 0.01, i, f"{v:,}", va="center", fontsize=10)

    fig.tight_layout()
    _save_figure(fig, output_path)
    return fig


# ---------------------------------------------------------------------------
# Attenuation scatter
# ---------------------------------------------------------------------------


def plot_attenuation(
    df: pd.DataFrame,
    output_path: Optional[str] = None,
    title: str = "eQTL vs pQTL Effect Sizes",
    figsize: tuple = (7, 7),
) -> "plt.Figure":
    """Scatter plot of eQTL beta vs pQTL beta for concordant/discordant pairs.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain ``eqtl_beta``, ``pqtl_beta``, ``cascade_class``.
    """
    _require_matplotlib()

    both = df[df["cascade_class"].isin(["eqtl_mediated", "discordant"])].copy()
    if both.empty:
        logger.warning("No concordant/discordant pairs for attenuation plot")
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(
            0.5,
            0.5,
            "No data",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=14,
        )
        _save_figure(fig, output_path)
        return fig

    fig, ax = plt.subplots(figsize=figsize)
    for cls in ("eqtl_mediated", "discordant"):
        subset = both[both["cascade_class"] == cls]
        if subset.empty:
            continue
        ax.scatter(
            subset["eqtl_beta"],
            subset["pqtl_beta"],
            c=CASCADE_CLASS_COLORS[cls],
            label=CASCADE_CLASS_LABELS[cls],
            alpha=0.6,
            edgecolor="black",
            linewidth=0.3,
            s=20,
        )

    lim = (
        max(
            abs(both["eqtl_beta"].max()),
            abs(both["pqtl_beta"].max()),
            abs(both["eqtl_beta"].min()),
            abs(both["pqtl_beta"].min()),
        )
        * 1.1
    )
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.axhline(0, color="grey", linewidth=0.5, linestyle="--")
    ax.axvline(0, color="grey", linewidth=0.5, linestyle="--")
    ax.plot(
        [-lim, lim],
        [-lim, lim],
        color="black",
        linewidth=0.8,
        linestyle=":",
        label="y = x",
    )
    ax.set_xlabel("eQTL beta (slope)", fontsize=12)
    ax.set_ylabel("pQTL beta", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.legend(fontsize=9, loc="upper left")
    fig.tight_layout()
    _save_figure(fig, output_path)
    return fig


# ---------------------------------------------------------------------------
# Coloc posterior histogram
# ---------------------------------------------------------------------------


def plot_coloc_posteriors(
    coloc_df: pd.DataFrame,
    output_path: Optional[str] = None,
    title: str = "Colocalization Posterior P(H4)",
    figsize: tuple = (8, 5),
) -> "plt.Figure":
    """Histogram of P(H4) across cascade genes.

    Parameters
    ----------
    coloc_df : pd.DataFrame
        Must contain ``H4`` column.
    """
    _require_matplotlib()

    fig, ax = plt.subplots(figsize=figsize)
    ax.hist(
        coloc_df["H4"],
        bins=30,
        color="#2196F3",
        edgecolor="black",
        linewidth=0.5,
        alpha=0.85,
    )
    ax.axvline(
        DEFAULT_COLOC_H4_THRESHOLD,
        color="#F44336",
        linewidth=1.5,
        linestyle="--",
        label=f"Threshold = {DEFAULT_COLOC_H4_THRESHOLD}",
    )
    n_pass = (coloc_df["H4"] > DEFAULT_COLOC_H4_THRESHOLD).sum()
    ax.set_xlabel("P(H4) — shared causal variant", fontsize=12)
    ax.set_ylabel("Number of genes", fontsize=12)
    ax.set_title(
        f"{title}  ({n_pass}/{len(coloc_df)} genes pass threshold)",
        fontsize=14,
        fontweight="bold",
    )
    ax.legend(fontsize=10)
    fig.tight_layout()
    _save_figure(fig, output_path)
    return fig


# ---------------------------------------------------------------------------
# Cross-tissue heatmap
# ---------------------------------------------------------------------------


def plot_cross_tissue_heatmap(
    gene_tissue_df: pd.DataFrame,
    index_col: str = "gene_id",
    value_col: str = "n_concordant",
    output_path: Optional[str] = None,
    title: str = "Cascade Evidence Across Tissues",
    figsize: tuple = (10, 8),
    top_n: int = 30,
) -> "plt.Figure":
    """Heatmap of rows (genes or cascade classes) × tissues (columns).

    Parameters
    ----------
    gene_tissue_df : pd.DataFrame
        Long-format with columns *index_col*, ``tissue``, and *value_col*.
    index_col : str
        Column to use as row index (e.g. ``gene_id`` or ``cascade_class``).
    value_col : str
        Column to display in the heatmap cells.
    top_n : int
        Show only the top-N rows by total across tissues.
    """
    _require_matplotlib()

    pivot = gene_tissue_df.pivot_table(
        index=index_col,
        columns="tissue",
        values=value_col,
        aggfunc="sum",
        fill_value=0,
    )
    # Select top genes
    pivot["_total"] = pivot.sum(axis=1)
    pivot = pivot.nlargest(top_n, "_total").drop(columns="_total")

    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(pivot.values, aspect="auto", cmap="YlOrRd")
    ax.set_xticks(np.arange(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(np.arange(pivot.shape[0]))
    ax.set_yticklabels(pivot.index, fontsize=8)
    ax.set_title(title, fontsize=14, fontweight="bold")
    fig.colorbar(im, ax=ax, label=value_col, shrink=0.6)
    fig.tight_layout()
    _save_figure(fig, output_path)
    return fig


# ---------------------------------------------------------------------------
# LOEUF boxplots
# ---------------------------------------------------------------------------


def plot_loeuf_by_cascade_class(
    df: pd.DataFrame,
    output_path: Optional[str] = None,
    title: str = "LOEUF by Cascade Class",
    figsize: tuple = (8, 5),
) -> "plt.Figure":
    """Boxplots of LOEUF (oe_lof_upper) stratified by cascade class.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain ``cascade_class`` and ``oe_lof_upper``.
    """
    _require_matplotlib()

    classes = [c for c in CASCADE_CLASSES if c in df["cascade_class"].unique()]
    data = [df.loc[df["cascade_class"] == c, "oe_lof_upper"].dropna() for c in classes]
    colors = [CASCADE_CLASS_COLORS[c] for c in classes]

    fig, ax = plt.subplots(figsize=figsize)
    bp = ax.boxplot(
        data,
        vert=True,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="black", linewidth=1.5),
    )
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax.set_xticklabels(
        [CASCADE_CLASS_LABELS[c] for c in classes],
        rotation=30,
        ha="right",
        fontsize=9,
    )
    ax.set_ylabel("LOEUF (oe_lof_upper)", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    fig.tight_layout()
    _save_figure(fig, output_path)
    return fig
