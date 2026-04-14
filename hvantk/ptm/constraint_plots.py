"""Visualization panels for PTM constraint analysis (M4).

Renders the four diagnostic panels produced by :mod:`hvantk.ptm.constraint`:

    A. Per-group ranking (forest plot of log2(non-PTM / PTM AF ratio))
    B. τ quartile stratification (depletion ratio by tissue specificity)
    C. τ × LOEUF factorial (2x2 heatmap of log2 ratios)
    D. PTM category × group heatmap (diverging, centered at 0)

This module mirrors the tone of :mod:`hvantk.ptm.plot`: a single public entry
point, lazy imports for matplotlib/seaborn, and empty-data panels are skipped
rather than crashing.
"""

from __future__ import annotations

import logging
import os
from typing import Dict

import matplotlib

if matplotlib.get_backend().lower() not in {"agg", "module://matplotlib_inline.backend_inline"}:
    try:
        matplotlib.use("Agg")
    except Exception:  # pragma: no cover - backend already locked in
        pass

import pandas as pd

logger = logging.getLogger(__name__)

__all__ = ["render_panels"]


_TAU_QUARTILE_ORDER = ["Q1_housekeeping", "Q2", "Q3", "Q4_tissue_specific"]
_TAU_BIN_ORDER = ["broad", "tissue_specific"]
_LOEUF_BIN_ORDER = ["constrained", "unconstrained"]
_MAX_GROUPS_PANEL_A = 25


def render_panels(
    results: Dict[str, pd.DataFrame], output_dir: str
) -> Dict[str, str]:
    """Render the 4 constraint panels as PNGs under ``output_dir/plots/``.

    Parameters
    ----------
    results
        Mapping with keys ``per_group_ranking``, ``tau_quartile``,
        ``loeuf_factorial``, ``category_group_heatmap``.
    output_dir
        Base directory; PNGs are written to ``<output_dir>/plots/``.

    Returns
    -------
    dict
        Mapping panel name -> absolute PNG path. Panels with empty/missing
        input DataFrames are omitted (a warning is logged).
    """
    import matplotlib.pyplot as plt  # noqa: F401

    from hvantk.visualization.base import set_default_style

    set_default_style(style="publication")

    plots_dir = os.path.abspath(os.path.join(output_dir, "plots"))
    os.makedirs(plots_dir, exist_ok=True)

    written: Dict[str, str] = {}

    panel_specs = [
        ("per_group_ranking", _render_per_group_ranking),
        ("tau_quartile", _render_tau_quartile),
        ("loeuf_factorial", _render_loeuf_factorial),
        ("category_group_heatmap", _render_category_group_heatmap),
    ]

    for name, render_fn in panel_specs:
        df = results.get(name)
        if df is None or not isinstance(df, pd.DataFrame) or df.empty:
            logger.warning(
                "Panel '%s' skipped: input DataFrame is empty or missing.", name
            )
            continue

        out_path = os.path.join(plots_dir, f"{name}.png")
        try:
            render_fn(df, out_path)
        except Exception as exc:  # pragma: no cover - defensive
            logger.warning("Failed to render panel '%s': %s", name, exc)
            continue
        logger.info("Wrote %s", out_path)
        written[name] = out_path

    return written


# ---------------------------------------------------------------------------
# Panel A: per-group ranking (forest plot of log2 ratios)
# ---------------------------------------------------------------------------


def _render_per_group_ranking(df: pd.DataFrame, output_path: str) -> None:
    import matplotlib.pyplot as plt

    data = df.dropna(subset=["log2_ratio"]).copy()
    if data.empty:
        logger.warning(
            "per_group_ranking has no finite log2_ratio values; skipping."
        )
        return

    n_total_groups = len(data)
    data = data.sort_values("log2_ratio", ascending=False).head(_MAX_GROUPS_PANEL_A)
    data = data.iloc[::-1]  # barh plots bottom-up; reverse so top bar is largest

    fig_height = max(3.0, 0.28 * len(data) + 1.2)
    fig, ax = plt.subplots(figsize=(7.5, fig_height))

    y = range(len(data))
    colors = ["#d62728" if v >= 0 else "#1f77b4" for v in data["log2_ratio"]]
    ax.barh(
        list(y),
        data["log2_ratio"].to_numpy(),
        color=colors,
        edgecolor="black",
        linewidth=0.4,
    )
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_yticks(list(y))
    ax.set_yticklabels(data["group"].astype(str).tolist(), fontsize=8)
    ax.set_xlabel("log2(non-PTM / PTM AF ratio)")

    title_suffix = (
        f"top {len(data)} of {n_total_groups} groups"
        if n_total_groups > len(data)
        else f"n_groups={n_total_groups}"
    )
    ax.set_title(f"Per-group PTM constraint ranking ({title_suffix})")
    ax.grid(axis="x", linestyle="--", alpha=0.4)

    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Panel B: τ quartile stratification
# ---------------------------------------------------------------------------


def _render_tau_quartile(df: pd.DataFrame, output_path: str) -> None:
    import matplotlib.pyplot as plt

    data = df.copy()
    order = [q for q in _TAU_QUARTILE_ORDER if q in set(data["tau_quartile"].astype(str))]
    extras = [q for q in data["tau_quartile"].astype(str).unique() if q not in order]
    order = order + extras
    data["tau_quartile"] = pd.Categorical(
        data["tau_quartile"].astype(str), categories=order, ordered=True
    )
    data = data.sort_values("tau_quartile")

    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    ratios = data["ratio"].to_numpy()
    log2_vals = data["log2_ratio"].to_numpy()

    x = range(len(data))
    bars = ax.bar(
        list(x),
        ratios,
        color="#4c72b0",
        edgecolor="black",
        linewidth=0.5,
    )
    ax.set_xticks(list(x))
    ax.set_xticklabels(data["tau_quartile"].astype(str).tolist(), rotation=15, ha="right")
    ax.set_ylabel("non-PTM / PTM AF ratio")
    ax.set_title("PTM AF depletion by τ quartile")
    ax.axhline(1.0, color="gray", linestyle="--", linewidth=0.8, alpha=0.7)
    ax.grid(axis="y", linestyle="--", alpha=0.4)

    for bar, log2_v in zip(bars, log2_vals):
        if pd.isna(log2_v):
            continue
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            f"log2={log2_v:.2f}",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Panel C: τ × LOEUF factorial heatmap
# ---------------------------------------------------------------------------


def _render_loeuf_factorial(df: pd.DataFrame, output_path: str) -> None:
    import matplotlib.pyplot as plt

    try:
        import seaborn as sns
    except ImportError:  # pragma: no cover
        logger.warning("seaborn not available; loeuf_factorial panel skipped.")
        return

    data = df.copy()
    pivot = data.pivot_table(
        index="tau_bin",
        columns="loeuf_bin",
        values="log2_ratio",
        aggfunc="first",
    )

    row_order = [r for r in _TAU_BIN_ORDER if r in pivot.index] + [
        r for r in pivot.index if r not in _TAU_BIN_ORDER
    ]
    col_order = [c for c in _LOEUF_BIN_ORDER if c in pivot.columns] + [
        c for c in pivot.columns if c not in _LOEUF_BIN_ORDER
    ]
    pivot = pivot.reindex(index=row_order, columns=col_order)

    vmax = float(pivot.abs().max().max()) if pivot.size else 1.0
    if not (vmax > 0):
        vmax = 1.0

    fig, ax = plt.subplots(figsize=(5.5, 4.0))
    sns.heatmap(
        pivot,
        annot=True,
        fmt=".2f",
        cmap="RdBu_r",
        center=0.0,
        vmin=-vmax,
        vmax=vmax,
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": "log2(non-PTM / PTM AF ratio)"},
        ax=ax,
    )
    ax.set_title("τ × LOEUF factorial (log2 ratio)")
    ax.set_xlabel("LOEUF bin")
    ax.set_ylabel("τ bin")

    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Panel D: PTM category × group heatmap
# ---------------------------------------------------------------------------


def _render_category_group_heatmap(df: pd.DataFrame, output_path: str) -> None:
    import matplotlib.pyplot as plt

    try:
        import seaborn as sns
    except ImportError:  # pragma: no cover
        logger.warning("seaborn not available; category_group_heatmap panel skipped.")
        return

    data = df.copy()
    pivot = data.pivot_table(
        index="category",
        columns="group",
        values="log2_ratio",
        aggfunc="mean",
    )

    if pivot.empty:
        logger.warning("category_group_heatmap pivot is empty; skipping.")
        return

    pivot = pivot.reindex(sorted(pivot.index))
    col_order = pivot.mean(axis=0, skipna=True).sort_values(ascending=False).index.tolist()
    pivot = pivot.reindex(columns=col_order)

    vmax = float(pivot.abs().max().max()) if pivot.size else 1.0
    if not (vmax > 0):
        vmax = 1.0

    n_cols = max(pivot.shape[1], 1)
    n_rows = max(pivot.shape[0], 1)
    fig_width = min(18.0, max(6.0, 0.35 * n_cols + 2.0))
    fig_height = min(14.0, max(3.5, 0.45 * n_rows + 1.5))

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    sns.heatmap(
        pivot,
        annot=False,
        cmap="RdBu_r",
        center=0.0,
        vmin=-vmax,
        vmax=vmax,
        linewidths=0.3,
        linecolor="white",
        cbar_kws={"label": "log2(non-PTM / PTM AF ratio)"},
        ax=ax,
    )
    # Paint NaN cells grey so missing combinations are visually distinct.
    ax.set_facecolor("#dddddd")

    ax.set_xticklabels(ax.get_xticklabels(), rotation=60, ha="right", fontsize=7)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=8)
    ax.set_title("PTM category × group (log2 ratio)")
    ax.set_xlabel("Group")
    ax.set_ylabel("PTM category")

    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
