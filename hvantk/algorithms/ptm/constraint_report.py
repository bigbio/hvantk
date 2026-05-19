"""HTML report generation for PTM constraint analysis.

Builds a self-contained report from :class:`PTMConstraintResult`, embedding
PNG panels produced by :mod:`hvantk.ptm.constraint_plots` as base64 ``<img>``
tags. Follows :mod:`hvantk.ptm.report`: inline CSS, string-built sections,
no template engine.
"""

from __future__ import annotations

import base64
import html
import logging
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Sequence

import pandas as pd

if TYPE_CHECKING:  # pragma: no cover
    from hvantk.algorithms.ptm.constraint import PTMConstraintConfig, PTMConstraintResult

logger = logging.getLogger(__name__)

__all__ = ["render_html"]

_DEFAULT_COLORS = {
    "primary": "#1b2952",
    "secondary": "#4b6cb7",
    "accent": "#e377c2",
}
_NAN_STR = "—"

# Ordered columns per section (only those actually present are rendered).
_COLS_GROUP = (
    "group", "n_ptm", "n_non", "mean_ptm", "mean_non",
    "ratio", "log2_ratio", "p_value",
)
_COLS_TAU = (
    "tau_quartile", "n_ptm", "n_non", "mean_ptm", "mean_non",
    "ratio", "log2_ratio", "p_value",
)
_COLS_LOEUF = (
    "tau_bin", "loeuf_bin", "n_ptm", "n_non", "mean_ptm", "mean_non",
    "ratio", "log2_ratio", "p_value",
)
_COLS_CAT = (
    "group", "category", "n_ptm", "n_non", "mean_ptm", "mean_non",
    "ratio", "log2_ratio", "p_value",
)

_RATIO_COLS = {
    "ratio", "log2_ratio", "mean_ptm", "mean_non",
    "med_ptm", "med_non", "delta", "statistic",
}
_INT_COLS = {"n_ptm", "n_non", "n_genes"}


def render_html(
    result: "PTMConstraintResult",
    config: "PTMConstraintConfig",
    output_dir: str,
    title: str = "PTM Constraint Analysis Report",
) -> str:
    """Write a self-contained HTML report into ``output_dir/report.html``.

    PNG panels are expected at ``output_dir/plots/*.png`` (created by
    :func:`hvantk.ptm.constraint_plots.render_panels`). Missing panels fall
    back to a short notice rather than raising.

    Returns the absolute path to the written HTML file.
    """
    # Lazy import to avoid circular dependency.
    from hvantk.algorithms.ptm.constraint import (  # noqa: F401
        PTMConstraintConfig,
        PTMConstraintResult,
    )

    out_path = Path(output_dir).expanduser().resolve()
    out_path.mkdir(parents=True, exist_ok=True)
    plots_dir = out_path / "plots"
    date_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    sections: List[str] = [
        _build_header(title, date_str, config),
        _build_key_findings(result),
        _build_per_group_section(result, plots_dir),
        _build_tau_section(result, plots_dir),
        _build_loeuf_section(result, plots_dir),
        _build_category_section(result, plots_dir),
        _build_within_gene_section(result),
        "<footer><p>Generated with hvantk PTM constraint.</p></footer>",
    ]

    body = "\n".join(s for s in sections if s)
    html_doc = (
        "<!DOCTYPE html><html><head><meta charset='utf-8'>"
        f"<title>{html.escape(title)}</title>"
        f"<style>{_get_css(_DEFAULT_COLORS)}</style></head>"
        f"<body>{body}</body></html>"
    )

    report_path = out_path / "report.html"
    report_path.write_text(html_doc, encoding="utf-8")
    abs_path = str(report_path.resolve())
    logger.info("Wrote PTM constraint report to %s", abs_path)
    return abs_path


# ---------------------------------------------------------------------------
# Section builders
# ---------------------------------------------------------------------------


def _build_header(
    title: str, date: str, config: "PTMConstraintConfig"
) -> str:
    items = [
        ("grouping", getattr(config, "grouping", "")),
        ("expression source", getattr(config, "expression_source", "")),
        ("label filter", getattr(config, "label_filter", "")),
        ("LOEUF field", getattr(config, "loeuf_field", "")),
        ("min variants/group", getattr(config, "min_variants_per_group", "")),
        ("min cells/group", getattr(config, "min_cells_per_group", "")),
        ("expressed threshold", getattr(config, "expressed_threshold", "")),
        ("flanking codons", getattr(config, "flanking_codons", "")),
    ]
    config_lines = "<br/>".join(
        f"<span class='k'>{html.escape(str(k))}</span>: "
        f"<span class='v'>{html.escape(str(v))}</span>"
        for k, v in items
    )
    return (
        "<header>"
        f"<h1>{html.escape(title)}</h1>"
        f"<p class='date'>Generated: {html.escape(date)}</p>"
        f"<div class='config-echo'>{config_lines}</div>"
        "</header>"
    )


def _build_key_findings(result: "PTMConstraintResult") -> str:
    top_groups = list(getattr(result, "top_groups", []) or [])
    top_card = ""
    if top_groups:
        top = top_groups[0]
        group = html.escape(str(top.get("group", "?")))
        ratio_str = _format_ratio(top.get("ratio"))
        log2_str = _format_ratio(top.get("log2_ratio"))
        top_card = (
            f"<div class='card'><h3>Top group</h3>"
            f"<p><strong>{group}</strong></p>"
            f"<p>ratio = {ratio_str} (log2 = {log2_str})</p></div>"
        )
    return (
        "<section class='key-findings'><h2>Key Findings</h2>"
        "<div class='card-grid'>"
        f"<div class='card'><h3>Variants</h3>"
        f"<p>{int(result.n_variants):,} total</p>"
        f"<p>{int(result.n_variants_ptm):,} PTM · "
        f"{int(result.n_variants_non_ptm):,} non-PTM</p></div>"
        f"<div class='card'><h3>Coverage</h3>"
        f"<p>{int(result.n_groups):,} groups</p>"
        f"<p>{int(result.n_genes):,} genes</p></div>"
        f"{top_card}"
        "</div></section>"
    )


def _build_per_group_section(
    result: "PTMConstraintResult", plots_dir: Path
) -> str:
    rows = list(getattr(result, "top_groups", []) or [])[:20]
    return _build_table_section(
        heading="Per-group ranking",
        note=(
            "Mann-Whitney U per group, sorted by log2 depletion ratio "
            "(non-PTM mean AF / PTM mean AF). Top 20 rows shown."
        ),
        rows=rows,
        columns=_COLS_GROUP,
        image=_embed_png(plots_dir / "per_group_ranking.png", "Per-group ranking"),
    )


def _build_tau_section(result: "PTMConstraintResult", plots_dir: Path) -> str:
    rows = list(getattr(result, "tau_quartile", []) or [])
    return _build_table_section(
        heading="τ quartile stratification",
        note=(
            "Depletion ratios binned by gene-level τ quartile "
            "(Q1 housekeeping → Q4 tissue-specific)."
        ),
        rows=rows,
        columns=_COLS_TAU,
        image=_embed_png(plots_dir / "tau_quartile.png", "tau quartile"),
    )


def _build_loeuf_section(result: "PTMConstraintResult", plots_dir: Path) -> str:
    rows = list(getattr(result, "loeuf_factorial", []) or [])
    return _build_table_section(
        heading="LOEUF × τ factorial",
        note=(
            "Two-way factorial: gene τ bin × LOEUF bin "
            "(constrained vs unconstrained)."
        ),
        rows=rows,
        columns=_COLS_LOEUF,
        image=_embed_png(plots_dir / "loeuf_factorial.png", "LOEUF factorial"),
    )


def _build_category_section(
    result: "PTMConstraintResult", plots_dir: Path
) -> str:
    rows = list(getattr(result, "category_heatmap", []) or [])
    if rows:
        try:
            df = pd.DataFrame(rows)
            if "log2_ratio" in df.columns:
                df = df.assign(_abs=df["log2_ratio"].abs()).sort_values(
                    "_abs", ascending=False, na_position="last"
                ).drop(columns=["_abs"])
            rows = df.head(30).to_dict(orient="records")
        except Exception:  # pragma: no cover - defensive
            logger.warning("Could not sort category heatmap rows.", exc_info=True)
            rows = rows[:30]
    return _build_table_section(
        heading="PTM category × group heatmap",
        note=(
            "Depletion per PTM category within each group; top 30 rows "
            "ranked by |log2 ratio|."
        ),
        rows=rows,
        columns=_COLS_CAT,
        image=_embed_png(
            plots_dir / "category_group_heatmap.png", "Category heatmap"
        ),
    )


def _build_within_gene_section(result: "PTMConstraintResult") -> str:
    s: Dict[str, Any] = dict(getattr(result, "within_gene_paired", {}) or {})
    n_genes_str = _format_int(s.get("n_genes"))
    pct_str = _format_percent(s.get("pct_genes_depleted"))
    cards = (
        "<div class='card-grid'>"
        f"<div class='card'><h3>Genes tested</h3><p>{n_genes_str}</p></div>"
        f"<div class='card'><h3>% genes depleted</h3><p>{pct_str}</p></div>"
        f"<div class='card'><h3>Median Δ (non − PTM)</h3>"
        f"<p>{_format_ratio(s.get('median_delta'))}</p></div>"
        f"<div class='card'><h3>Wilcoxon</h3>"
        f"<p>statistic = {_format_ratio(s.get('statistic'))}</p>"
        f"<p>p = {_format_pvalue(s.get('p_value'))}</p></div>"
        "</div>"
    )
    return (
        "<section><h2>Within-gene paired Wilcoxon</h2>"
        "<p class='note'>Per-gene paired comparison of median gnomAD AF "
        "(PTM vs non-PTM); Wilcoxon signed-rank one-sided "
        "(alternative = 'greater').</p>"
        f"{cards}</section>"
    )


def _build_table_section(
    heading: str,
    note: str,
    rows: Sequence[Dict[str, Any]],
    columns: Sequence[str],
    image: str,
) -> str:
    return (
        f"<section><h2>{html.escape(heading)}</h2>"
        f"<p class='note'>{html.escape(note)}</p>"
        f"{image}{_render_table(rows, columns)}"
        "</section>"
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _embed_png(path: Path, alt: str = "") -> str:
    """Return an ``<img>`` tag with base64-embedded PNG, or a fallback notice."""
    try:
        if not path.exists():
            raise FileNotFoundError(str(path))
        b64 = base64.b64encode(path.read_bytes()).decode("ascii")
        safe_alt = html.escape(alt or path.stem)
        return (
            f"<img src='data:image/png;base64,{b64}' alt='{safe_alt}' "
            "class='embedded-image'/>"
        )
    except Exception:
        logger.warning("Plot not available: %s", path)
        return "<p class='plot-missing'>Plot not available</p>"


def _render_table(rows: Sequence[Dict[str, Any]], columns: Sequence[str]) -> str:
    """Render a list of row dicts as an HTML table restricted to ``columns``."""
    if not rows:
        return "<p class='empty'>No rows to display.</p>"
    try:
        df = pd.DataFrame(list(rows))
    except Exception:  # pragma: no cover
        logger.warning("Could not build DataFrame for table.", exc_info=True)
        return "<p class='empty'>Table unavailable.</p>"

    present = [c for c in columns if c in df.columns] or list(df.columns)
    df = df[present].copy()
    for col in df.columns:
        df[col] = df[col].apply(lambda v, c=col: _format_cell(v, c))

    thead = "".join(f"<th>{html.escape(str(c))}</th>" for c in df.columns)
    body_rows = [
        "<tr>" + "".join(f"<td>{row[c]}</td>" for c in df.columns) + "</tr>"
        for _, row in df.iterrows()
    ]
    return (
        f"<table><thead><tr>{thead}</tr></thead>"
        f"<tbody>{''.join(body_rows)}</tbody></table>"
    )


def _format_cell(value: Any, column: str) -> str:
    """Format a single cell for HTML output with column-aware rules."""
    if _is_nan(value):
        return _NAN_STR
    col = column.lower() if isinstance(column, str) else ""
    if col == "p_value":
        return _format_pvalue(value)
    if col in _RATIO_COLS:
        return _format_ratio(value)
    if col in _INT_COLS:
        return _format_int(value)
    if isinstance(value, float):
        return html.escape(f"{value:.3f}")
    return html.escape(str(value))


def _is_nan(value: Any) -> bool:
    if value is None:
        return True
    try:
        if isinstance(value, float) and pd.isna(value):
            return True
    except Exception:  # pragma: no cover
        return False
    return False


def _format_ratio(value: Any) -> str:
    """Format a floating-point ratio to 3 decimals; NaN → em dash."""
    try:
        if _is_nan(value):
            return _NAN_STR
        fv = float(value)
        return _NAN_STR if pd.isna(fv) else f"{fv:.3f}"
    except (TypeError, ValueError):
        return _NAN_STR


def _format_pvalue(value: Any) -> str:
    """Format a p-value in scientific notation (2 digits after the decimal)."""
    try:
        if _is_nan(value):
            return _NAN_STR
        fv = float(value)
        return _NAN_STR if pd.isna(fv) else f"{fv:.2e}"
    except (TypeError, ValueError):
        return _NAN_STR


def _format_int(value: Any) -> str:
    try:
        if _is_nan(value):
            return _NAN_STR
        return f"{int(value):,}"
    except (TypeError, ValueError):
        return _NAN_STR


def _format_percent(value: Any) -> str:
    try:
        if _is_nan(value):
            return _NAN_STR
        fv = float(value)
        return _NAN_STR if pd.isna(fv) else f"{fv * 100:.1f}%"
    except (TypeError, ValueError):
        return _NAN_STR


def _get_css(colors: Dict[str, str]) -> str:
    p, s, a = colors["primary"], colors["secondary"], colors["accent"]
    return f"""
        body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI',
            Roboto, sans-serif; line-height: 1.55; color: #333;
            max-width: 1060px; margin: 0 auto; padding: 20px;
            background-color: #f7f7f7; }}
        header {{ background: linear-gradient(90deg, {p}, {s}); color: white;
            padding: 26px 30px; border-radius: 10px; margin-bottom: 22px; }}
        header h1 {{ margin: 0 0 8px 0; }}
        header .date {{ margin: 0 0 12px 0; opacity: 0.9; font-size: 0.95em; }}
        header .config-echo {{ font-family: ui-monospace, SFMono-Regular, Menlo,
            Monaco, Consolas, monospace; font-size: 0.85em;
            background: rgba(255,255,255,0.08); padding: 10px 14px;
            border-radius: 6px; }}
        header .config-echo .k {{ opacity: 0.85; }}
        header .config-echo .v {{ font-weight: 600; }}
        section {{ margin-bottom: 26px; background-color: white; padding: 20px;
            border-radius: 8px; box-shadow: 0 2px 6px rgba(0,0,0,0.05); }}
        section h2 {{ margin-top: 0; color: {p};
            border-bottom: 2px solid {s}; padding-bottom: 6px; }}
        .note {{ color: #555; font-size: 0.92em; margin: 4px 0 12px 0; }}
        .empty {{ color: #888; font-style: italic; }}
        .card-grid {{ display: grid;
            grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
            gap: 16px; margin-bottom: 16px; }}
        .card {{ background-color: #fafbff; border: 1px solid #e4e8f3;
            padding: 16px; border-radius: 8px; }}
        .card h3 {{ margin: 0 0 6px 0; font-size: 0.95em; color: {p};
            text-transform: uppercase; letter-spacing: 0.03em; }}
        .card p {{ margin: 2px 0; }}
        table {{ width: 100%; border-collapse: collapse; margin-top: 10px;
            font-size: 0.9em; }}
        th, td {{ border-bottom: 1px solid #e0e0e0; text-align: left;
            padding: 8px 10px; }}
        th {{ background-color: {p}; color: white; }}
        tbody tr:nth-child(even) {{ background-color: #fafafa; }}
        .embedded-image {{ display: block; max-width: 100%; margin: 12px auto;
            border: 1px solid #eee; border-radius: 4px; }}
        .plot-missing {{ color: #a33; font-style: italic; padding: 10px;
            background: #fdf3f3; border-left: 3px solid {a};
            border-radius: 4px; }}
        .key-findings {{ background-color: #eef4ff;
            border-left: 4px solid {s}; }}
        footer {{ text-align: center; margin-top: 30px; color: #666;
            font-size: 0.9em; }}
    """
