"""HTML report generation for PTM-variant analysis.

Generates a static HTML report combining landscape (Q1) and population (Q3)
results with embedded matplotlib plots. Follows the enrichex/report.py pattern:
inline CSS, base64 PNGs, lightweight string formatting (no template engines).
"""

from __future__ import annotations

import html
import logging
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import matplotlib.pyplot as plt

from hvantk.ptm.analysis import PTMLandscapeResult, PTMPopulationResult
from hvantk.ptm.plot import (
    encode_figure_to_base64,
    plot_distance_distribution,
    plot_landscape_summary,
    plot_overlap_by_category,
    plot_population_af,
)

logger = logging.getLogger(__name__)

_DEFAULT_COLORS = {
    "primary": "#1b2952",
    "secondary": "#4b6cb7",
    "accent": "#e377c2",
}


def generate_report(
    output_path: str,
    landscape_result: Optional[PTMLandscapeResult] = None,
    population_result: Optional[PTMPopulationResult] = None,
    title: str = "PTM-Variant Analysis Report",
    description: Optional[str] = None,
    analysis_date: Optional[str] = None,
    embed_plots: bool = True,
) -> None:
    """Generate a static HTML report for PTM analysis results.

    Parameters
    ----------
    output_path : str
        Destination HTML file path.
    landscape_result : Optional[PTMLandscapeResult]
        Q1 landscape analysis result.
    population_result : Optional[PTMPopulationResult]
        Q3 population analysis result.
    title : str
        Report title.
    description : Optional[str]
        Short paragraph describing the analysis context.
    analysis_date : Optional[str]
        Timestamp string. Defaults to current time.
    embed_plots : bool
        Embed plots as base64 PNGs (True) or save alongside HTML (False).
    """
    if not landscape_result and not population_result:
        raise ValueError(
            "At least one of landscape_result or population_result is required."
        )

    output_path = Path(output_path).expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Generating PTM report at %s", output_path)

    sections: List[str] = []
    sections.append(_build_header(
        title=title,
        description=description,
        date=analysis_date or datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    ))

    if landscape_result:
        sections.append(_build_landscape_section(
            landscape_result, output_path.parent, embed_plots,
        ))

    if population_result:
        sections.append(_build_population_section(
            population_result, output_path.parent, embed_plots,
        ))

    sections.append(_build_methods_section(
        has_landscape=landscape_result is not None,
        has_population=population_result is not None,
    ))
    sections.append(_build_footer())

    css = _get_css(_DEFAULT_COLORS)
    body = "\n".join(s for s in sections if s)
    html_out = (
        f"<!DOCTYPE html><html><head><meta charset='utf-8'>"
        f"<title>{html.escape(title)}</title>"
        f"<style>{css}</style></head><body>{body}</body></html>"
    )
    output_path.write_text(html_out, encoding="utf-8")
    logger.info("Report saved to %s", output_path)


# ---------------------------------------------------------------------------
# Section builders
# ---------------------------------------------------------------------------


def _build_header(title: str, description: Optional[str], date: str) -> str:
    desc_html = (
        f"<p class='description'>{html.escape(description)}</p>"
        if description else ""
    )
    return (
        "<header>"
        f"<h1>{html.escape(title)}</h1>"
        f"{desc_html}"
        f"<p class='date'>Generated: {html.escape(date)}</p>"
        "</header>"
    )


def _build_landscape_section(
    result: PTMLandscapeResult,
    output_dir: Path,
    embed: bool,
) -> str:
    # Overview card
    overview = (
        "<div class='card-grid'>"
        f"<div class='card'><h3>Variants</h3>"
        f"<p>{result.n_variants:,} total</p>"
        f"<p>{result.n_pathogenic:,} P/LP, {result.n_benign:,} B/LB</p></div>"
        f"<div class='card'><h3>PTM Enrichment</h3>"
        f"<p>OR = {result.enrichment_odds_ratio:.2f}</p>"
        f"<p>p = {result.enrichment_p_value:.2e}</p></div>"
        "</div>"
    )

    # Landscape summary plot
    summary_img = _render_plot(
        plot_landscape_summary, result,
        output_dir / "ptm_landscape_summary.png", embed,
    )

    # Category overlap plot
    category_img = _render_plot(
        plot_overlap_by_category, result,
        output_dir / "ptm_overlap_by_category.png", embed,
    )

    # Distance distribution plot
    distance_img = _render_plot(
        plot_distance_distribution, result,
        output_dir / "ptm_distance_distribution.png", embed,
    )

    # Category table
    table_html = ""
    if result.overlap_by_category:
        rows = ""
        for cat, n in sorted(result.overlap_by_category.items(), key=lambda x: -x[1]):
            rows += f"<tr><td>{html.escape(cat)}</td><td>{n:,}</td></tr>"
        table_html = (
            "<table><thead><tr><th>PTM Category</th>"
            "<th>P/LP Count</th></tr></thead>"
            f"<tbody>{rows}</tbody></table>"
        )

    return (
        "<section><h2>Landscape Analysis (Q1)</h2>"
        f"{overview}{summary_img}{category_img}{table_html}{distance_img}"
        "</section>"
    )


def _build_population_section(
    result: PTMPopulationResult,
    output_dir: Path,
    embed: bool,
) -> str:
    overview = (
        "<div class='card-grid'>"
        f"<div class='card'><h3>Variants</h3>"
        f"<p>{result.n_variants:,} total</p>"
        f"<p>PTM site: {result.n_ptm_site:,}</p>"
        f"<p>Proximal: {result.n_ptm_proximal:,}</p>"
        f"<p>Non-PTM: {result.n_non_ptm:,}</p></div>"
        f"<div class='card'><h3>Mean AF</h3>"
        f"<p>PTM site: {result.mean_af_ptm_site:.2e}</p>"
        f"<p>Proximal: {result.mean_af_ptm_proximal:.2e}</p>"
        f"<p>Non-PTM: {result.mean_af_non_ptm:.2e}</p></div>"
        "</div>"
    )

    af_img = _render_plot(
        plot_population_af, result,
        output_dir / "ptm_population_af.png", embed,
    )

    ccr_html = ""
    if result.ccr_mean_ptm is not None:
        ccr_html = (
            "<div class='card'>"
            "<h3>CCR Comparison</h3>"
            f"<p>Mean CCR at PTM sites: {result.ccr_mean_ptm:.1f}</p>"
            f"<p>Mean CCR at non-PTM: {result.ccr_mean_non_ptm:.1f}</p>"
            "</div>"
        )

    return (
        "<section><h2>Population Analysis (Q3)</h2>"
        f"{overview}{af_img}{ccr_html}"
        "</section>"
    )


def _build_methods_section(has_landscape: bool, has_population: bool) -> str:
    paragraphs = []
    if has_landscape:
        paragraphs.append(
            "<p><strong>Landscape (Q1):</strong> ClinVar P/LP and B/LB variants "
            "were cross-referenced with UniProt PTM sites using position-based "
            "annotation. Enrichment was tested with Fisher's exact test on the "
            "2x2 table of (P/B) x (PTM/non-PTM).</p>"
        )
    if has_population:
        paragraphs.append(
            "<p><strong>Population (Q3):</strong> gnomAD variant allele frequencies "
            "were compared between PTM-site, proximal, and non-PTM coding positions. "
            "Lower mean AF at PTM sites suggests purifying selection.</p>"
        )
    return "<section><h2>Methods</h2>" + "".join(paragraphs) + "</section>"


def _build_footer() -> str:
    return "<footer><p>Generated with hvantk PTM.</p></footer>"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _render_plot(plot_func, result, output_path: Path, embed: bool) -> str:
    """Call a plot function, return an <img> tag (base64 or file reference)."""
    try:
        fig = plot_func(result, str(output_path))
        if embed:
            src = f"data:image/png;base64,{encode_figure_to_base64(fig)}"
        else:
            src = output_path.name
        plt.close(fig)
        return f"<img src='{src}' alt='{output_path.stem}' class='embedded-image'/>"
    except Exception:
        logger.warning("Failed to generate plot %s, skipping.", output_path.stem)
        return ""


def _get_css(colors: Dict[str, str]) -> str:
    return f"""
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            line-height: 1.6;
            color: #333;
            max-width: 960px;
            margin: 0 auto;
            padding: 20px;
            background-color: #f7f7f7;
        }}
        header {{
            background: linear-gradient(90deg, {colors['primary']}, {colors['secondary']});
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 25px;
        }}
        header h1 {{ margin: 0 0 10px 0; }}
        .card-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
            gap: 20px;
            margin-bottom: 20px;
        }}
        .card {{
            background-color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 6px rgba(0,0,0,0.08);
        }}
        section {{
            margin-bottom: 30px;
            background-color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 6px rgba(0,0,0,0.05);
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 15px;
        }}
        th, td {{
            border-bottom: 1px solid #e0e0e0;
            text-align: left;
            padding: 10px;
        }}
        th {{
            background-color: {colors['primary']};
            color: white;
        }}
        .embedded-image {{
            display: block;
            max-width: 100%;
            margin: 15px auto;
            border: 1px solid #eee;
            border-radius: 4px;
        }}
        footer {{
            text-align: center;
            margin-top: 40px;
            color: #666;
        }}
    """
