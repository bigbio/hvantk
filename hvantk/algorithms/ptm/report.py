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
from typing import Any, Dict, List, Optional, Sequence, TYPE_CHECKING

import matplotlib.pyplot as plt

from hvantk.algorithms.ptm.analysis import PTMLandscapeResult, PTMPopulationResult
from hvantk.algorithms.visualization.base import encode_figure_to_base64
from hvantk.algorithms.ptm.plot import (
    plot_distance_distribution,
    plot_landscape_summary,
    plot_overlap_by_category,
    plot_population_af,
)

if TYPE_CHECKING:  # avoid import-time Hail/statsmodels load for type hints
    from hvantk.algorithms.ptm.atlas import PTMAtlasResult
    from hvantk.algorithms.ptm.lmm import BinnedLMMResult, LMMResult

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
    sections.append(
        _build_header(
            title=title,
            description=description,
            date=analysis_date or datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        )
    )

    sections.append(_build_key_findings(landscape_result, population_result))

    if landscape_result:
        sections.append(
            _build_landscape_section(
                landscape_result,
                output_path.parent,
                embed_plots,
            )
        )

    if population_result:
        sections.append(
            _build_population_section(
                population_result,
                output_path.parent,
                embed_plots,
            )
        )

    sections.append(
        _build_methods_section(
            has_landscape=landscape_result is not None,
            has_population=population_result is not None,
        )
    )
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
        f"<p class='description'>{html.escape(description)}</p>" if description else ""
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
    ci_hi_str = (
        f"{result.enrichment_ci_high:.2f}" if result.enrichment_ci_high < 1e6 else "∞"
    )
    overview = (
        "<div class='card-grid'>"
        f"<div class='card'><h3>Variants</h3>"
        f"<p>{result.n_variants:,} total</p>"
        f"<p>{result.n_pathogenic:,} P/LP, {result.n_benign:,} B/LB</p></div>"
        f"<div class='card'><h3>PTM Enrichment</h3>"
        f"<p>OR = {result.enrichment_odds_ratio:.2f} "
        f"(95% CI: {result.enrichment_ci_low:.2f}–{ci_hi_str})</p>"
        f"<p>p = {result.enrichment_p_value:.2e}</p></div>"
        "</div>"
    )

    # Landscape summary plot
    summary_img = _render_plot(
        plot_landscape_summary,
        result,
        output_dir / "ptm_landscape_summary.png",
        embed,
    )

    # Category overlap plot
    category_img = _render_plot(
        plot_overlap_by_category,
        result,
        output_dir / "ptm_overlap_by_category.png",
        embed,
    )

    # Distance distribution plot
    distance_img = _render_plot(
        plot_distance_distribution,
        result,
        output_dir / "ptm_distance_distribution.png",
        embed,
    )

    # Per-category enrichment table
    table_html = ""
    if result.category_enrichment:
        rows = ""
        for cat, info in sorted(
            result.category_enrichment.items(), key=lambda x: x[1]["p_value"]
        ):
            sig = "*" if info["p_value"] < 0.05 else ""
            ci_lo = info.get("ci_low", 0)
            ci_hi = info.get("ci_high", float("inf"))
            ci_hi_str = f"{ci_hi:.2f}" if ci_hi < 1e6 else "∞"
            rows += (
                f"<tr><td>{html.escape(cat)}</td>"
                f"<td>{info['n_pathogenic']:,}</td>"
                f"<td>{info['n_benign']:,}</td>"
                f"<td>{info['odds_ratio']:.2f} ({ci_lo:.2f}–{ci_hi_str})</td>"
                f"<td>{info['p_value']:.2e} {sig}</td></tr>"
            )
        table_html = (
            "<table><thead><tr><th>PTM Category</th>"
            "<th>P/LP</th><th>B/LB</th>"
            "<th>OR (95% CI)</th><th>p-value</th></tr></thead>"
            f"<tbody>{rows}</tbody></table>"
        )
    elif result.overlap_by_category:
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
        f"{overview}"
        f"{summary_img}{category_img}{table_html}{distance_img}"
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
        plot_population_af,
        result,
        output_dir / "ptm_population_af.png",
        embed,
    )

    return (
        "<section><h2>Population Analysis (Q3)</h2>"
        f"{overview}{af_img}"
        "</section>"
    )


def _build_key_findings(
    landscape: Optional[PTMLandscapeResult],
    population: Optional[PTMPopulationResult],
) -> str:
    """Build a top-level Key Findings summary card."""
    bullets = []
    if landscape and landscape.n_variants > 0:
        sig = (
            "significant" if landscape.enrichment_p_value < 0.05 else "non-significant"
        )
        ci_hi_str = (
            f"{landscape.enrichment_ci_high:.1f}"
            if landscape.enrichment_ci_high < 1e6
            else "∞"
        )
        if landscape.enrichment_odds_ratio > 1:
            effect = f"{landscape.enrichment_odds_ratio:.1f}x enrichment"
        elif 0 < landscape.enrichment_odds_ratio < 1:
            effect = f"{1 / landscape.enrichment_odds_ratio:.1f}x depletion"
        else:
            effect = "no measurable enrichment"
        bullets.append(
            f"Pathogenic variants show <strong>{effect}</strong> at/near PTM sites "
            f"(95% CI: {landscape.enrichment_ci_low:.1f}–{ci_hi_str}; "
            f"{sig}, p={landscape.enrichment_p_value:.1e})."
        )
        if landscape.category_enrichment:
            top_cat = min(
                landscape.category_enrichment.items(),
                key=lambda x: (x[1]["p_value"], -x[1]["odds_ratio"]),
            )
            if top_cat[1]["p_value"] < 0.05:
                bullets.append(
                    f"Strongest per-category signal: <strong>{html.escape(top_cat[0])}</strong> "
                    f"(OR={top_cat[1]['odds_ratio']:.1f}, "
                    f"p={top_cat[1]['p_value']:.1e})."
                )
    if (
        population
        and population.n_ptm_site > 0
        and population.n_non_ptm > 0
        and population.mean_af_non_ptm > 0
    ):
        if population.mean_af_ptm_site > 0:
            if population.mean_af_non_ptm > population.mean_af_ptm_site:
                fold = population.mean_af_non_ptm / population.mean_af_ptm_site
                bullets.append(
                    f"PTM-site variants are <strong>{fold:.1f}x rarer</strong> "
                    f"in the population than non-PTM coding variants "
                    f"(mean AF {population.mean_af_ptm_site:.1e} vs "
                    f"{population.mean_af_non_ptm:.1e})."
                )
            elif population.mean_af_non_ptm < population.mean_af_ptm_site:
                fold = population.mean_af_ptm_site / population.mean_af_non_ptm
                bullets.append(
                    f"PTM-site variants are <strong>{fold:.1f}x more common</strong> "
                    f"in the population than non-PTM coding variants "
                    f"(mean AF {population.mean_af_ptm_site:.1e} vs "
                    f"{population.mean_af_non_ptm:.1e})."
                )
            else:
                bullets.append(
                    f"PTM-site and non-PTM coding variants have similar mean AF "
                    f"({population.mean_af_ptm_site:.1e})."
                )
        else:
            bullets.append(
                f"No non-zero AFs observed at PTM sites "
                f"(non-PTM mean AF {population.mean_af_non_ptm:.1e})."
            )
    if not bullets:
        return ""
    items = "".join(f"<li>{b}</li>" for b in bullets)
    return (
        "<section class='key-findings'>"
        "<h2>Key Findings</h2>"
        f"<ul>{items}</ul>"
        "</section>"
    )


def _build_methods_section(has_landscape: bool, has_population: bool) -> str:
    paragraphs = []
    if has_landscape:
        paragraphs.append(
            "<p><strong>Landscape (Q1):</strong> ClinVar P/LP and B/LB variants "
            "were cross-referenced with UniProt PTM sites mapped to GRCh38 genomic "
            "coordinates via Ensembl GTF (release 113). A variant is classified as "
            "'PTM site' if it overlaps a PTM-modified codon, 'proximal' if within "
            "the flanking window (default: 5 codons / 15 bp), or 'non-PTM' otherwise. "
            "Overall enrichment was tested with Fisher's exact test on the 2×2 table "
            "of (P/B) × (PTM/non-PTM). Per-category enrichment tests each PTM type "
            "against all other PTM-proximal variants.</p>"
        )
    if has_population:
        paragraphs.append(
            "<p><strong>Population (Q3):</strong> gnomAD variant allele frequencies "
            "were compared between PTM-site, proximal, and non-PTM coding positions. "
            "Lower mean AF at PTM sites suggests purifying selection. The '% ultra-rare' "
            "metric shows the fraction of variants with AF &lt; 10<sup>-4</sup>.</p>"
        )
    paragraphs.append(
        "<p><strong>PTM categories:</strong> UniProt MOD_RES descriptions are mapped "
        "to categories (phosphorylation, acetylation, methylation, ubiquitination, "
        "glycosylation, sumoylation) by prefix matching. Unrecognized descriptions "
        "are classified as 'other'.</p>"
    )
    paragraphs.append(
        "<p><strong>Data sources:</strong> PTM sites from UniProt/Swiss-Prot "
        "(reviewed human proteins). Gene models from Ensembl GRCh38 release 113. "
        "Transcript resolution prioritizes MANE Select cross-references.</p>"
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
        logger.warning(
            "Failed to generate plot %s, skipping.", output_path.stem, exc_info=True
        )
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
        .key-findings {{
            background-color: #eef4ff;
            border-left: 4px solid {colors['secondary']};
        }}
        .key-findings ul {{
            margin: 10px 0;
            padding-left: 20px;
        }}
        .key-findings li {{
            margin-bottom: 8px;
        }}
        footer {{
            text-align: center;
            margin-top: 40px;
            color: #666;
        }}
    """


# ---------------------------------------------------------------------------
# Phase-2 report (atlas + SYMBOL annotation + LMM results)
# ---------------------------------------------------------------------------


def generate_phase2_report(
    output_path: str,
    *,
    atlas_result: Optional["PTMAtlasResult"] = None,
    annotation_summary: Optional[Dict[str, Any]] = None,
    lmm_results: Optional[Sequence["LMMResult"]] = None,
    binned_lmm_results: Optional[Sequence["BinnedLMMResult"]] = None,
    title: str = "PTM Phase-2 Analysis Report",
    description: Optional[str] = None,
) -> str:
    """Write a Phase-2 HTML summary (no plots; plots added in Phase 4).

    Each section renders only if its corresponding input is non-None so the
    same report writer serves atlas-only, annotation-only, or test-only runs.

    Parameters
    ----------
    output_path : str
        Destination HTML file path.
    atlas_result : PTMAtlasResult, optional
        Output of :func:`hvantk.ptm.atlas.build_atlas`.
    annotation_summary : dict, optional
        Counters for the SYMBOL-based annotation; expected keys are
        ``n_total``, ``n_ptm_site``, ``n_ptm_proximal``, ``n_both``,
        ``n_neither`` but any subset is accepted.
    lmm_results : sequence of LMMResult, optional
        Per-stratum constraint LMM results (notebook M style).
    binned_lmm_results : sequence of BinnedLMMResult, optional
        Per-stratum binned-interaction LMM results (notebook K style).
    title : str
        Report title.
    description : str, optional
        Short description paragraph rendered beneath the title.

    Returns
    -------
    str
        The ``output_path`` (absolute or relative, whatever the caller passed).
    """
    out_path = Path(output_path).expanduser()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Generating PTM Phase-2 report at %s", out_path)

    sections: List[str] = [
        _build_header(
            title=title,
            description=description,
            date=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        )
    ]

    if atlas_result is not None:
        sections.append(_build_phase2_atlas_section(atlas_result))
    if annotation_summary is not None:
        sections.append(_build_phase2_annotation_section(annotation_summary))
    if lmm_results is not None:
        sections.append(_build_phase2_lmm_section(lmm_results))
    if binned_lmm_results is not None:
        sections.append(_build_phase2_binned_section(binned_lmm_results))

    sections.append(_build_footer())

    css = _get_css(_DEFAULT_COLORS)
    body = "\n".join(s for s in sections if s)
    html_out = (
        f"<!DOCTYPE html><html><head><meta charset='utf-8'>"
        f"<title>{html.escape(title)}</title>"
        f"<style>{css}</style></head><body>{body}</body></html>"
    )
    out_path.write_text(html_out, encoding="utf-8")
    logger.info("Phase-2 report saved to %s", out_path)
    return str(output_path)


def _build_phase2_atlas_section(result: "PTMAtlasResult") -> str:
    sources = ", ".join(html.escape(s) for s in result.sources_used) or "-"
    return (
        "<section><h2>Atlas Summary</h2>"
        "<div class='card-grid'>"
        f"<div class='card'><h3>Sources</h3><p>{sources}</p></div>"
        f"<div class='card'><h3>Sites Mapped</h3>"
        f"<p>{result.n_sites:,}</p></div>"
        "</div>"
        "<table><thead><tr><th>Output</th><th>Path</th></tr></thead>"
        "<tbody>"
        f"<tr><td>Combined TSV</td><td>{html.escape(result.combined_tsv)}</td></tr>"
        f"<tr><td>Hail Table</td><td>{html.escape(result.output_ht)}</td></tr>"
        "</tbody></table>"
        "</section>"
    )


def _build_phase2_annotation_section(summary: Dict[str, Any]) -> str:
    keys = ("n_total", "n_ptm_site", "n_ptm_proximal", "n_both", "n_neither")
    rows = ""
    for k in keys:
        if k not in summary:
            continue
        val = summary[k]
        try:
            val_str = f"{int(val):,}"
        except (TypeError, ValueError):
            val_str = html.escape(str(val))
        rows += f"<tr><td>{html.escape(k)}</td><td>{val_str}</td></tr>"
    # Include any extra keys the caller supplied.
    extras = {k: v for k, v in summary.items() if k not in keys}
    for k, v in extras.items():
        try:
            v_str = f"{int(v):,}"
        except (TypeError, ValueError):
            v_str = html.escape(str(v))
        rows += f"<tr><td>{html.escape(str(k))}</td><td>{v_str}</td></tr>"
    if not rows:
        rows = "<tr><td colspan='2'>(no counts provided)</td></tr>"
    return (
        "<section><h2>SYMBOL Annotation Summary</h2>"
        "<table><thead><tr><th>Metric</th><th>Count</th></tr></thead>"
        f"<tbody>{rows}</tbody></table>"
        "</section>"
    )


def _fmt_float(val: float, fmt: str = "{:.3g}") -> str:
    try:
        if val is None:
            return "-"
        fv = float(val)
    except (TypeError, ValueError):
        return html.escape(str(val))
    if fv != fv:  # NaN check
        return "NaN"
    return fmt.format(fv)


def _build_phase2_lmm_section(results: Sequence["LMMResult"]) -> str:
    rows = ""
    for r in results:
        rows += (
            "<tr>"
            f"<td>{html.escape(str(r.stratum))}</td>"
            f"<td>{r.n_variants:,}</td>"
            f"<td>{r.n_ptm:,}</td>"
            f"<td>{r.n_nonptm:,}</td>"
            f"<td>{r.n_genes:,}</td>"
            f"<td>{r.n_mixed_genes:,}</td>"
            f"<td>{_fmt_float(r.beta_ptm)}</td>"
            f"<td>{_fmt_float(r.se_ptm)}</td>"
            f"<td>{_fmt_float(r.p_ptm, '{:.2e}')}</td>"
            f"<td>{'yes' if r.converged else 'no'}</td>"
            f"<td>{html.escape(r.note or '')}</td>"
            "</tr>"
        )
    if not rows:
        rows = "<tr><td colspan='11'>(no results)</td></tr>"
    return (
        "<section><h2>Constraint LMM</h2>"
        "<p>Per-stratum <code>log_af ~ is_ptm + (1|gene)</code> "
        "(notebook M).</p>"
        "<table><thead><tr>"
        "<th>Stratum</th><th>N</th><th>N PTM</th><th>N non-PTM</th>"
        "<th>N genes</th><th>N mixed</th>"
        "<th>beta</th><th>SE</th><th>p</th><th>converged</th><th>note</th>"
        f"</tr></thead><tbody>{rows}</tbody></table>"
        "</section>"
    )


def _build_phase2_binned_section(results: Sequence["BinnedLMMResult"]) -> str:
    # Union of bin labels across strata (preserve first-seen order).
    all_bins: List[str] = []
    for r in results:
        for b in r.bin_levels:
            if b not in all_bins:
                all_bins.append(b)

    header = "<tr><th>Stratum</th><th>N</th><th>N genes</th>"
    for b in all_bins:
        header += f"<th>{html.escape(b)} beta</th><th>SE</th><th>p</th>"
    header += "<th>converged</th><th>note</th></tr>"

    body_rows = ""
    for r in results:
        row = (
            f"<tr><td>{html.escape(str(r.stratum))}</td>"
            f"<td>{r.n_variants:,}</td>"
            f"<td>{r.n_genes:,}</td>"
        )
        for b in all_bins:
            if b in r.bin_betas:
                row += (
                    f"<td>{_fmt_float(r.bin_betas[b])}</td>"
                    f"<td>{_fmt_float(r.bin_ses.get(b, float('nan')))}</td>"
                    f"<td>{_fmt_float(r.bin_pvalues.get(b, float('nan')), '{:.2e}')}</td>"
                )
            else:
                row += "<td>-</td><td>-</td><td>-</td>"
        row += (
            f"<td>{'yes' if r.converged else 'no'}</td>"
            f"<td>{html.escape(r.note or '')}</td></tr>"
        )
        body_rows += row

    if not body_rows:
        body_rows = (
            f"<tr><td colspan='{3 + 3 * len(all_bins) + 2}'>" "(no results)</td></tr>"
        )

    return (
        "<section><h2>Binned-Interaction LMM</h2>"
        "<p>Per-stratum <code>log_af ~ is_ptm * C(expr_bin) + (1|gene)</code> "
        "(notebook K). Reference bin <code>b0_none</code> is zero-expression; "
        "<code>b1..bK</code> are quantiles of <code>log2(expr + 1)</code>. "
        "Missing cells indicate bins not realized for that stratum.</p>"
        f"<table><thead>{header}</thead><tbody>{body_rows}</tbody></table>"
        "</section>"
    )
