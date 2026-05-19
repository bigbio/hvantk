"""
HTML report generation for QTL cascade analysis.

Produces a self-contained static HTML file with embedded CSS, summary
statistics, gene tables, and base64-encoded plots.
"""

import html as html_mod
import logging
from datetime import datetime
from pathlib import Path
from typing import Optional

import pandas as pd

logger = logging.getLogger(__name__)


def generate_report(
    output_path: str,
    cascade_summary_df: Optional[pd.DataFrame] = None,
    gene_summary_df: Optional[pd.DataFrame] = None,
    coloc_df: Optional[pd.DataFrame] = None,
    class_counts: Optional[dict] = None,
    plot_paths: Optional[dict] = None,
    title: str = "QTL Cascade Analysis Report",
    tissues: Optional[list] = None,
    description: Optional[str] = None,
) -> None:
    """Generate a static HTML report.

    Parameters
    ----------
    output_path : str
        Where to write the HTML file.
    cascade_summary_df : pd.DataFrame, optional
        Variant-level cascade statistics per tissue.
    gene_summary_df : pd.DataFrame, optional
        Gene-level summary table.
    coloc_df : pd.DataFrame, optional
        Colocalization results.
    class_counts : dict, optional
        Cascade class → count mapping.
    plot_paths : dict, optional
        Mapping of plot names to file paths (embedded as base64).
    title : str
        Report title.
    tissues : list[str], optional
        Tissues analysed.
    description : str, optional
        Free-text description.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Generating QTL cascade report at %s", output_path)

    sections = []
    sections.append(
        _overview_section(
            class_counts,
            tissues,
            cascade_summary_df,
        )
    )

    if gene_summary_df is not None and not gene_summary_df.empty:
        sections.append(_gene_table_section(gene_summary_df))

    if coloc_df is not None and not coloc_df.empty:
        sections.append(_coloc_section(coloc_df))

    if plot_paths:
        sections.append(_plots_section(plot_paths))

    sections.append(_methods_section())

    html = _render_html(title, description, sections)
    output_path.write_text(html, encoding="utf-8")
    logger.info("Report saved to %s", output_path)


# ---------------------------------------------------------------------------
# Section builders
# ---------------------------------------------------------------------------


def _overview_section(class_counts, tissues, cascade_summary_df):
    rows = []
    if tissues:
        safe_tissues = html_mod.escape(", ".join(tissues))
        rows.append(f"<tr><td>Tissues</td><td>{safe_tissues}</td></tr>")
    if class_counts:
        total = sum(class_counts.values())
        rows.append(f"<tr><td>Total variant-gene pairs</td><td>{total:,}</td></tr>")
        for cls, cnt in class_counts.items():
            pct = cnt / total * 100 if total else 0
            safe_cls = html_mod.escape(str(cls))
            rows.append(
                f"<tr><td>&nbsp;&nbsp;{safe_cls}</td>"
                f"<td>{cnt:,} ({pct:.1f}%)</td></tr>"
            )
    if cascade_summary_df is not None:
        n_genes = (
            cascade_summary_df["gene_id"].nunique()
            if "gene_id" in cascade_summary_df.columns
            else "?"
        )
        rows.append(f"<tr><td>Unique genes</td><td>{n_genes}</td></tr>")

    table = f"<table class='summary'>{''.join(rows)}</table>" if rows else ""
    return f"<h2>Overview</h2>\n{table}"


def _gene_table_section(df, top_n=50):
    cols_display = [
        c
        for c in [
            "gene_id",
            "gene_symbol",
            "n_concordant",
            "n_eqtl_variants",
            "n_pqtl_variants",
            "best_eqtl_pvalue",
            "best_pqtl_pvalue",
            "oe_lof_upper",
            "has_complete_cascade",
            "coloc_max_h4",
        ]
        if c in df.columns
    ]

    display = (
        df.nlargest(top_n, "n_concordant")
        if "n_concordant" in df.columns
        else df.head(top_n)
    )
    display = display[cols_display]

    header = "".join(f"<th>{html_mod.escape(str(c))}</th>" for c in display.columns)
    body_rows = []
    for _, row in display.iterrows():
        cells = []
        for c in display.columns:
            val = row[c]
            if isinstance(val, float):
                cells.append(f"<td>{val:.4g}</td>")
            else:
                cells.append(f"<td>{html_mod.escape(str(val))}</td>")
        body_rows.append(f"<tr>{''.join(cells)}</tr>")
    body = "".join(body_rows)

    return (
        f"<h2>Top Cascade Genes (by concordant count)</h2>\n"
        f"<table class='genes'><thead><tr>{header}</tr></thead>"
        f"<tbody>{body}</tbody></table>"
    )


def _coloc_section(coloc_df):
    from hvantk.core.qtl_constants import DEFAULT_COLOC_H4_THRESHOLD

    n_total = len(coloc_df)
    n_pass = (coloc_df["H4"] > DEFAULT_COLOC_H4_THRESHOLD).sum()
    median_h4 = coloc_df["H4"].median()

    return (
        f"<h2>Colocalization Results</h2>\n"
        f"<table class='summary'>"
        f"<tr><td>Genes tested</td><td>{n_total}</td></tr>"
        f"<tr><td>P(H4) &gt; {DEFAULT_COLOC_H4_THRESHOLD}</td>"
        f"<td>{n_pass} ({n_pass / n_total * 100:.1f}%)</td></tr>"
        f"<tr><td>Median P(H4)</td><td>{median_h4:.3f}</td></tr>"
        f"</table>"
    )


def _plots_section(plot_paths):
    imgs = []
    for name, path in plot_paths.items():
        p = Path(path)
        if not p.exists():
            continue
        data = p.read_bytes()
        b64 = __import__("base64").b64encode(data).decode("utf-8")
        safe_name = html_mod.escape(name)
        imgs.append(
            f'<div class="plot">'
            f"<h3>{safe_name}</h3>"
            f'<img src="data:image/png;base64,{b64}" alt="{safe_name}"/>'
            f"</div>"
        )
    return f"<h2>Plots</h2>\n{''.join(imgs)}" if imgs else ""


def _methods_section():
    return (
        "<h2>Methods</h2>\n"
        "<p>The QTL cascade joins eQTL and pQTL summary statistics on "
        "<code>(locus, alleles, gene_id)</code> and classifies each "
        "variant-gene pair as <em>eqtl_mediated</em> (concordant direction), "
        "<em>discordant</em>, <em>eqtl_only</em>, or <em>pqtl_only</em>.</p>"
        "<p>Colocalization uses the ABF method "
        "(Giambartolomei et al., 2014, PLoS Genet 10(5):e1004383) "
        "to test whether eQTL and pQTL signals share a causal variant. "
        "P(H4) &gt; 0.8 indicates strong evidence for colocalization.</p>"
        "<p>Gene-level summaries aggregate variant evidence and overlay "
        "LOEUF constraint scores (gnomAD) and disease-gene annotations.</p>"
    )


# ---------------------------------------------------------------------------
# HTML renderer
# ---------------------------------------------------------------------------


def _render_html(title, description, sections):
    date_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    safe_title = html_mod.escape(title) if title else ""
    safe_desc = html_mod.escape(description) if description else ""
    desc_html = f"<p class='desc'>{safe_desc}</p>" if description else ""
    body = "\n".join(sections)

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<title>{safe_title}</title>
<style>
body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI",
        Roboto, Helvetica, sans-serif; max-width: 1100px;
        margin: 2em auto; padding: 0 1em; color: #222; }}
h1 {{ border-bottom: 2px solid #2196F3; padding-bottom: 0.3em; }}
h2 {{ color: #1565C0; margin-top: 2em; }}
table {{ border-collapse: collapse; margin: 1em 0; width: 100%; }}
th, td {{ border: 1px solid #ddd; padding: 6px 10px; text-align: left; }}
th {{ background: #f5f5f5; }}
tr:nth-child(even) {{ background: #fafafa; }}
table.summary {{ width: auto; }}
table.summary td:first-child {{ font-weight: bold; padding-right: 2em; }}
.plot {{ margin: 1.5em 0; }}
.plot img {{ max-width: 100%; border: 1px solid #eee; }}
.desc {{ color: #555; font-style: italic; }}
.footer {{ margin-top: 3em; color: #999; font-size: 0.85em; }}
</style>
</head>
<body>
<h1>{safe_title}</h1>
{desc_html}
{body}
<div class="footer">
Generated on {date_str} by <code>hvantk qtlcascade</code>
</div>
</body>
</html>"""
