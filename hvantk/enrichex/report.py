"""
HTML report generation for EnrichEx visualizations.

The report builder follows the simplified approach documented in
docs/planning/ENRICHEX_VISUALIZATION_DESIGN.md: static matplotlib plots,
inline CSS, and lightweight string formatting (no template engines).
"""

from __future__ import annotations

import html
import math
import logging
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import matplotlib.pyplot as plt
import pandas as pd

from hvantk.enrichex.gene_sets import GeneSetCollection
from hvantk.enrichex.plot import (
    encode_figure_to_base64,
    plot_burden_forest,
    plot_enrichment_dotplot,
)

logger = logging.getLogger(__name__)

_DEFAULT_COLORS = {
    "primary": "#1b2952",
    "secondary": "#4b6cb7",
    "accent": "#4ecdc4",
    "success": "#2ca02c",
    "danger": "#d62728",
}


def generate_report(
    output_path: str,
    overlap_results: Optional[str] = None,
    burden_results: Optional[str] = None,
    gene_sets_path: Optional[str] = None,
    title: str = "EnrichEx Analysis Report",
    description: Optional[str] = None,
    analysis_date: Optional[str] = None,
    analyst_name: Optional[str] = None,
    top_n: int = 20,
    include_methods: bool = True,
    include_gene_lists: bool = True,
    embed_static_plots: bool = True,
) -> None:
    """
    Generate a static HTML report summarizing EnrichEx analyses.

    Parameters
    ----------
    output_path : str
        Destination HTML path. Parent directories will be created automatically.
    overlap_results : Optional[str]
        Path to overlap enrichment TSV file (as produced by the CLI/Python API).
    burden_results : Optional[str]
        Path to burden testing TSV file.
    gene_sets_path : Optional[str]
        Path to GeneSetCollection JSON file (optional metadata section).
    title : str
        Report title shown in the header.
    description : Optional[str]
        Short paragraph describing the analysis context.
    analysis_date : Optional[str]
        Timestamp string. Defaults to the current time if omitted.
    analyst_name : Optional[str]
        Optional analyst attribution shown under the header.
    top_n : int
        Number of rows to highlight in tables/plots.
    include_methods : bool
        Include a short Methods section describing the analyses performed.
    include_gene_lists : bool
        Include <details> blocks with overlapping genes per gene set.
    embed_static_plots : bool
        When True, embed inline base64 PNGs. Otherwise, reference PNGs written
        next to the HTML report.
    """
    if not overlap_results and not burden_results:
        raise ValueError(
            "At least one of overlap_results or burden_results is required."
        )

    output_path = Path(output_path).expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Generating EnrichEx report at %s", output_path)

    overlap_df = _load_results(overlap_results) if overlap_results else None
    burden_df = _load_results(burden_results) if burden_results else None
    gene_sets = GeneSetCollection.load(gene_sets_path) if gene_sets_path else None

    report_data = {
        "title": title,
        "description": description,
        "date": analysis_date or datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "analyst": analyst_name,
        "overview": _create_overview(overlap_df, burden_df),
        "overlap_section": (
            _create_overlap_section(
                overlap_df,
                top_n=top_n,
                include_gene_lists=include_gene_lists,
                embed_static_plots=embed_static_plots,
                output_dir=output_path.parent,
            )
            if overlap_df is not None and not overlap_df.empty
            else None
        ),
        "burden_section": (
            _create_burden_section(
                burden_df,
                top_n=top_n,
                embed_static_plots=embed_static_plots,
                output_dir=output_path.parent,
            )
            if burden_df is not None and not burden_df.empty
            else None
        ),
        "gene_sets_section": (
            _create_gene_sets_section(gene_sets) if gene_sets else None
        ),
        "methods_section": (
            _create_methods_section(
                has_overlap=overlap_df is not None and not overlap_df.empty,
                has_burden=burden_df is not None and not burden_df.empty,
            )
            if include_methods
            else None
        ),
    }

    html_output = _render_report(report_data)
    output_path.write_text(html_output, encoding="utf-8")
    logger.info("Report saved to %s", output_path)


def _load_results(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    logger.info("Loaded %s rows from %s", len(df), path)
    return df


def _create_overview(
    overlap_df: Optional[pd.DataFrame],
    burden_df: Optional[pd.DataFrame],
) -> Optional[Dict[str, Any]]:
    if (overlap_df is None or overlap_df.empty) and (
        burden_df is None or burden_df.empty
    ):
        return None

    overview: Dict[str, Any] = {}
    if overlap_df is not None and not overlap_df.empty:
        p_col = _resolve_pvalue_column(overlap_df)
        top_hits = (
            overlap_df.nsmallest(3, p_col)
            .loc[:, ["gene_set_name", "odds_ratio", p_col]]
            .to_dict(orient="records")
        )
        overview["overlap"] = {
            "total_tests": len(overlap_df),
            "significant": int(_resolve_significance(overlap_df, p_col=p_col).sum()),
            "top_hits": top_hits,
            "effect_label": "Odds Ratio",
            "effect_col": "odds_ratio",
            "p_col": p_col,
        }

    if burden_df is not None and not burden_df.empty:
        p_col = _resolve_pvalue_column(burden_df)
        effect_col = "odds_ratio" if "odds_ratio" in burden_df.columns else "beta"
        top_hits = (
            burden_df.nsmallest(3, p_col)
            .loc[:, ["gene_set_name", effect_col, p_col]]
            .to_dict(orient="records")
        )
        overview["burden"] = {
            "total_tests": len(burden_df),
            "significant": int(_resolve_significance(burden_df, p_col=p_col).sum()),
            "effect_label": (
                "Odds Ratio" if effect_col == "odds_ratio" else "Effect Size"
            ),
            "effect_col": effect_col,
            "top_hits": top_hits,
            "p_col": p_col,
        }

    return overview


def _create_overlap_section(
    df: pd.DataFrame,
    top_n: int,
    include_gene_lists: bool,
    embed_static_plots: bool,
    output_dir: Path,
) -> Dict[str, Any]:
    p_col = _resolve_pvalue_column(df)
    df = df.sort_values(p_col).head(top_n)
    table_rows = []
    for _, row in df.iterrows():
        table_rows.append(
            {
                "gene_set_name": row.get("gene_set_name"),
                "odds_ratio": row.get("odds_ratio"),
                "p_value": row.get("p_value"),
                "p_adjusted": row.get("p_adjusted"),
                "n_overlap": row.get("n_overlap"),
                "significant": bool(row.get("significant", False)),
            }
        )

    plot_path = output_dir / "enrichex_overlap.png"
    fig = plot_enrichment_dotplot(
        df,
        output_path=str(plot_path),
        top_n=top_n,
        sort_by=p_col,
        title="Overlap Enrichment",
    )
    plot_src = plot_path.name
    if embed_static_plots:
        plot_src = f"data:image/png;base64,{encode_figure_to_base64(fig)}"
    plt.close(fig)

    gene_lists = None
    if include_gene_lists and "overlap_genes" in df.columns:
        gene_lists = [
            {
                "gene_set_name": row.get("gene_set_name"),
                "genes": _parse_gene_list(row.get("overlap_genes")),
            }
            for _, row in df.iterrows()
        ]

    return {
        "title": "Overlap Enrichment",
        "plot_src": plot_src,
        "plot_is_inline": embed_static_plots,
        "results": table_rows,
        "gene_lists": gene_lists,
    }


def _create_burden_section(
    df: pd.DataFrame,
    top_n: int,
    embed_static_plots: bool,
    output_dir: Path,
) -> Dict[str, Any]:
    df = df.copy()
    effect_col = "odds_ratio" if "odds_ratio" in df.columns else "beta"
    if "ci_lower" not in df.columns and {"beta", "standard_error"}.issubset(df.columns):
        if effect_col == "odds_ratio":
            df["ci_lower"] = (df["beta"] - 1.96 * df["standard_error"]).map(math.exp)
            df["ci_upper"] = (df["beta"] + 1.96 * df["standard_error"]).map(math.exp)
        else:
            df["ci_lower"] = df["beta"] - 1.96 * df["standard_error"]
            df["ci_upper"] = df["beta"] + 1.96 * df["standard_error"]

    p_col = _resolve_pvalue_column(df)
    df = df.sort_values(p_col).head(top_n)
    table_rows = []
    for _, row in df.iterrows():
        effect_value = row.get(effect_col)
        if effect_col == "odds_ratio":
            if pd.isna(effect_value):
                beta = row.get("beta")
                effect_value = math.exp(beta) if beta is not None and not pd.isna(beta) else beta
        table_rows.append(
            {
                "gene_set_name": row.get("gene_set_name"),
                "effect": effect_value,
                "ci_lower": row.get("ci_lower"),
                "ci_upper": row.get("ci_upper"),
                "p_value": row.get("p_value"),
                "p_adjusted": row.get("p_adjusted"),
                "significant": bool(row.get("significant", False)),
            }
        )

    plot_path = output_dir / "enrichex_burden.png"
    phenotype_type = "binary" if "odds_ratio" in df.columns else "continuous"
    if phenotype_type == "continuous":
        fig = plot_burden_forest(
            df,
            output_path=str(plot_path),
            top_n=top_n,
            phenotype_type=phenotype_type,
            sort_by=p_col,
            ascending=True,
            title="Burden Testing",
        )
    else:
        fig = plot_burden_forest(
            df,
            output_path=str(plot_path),
            top_n=top_n,
            phenotype_type=phenotype_type,
            title="Burden Testing",
        )
    plot_src = plot_path.name
    if embed_static_plots:
        plot_src = f"data:image/png;base64,{encode_figure_to_base64(fig)}"
    plt.close(fig)

    return {
        "title": "Burden Testing",
        "plot_src": plot_src,
        "plot_is_inline": embed_static_plots,
        "effect_label": "Odds Ratio" if phenotype_type == "binary" else "Effect Size",
        "results": table_rows,
    }


def _create_gene_sets_section(collection: GeneSetCollection) -> Dict[str, Any]:
    names = collection.names()
    return {
        "n_sets": len(collection),
        "n_background": len(collection.background_genes),
        "source": collection.source_description or "Custom collection",
        "examples": names[: min(len(names), 5)],
    }


def _create_methods_section(
    has_overlap: bool, has_burden: bool
) -> Optional[Dict[str, str]]:
    if not has_overlap and not has_burden:
        return None
    methods: Dict[str, str] = {}
    if has_overlap:
        methods["overlap"] = (
            "Overlap enrichment was computed using Fisher's exact test with multiple "
            "testing correction applied to adjusted p-values."
        )
    if has_burden:
        methods["burden"] = (
            "Burden testing results were generated using Hail-native regression "
            "(logistic for binary phenotypes or linear for continuous traits)."
        )
    return methods


def _render_report(data: Dict[str, Any]) -> str:
    sections: List[str] = []
    sections.append(
        _build_header(
            title=data["title"],
            description=data["description"],
            date=data["date"],
            analyst=data.get("analyst"),
        )
    )
    if data.get("overview"):
        sections.append(_build_overview(data["overview"]))
    if data.get("overlap_section"):
        sections.append(_build_overlap_section(data["overlap_section"]))
    if data.get("burden_section"):
        sections.append(_build_burden_section(data["burden_section"]))
    if data.get("gene_sets_section"):
        sections.append(_build_gene_sets_section(data["gene_sets_section"]))
    if data.get("methods_section"):
        sections.append(_build_methods_section(data["methods_section"]))
    sections.append(_build_footer())

    css = _get_css_styles(_DEFAULT_COLORS)
    body = "\n".join(section for section in sections if section)
    return f"<!DOCTYPE html><html><head><meta charset='utf-8'><title>{html.escape(data['title'])}</title><style>{css}</style></head><body>{body}</body></html>"


def _build_header(
    title: str,
    description: Optional[str],
    date: str,
    analyst: Optional[str],
) -> str:
    desc_html = (
        f"<p class='description'>{html.escape(description)}</p>" if description else ""
    )
    analyst_html = (
        f"<p class='analyst'>Analyst: {html.escape(analyst)}</p>" if analyst else ""
    )
    return (
        "<header>"
        f"<h1>{html.escape(title)}</h1>"
        f"{desc_html}"
        f"<p class='date'>Generated: {html.escape(date)}</p>"
        f"{analyst_html}"
        "</header>"
    )


def _build_overview(overview: Dict[str, Any]) -> str:
    cards = []
    if "overlap" in overview:
        section = overview["overlap"]
        card = "<div class='card'>"
        card += "<h3>Overlap Enrichment</h3>"
        card += f"<p><strong>{section['significant']}</strong> significant sets"
        card += f" out of {section['total_tests']} tested.</p>"
        if section["top_hits"]:
            card += "<ul>"
            for hit in section["top_hits"]:
                card += (
                    f"<li>{html.escape(str(hit['gene_set_name']))}: "
                    f"{section['effect_label']}={_format_float(hit.get(section['effect_col']))}, "
                    f"p={_format_float(hit.get(section['p_col']))}</li>"
                )
            card += "</ul>"
        card += "</div>"
        cards.append(card)
    if "burden" in overview:
        section = overview["burden"]
        effect_label = section.get("effect_label", "Effect")
        card = "<div class='card'>"
        card += "<h3>Burden Testing</h3>"
        card += f"<p><strong>{section['significant']}</strong> significant sets"
        card += f" out of {section['total_tests']} tested.</p>"
        if section["top_hits"]:
            card += "<ul>"
            for hit in section["top_hits"]:
                card += (
                    f"<li>{html.escape(str(hit['gene_set_name']))}: "
                    f"{effect_label}={_format_float(hit.get(section['effect_col']))}, "
                    f"p={_format_float(hit.get(section['p_col']))}</li>"
                )
            card += "</ul>"
        card += "</div>"
        cards.append(card)
    return f"<section><h2>Overview</h2><div class='card-grid'>{''.join(cards)}</div></section>"


def _build_overlap_section(section: Dict[str, Any]) -> str:
    rows = []
    for record in section["results"]:
        rows.append(
            "<tr>"
            f"<td>{html.escape(str(record['gene_set_name']))}</td>"
            f"<td>{_format_float(record['odds_ratio'])}</td>"
            f"<td>{record['n_overlap'] if record['n_overlap'] is not None else ''}</td>"
            f"<td>{_format_float(record['p_value'])}</td>"
            f"<td>{_format_float(record['p_adjusted'])}</td>"
            f"<td>{'✅' if record['significant'] else ''}</td>"
            "</tr>"
        )

    plot_html = ""
    if section["plot_src"]:
        plot_html = (
            f"<img src='{section['plot_src']}' alt='Overlap Enrichment Plot' "
            "class='embedded-image'/>"
        )

    gene_list_html = ""
    if section.get("gene_lists"):
        blocks = []
        for entry in section["gene_lists"]:
            genes = ", ".join(html.escape(gene) for gene in entry["genes"])
            blocks.append(
                f"<details><summary>{html.escape(entry['gene_set_name'])}</summary><p>{genes}</p></details>"
            )
        gene_list_html = "<div class='gene-lists'>" + "".join(blocks) + "</div>"

    table_html = (
        "<table><thead><tr>"
        "<th>Gene Set</th><th>Odds Ratio</th><th>Overlap</th>"
        "<th>p-value</th><th>Adjusted p</th><th>Significant</th>"
        "</tr></thead><tbody>"
        f"{''.join(rows)}</tbody></table>"
    )

    return (
        "<section id='overlap'>"
        "<h2>Overlap Enrichment</h2>"
        f"{plot_html}"
        f"{table_html}"
        f"{gene_list_html}"
        "</section>"
    )


def _build_burden_section(section: Dict[str, Any]) -> str:
    rows = []
    for record in section["results"]:
        rows.append(
            "<tr>"
            f"<td>{html.escape(str(record['gene_set_name']))}</td>"
            f"<td>{_format_float(record['effect'])}</td>"
            f"<td>{_format_float(record['ci_lower'])} – {_format_float(record['ci_upper'])}</td>"
            f"<td>{_format_float(record['p_value'])}</td>"
            f"<td>{_format_float(record['p_adjusted'])}</td>"
            f"<td>{'✅' if record['significant'] else ''}</td>"
            "</tr>"
        )

    plot_html = ""
    if section["plot_src"]:
        plot_html = f"<img src='{section['plot_src']}' alt='Burden Plot' class='embedded-image'/>"

    table_html = (
        "<table><thead><tr>"
        f"<th>Gene Set</th><th>{html.escape(section['effect_label'])}</th>"
        "<th>95% CI</th><th>p-value</th><th>Adjusted p</th><th>Significant</th>"
        "</tr></thead><tbody>"
        f"{''.join(rows)}</tbody></table>"
    )

    return (
        "<section id='burden'>"
        "<h2>Burden Testing</h2>"
        f"{plot_html}"
        f"{table_html}"
        "</section>"
    )


def _build_gene_sets_section(section: Dict[str, Any]) -> str:
    examples = ""
    if section["examples"]:
        example_list = ", ".join(html.escape(name) for name in section["examples"])
        examples = f"<p><strong>Example sets:</strong> {example_list}</p>"
    return (
        "<section id='gene-sets'>"
        "<h2>Gene Set Collection</h2>"
        f"<p><strong>Total Gene Sets:</strong> {section['n_sets']}</p>"
        f"<p><strong>Background Genes:</strong> {section['n_background']}</p>"
        f"<p><strong>Source:</strong> {html.escape(section['source'])}</p>"
        f"{examples}"
        "</section>"
    )


def _build_methods_section(section: Dict[str, str]) -> str:
    paragraphs = []
    for key, text in section.items():
        paragraphs.append(
            f"<p><strong>{html.escape(key.title())}:</strong> {html.escape(text)}</p>"
        )
    return "<section id='methods'><h2>Methods</h2>" + "".join(paragraphs) + "</section>"


def _build_footer() -> str:
    return "<footer><p>Generated with hvantk EnrichEx.</p></footer>"


def _get_css_styles(colors: Dict[str, str]) -> str:
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
        header h1 {{
            margin: 0 0 10px 0;
        }}
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
        .gene-lists details {{
            margin-top: 10px;
        }}
        footer {{
            text-align: center;
            margin-top: 40px;
            color: #666;
        }}
    """


def _resolve_pvalue_column(df: pd.DataFrame) -> str:
    for column in ("p_adjusted", "p_value"):
        if column in df.columns:
            return column
    raise ValueError("Results must include p_adjusted or p_value columns.")


def _resolve_significance(
    df: pd.DataFrame, p_col: str, threshold: float = 0.05
) -> pd.Series:
    if "significant" in df.columns:
        return df["significant"].astype(bool)
    return df[p_col] < threshold


def _format_float(value: Any) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    if isinstance(value, float):
        if abs(value) >= 0.1:
            return f"{value:.2f}"
        return f"{value:.2e}"
    return str(value)


def _parse_gene_list(value: Any) -> List[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [gene.strip() for gene in value.split(",") if gene.strip()]
    if isinstance(value, list):
        return [str(gene) for gene in value]
    return [str(value)]
