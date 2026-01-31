"""HTML report generation for ancestry inference.

This module provides functionality to generate comprehensive HTML reports
summarizing ancestry inference results, including visualizations, statistics,
and sample predictions.
"""

import logging
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, Optional, Union

import pandas as pd

from hvantk.ancestry.constants import (
    PREDICTED_ANCESTRY_COL,
    ANCESTRY_PROB_COL,
    SOURCE_COL,
    KNOWN_ANCESTRY_COL,
    POPULATION_NAMES,
)

if TYPE_CHECKING:
    from hvantk.ancestry.pipeline import AncestryInferenceResult

logger = logging.getLogger(__name__)


# HTML Report Template
REPORT_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            line-height: 1.6;
            color: #333;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background-color: #f8f9fa;
        }}
        h1 {{
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #34495e;
            border-bottom: 2px solid #bdc3c7;
            padding-bottom: 8px;
            margin-top: 30px;
        }}
        h3 {{
            color: #7f8c8d;
        }}
        .metadata {{
            background-color: #ecf0f1;
            padding: 15px;
            border-radius: 8px;
            margin-bottom: 20px;
        }}
        .metadata p {{
            margin: 5px 0;
        }}
        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin: 20px 0;
        }}
        .summary-card {{
            background: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            text-align: center;
        }}
        .summary-card .value {{
            font-size: 2em;
            font-weight: bold;
            color: #3498db;
        }}
        .summary-card .label {{
            color: #7f8c8d;
            font-size: 0.9em;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            background: white;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        th, td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        th {{
            background-color: #3498db;
            color: white;
        }}
        tr:hover {{
            background-color: #f5f5f5;
        }}
        .plot-container {{
            background: white;
            padding: 20px;
            border-radius: 8px;
            margin: 20px 0;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            text-align: center;
        }}
        .plot-container img {{
            max-width: 100%;
            height: auto;
        }}
        .config-table {{
            font-size: 0.9em;
        }}
        .config-table td:first-child {{
            font-weight: bold;
            width: 200px;
        }}
        .warning {{
            background-color: #fff3cd;
            border: 1px solid #ffc107;
            padding: 15px;
            border-radius: 8px;
            margin: 10px 0;
        }}
        .success {{
            background-color: #d4edda;
            border: 1px solid #28a745;
            padding: 15px;
            border-radius: 8px;
            margin: 10px 0;
        }}
        .footer {{
            margin-top: 40px;
            padding-top: 20px;
            border-top: 1px solid #ddd;
            text-align: center;
            color: #7f8c8d;
            font-size: 0.9em;
        }}
        .ancestry-badge {{
            display: inline-block;
            padding: 4px 12px;
            border-radius: 20px;
            color: white;
            font-weight: bold;
            font-size: 0.85em;
        }}
        .ancestry-AFR {{ background-color: #ff7f0e; }}
        .ancestry-AMR {{ background-color: #9467bd; }}
        .ancestry-EAS {{ background-color: #2ca02c; }}
        .ancestry-EUR {{ background-color: #1f77b4; }}
        .ancestry-SAS {{ background-color: #d62728; }}
        .ancestry-unassigned {{ background-color: #7f7f7f; }}
    </style>
</head>
<body>
    <h1>{title}</h1>

    <div class="metadata">
        <p><strong>Generated:</strong> {timestamp}</p>
        <p><strong>Pipeline Version:</strong> hvantk ancestry-inference</p>
    </div>

    <h2>Summary</h2>
    <div class="summary-grid">
        {summary_cards}
    </div>

    <h2>Ancestry Distribution</h2>
    {ancestry_section}

    <h2>PCA Visualization</h2>
    {pca_section}

    <h2>Model Performance</h2>
    {validation_section}

    <h2>Probability Distribution</h2>
    {probability_section}

    <h2>Sample Predictions</h2>
    {predictions_section}

    <h2>Pipeline Configuration</h2>
    {config_section}

    <div class="footer">
        <p>Generated by hvantk ancestry-inference pipeline</p>
        <p><a href="https://github.com/hvantk/hvantk">hvantk documentation</a></p>
    </div>
</body>
</html>
"""


def _create_summary_card(value: Any, label: str) -> str:
    """Create HTML for a summary card."""
    return f"""
    <div class="summary-card">
        <div class="value">{value}</div>
        <div class="label">{label}</div>
    </div>
    """


def _create_ancestry_table(predictions_df: pd.DataFrame) -> str:
    """Create ancestry distribution table."""
    query_df = predictions_df[predictions_df[SOURCE_COL] == "query"]
    counts = query_df[PREDICTED_ANCESTRY_COL].value_counts()
    total = len(query_df)

    rows = []
    for ancestry, count in counts.items():
        pct = 100 * count / total
        pop_name = POPULATION_NAMES.get(ancestry, ancestry)
        badge_class = f"ancestry-{ancestry}" if ancestry in ["AFR", "AMR", "EAS", "EUR", "SAS", "unassigned"] else ""
        rows.append(f"""
        <tr>
            <td><span class="ancestry-badge {badge_class}">{ancestry}</span></td>
            <td>{pop_name}</td>
            <td>{count:,}</td>
            <td>{pct:.1f}%</td>
        </tr>
        """)

    return f"""
    <table>
        <thead>
            <tr>
                <th>Code</th>
                <th>Population</th>
                <th>Count</th>
                <th>Percentage</th>
            </tr>
        </thead>
        <tbody>
            {''.join(rows)}
        </tbody>
    </table>
    """


def _create_predictions_table(predictions_df: pd.DataFrame, max_rows: int = 100) -> str:
    """Create sample predictions table."""
    query_df = predictions_df[predictions_df[SOURCE_COL] == "query"].copy()

    # Sort by probability descending
    if ANCESTRY_PROB_COL in query_df.columns:
        query_df = query_df.sort_values(ANCESTRY_PROB_COL, ascending=False)

    # Limit rows
    if len(query_df) > max_rows:
        query_df = query_df.head(max_rows)
        note = f"<p><em>Showing top {max_rows} samples by probability. Full results available in TSV export.</em></p>"
    else:
        note = ""

    # Select columns to display
    display_cols = ["s", PREDICTED_ANCESTRY_COL, ANCESTRY_PROB_COL]
    # Add PC columns if present
    pc_cols = [c for c in query_df.columns if c.startswith("PC") and c[2:].isdigit()]
    display_cols.extend(pc_cols[:3])  # Show first 3 PCs

    display_df = query_df[[c for c in display_cols if c in query_df.columns]].copy()

    # Format probability
    if ANCESTRY_PROB_COL in display_df.columns:
        # Convert to numeric if needed (e.g., when loaded from TSV)
        if not pd.api.types.is_numeric_dtype(display_df[ANCESTRY_PROB_COL]):
            display_df[ANCESTRY_PROB_COL] = pd.to_numeric(
                display_df[ANCESTRY_PROB_COL], errors="coerce"
            )
        display_df[ANCESTRY_PROB_COL] = display_df[ANCESTRY_PROB_COL].apply(
            lambda x: f"{x:.3f}" if pd.notna(x) else ""
        )

    # Format PC values
    for col in pc_cols[:3]:
        if col in display_df.columns:
            # Convert to numeric if needed (e.g., when loaded from TSV)
            if not pd.api.types.is_numeric_dtype(display_df[col]):
                display_df[col] = pd.to_numeric(display_df[col], errors="coerce")
            display_df[col] = display_df[col].apply(
                lambda x: f"{x:.4f}" if pd.notna(x) else ""
            )

    # Build table
    headers = "<tr>" + "".join(f"<th>{c}</th>" for c in display_df.columns) + "</tr>"
    rows = []
    for _, row in display_df.iterrows():
        cells = "".join(f"<td>{row[c]}</td>" for c in display_df.columns)
        rows.append(f"<tr>{cells}</tr>")

    return f"""
    {note}
    <table>
        <thead>{headers}</thead>
        <tbody>{''.join(rows)}</tbody>
    </table>
    """


def _create_config_table(config_dict: Dict[str, Any]) -> str:
    """Create configuration table."""
    rows = []
    for key, value in config_dict.items():
        rows.append(f"<tr><td>{key}</td><td>{value}</td></tr>")

    return f"""
    <table class="config-table">
        <thead>
            <tr><th>Parameter</th><th>Value</th></tr>
        </thead>
        <tbody>
            {''.join(rows)}
        </tbody>
    </table>
    """


def _create_validation_section(
    result: "AncestryInferenceResult",
    include_confusion_matrix: bool = True,
) -> str:
    """Create validation metrics section."""
    from hvantk.ancestry.plot import (
        plot_confusion_matrix,
        encode_figure_to_base64,
        close_figure,
    )

    accuracy = result.get_accuracy()

    if accuracy is None:
        return """
        <div class="warning">
            <p>Model validation was skipped. Enable with <code>--validate-model</code> to see performance metrics.</p>
        </div>
        """

    content = f"""
    <div class="success">
        <p><strong>Cross-validation Accuracy:</strong> {accuracy:.2%}</p>
    </div>
    """

    # Add confusion matrix if available
    if include_confusion_matrix and result.classification_result.confusion_matrix is not None:
        try:
            y_true = result.classification_result.confusion_matrix_labels[0]
            y_pred = result.classification_result.confusion_matrix_labels[1]
            labels = result.classification_result.classes

            fig = plot_confusion_matrix(y_true, y_pred, labels=labels, normalize=True)
            img_str = encode_figure_to_base64(fig)
            close_figure(fig)

            content += f"""
            <div class="plot-container">
                <h3>Confusion Matrix (Cross-validation)</h3>
                <img src="{img_str}" alt="Confusion Matrix">
            </div>
            """
        except Exception as e:
            logger.warning(f"Could not generate confusion matrix: {e}")

    return content


def generate_ancestry_report(
    result: "AncestryInferenceResult",
    output_path: Union[str, Path],
    title: str = "Ancestry Inference Report",
    max_table_rows: int = 100,
    include_pca_panel: bool = True,
) -> Path:
    """Generate comprehensive HTML report for ancestry inference results.

    Parameters
    ----------
    result : AncestryInferenceResult
        Complete ancestry inference results.
    output_path : str or Path
        Output file path for HTML report.
    title : str
        Report title. Default: "Ancestry Inference Report".
    max_table_rows : int
        Maximum rows to include in sample predictions table.
    include_pca_panel : bool
        Whether to include two-panel PCA plot. Default: True.

    Returns
    -------
    Path
        Path to generated report.
    """
    from hvantk.ancestry.plot import (
        plot_pca_scatter,
        plot_pca_panel,
        plot_ancestry_proportions,
        plot_probability_distribution,
        plot_variance_explained,
        encode_figure_to_base64,
        close_figure,
    )

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    predictions_df = result.get_predictions_df()
    stats = result.pipeline_stats
    config = result.config

    # Generate summary cards
    query_df = predictions_df[predictions_df[SOURCE_COL] == "query"]
    n_assigned = (query_df[PREDICTED_ANCESTRY_COL] != "unassigned").sum()
    n_unassigned = (query_df[PREDICTED_ANCESTRY_COL] == "unassigned").sum()

    summary_cards = "".join([
        _create_summary_card(f"{stats.get('n_query_samples', len(query_df)):,}", "Query Samples"),
        _create_summary_card(f"{stats.get('n_reference_samples', 0):,}", "Reference Samples"),
        _create_summary_card(f"{n_assigned:,}", "Assigned"),
        _create_summary_card(f"{n_unassigned:,}", "Unassigned"),
        _create_summary_card(f"{stats.get('n_populations', 0)}", "Populations"),
        _create_summary_card(f"{stats.get('n_shared_variants', 0):,}", "Shared Variants"),
    ])

    # Generate ancestry section
    ancestry_section = _create_ancestry_table(predictions_df)
    try:
        fig = plot_ancestry_proportions(predictions_df)
        img_str = encode_figure_to_base64(fig)
        close_figure(fig)
        ancestry_section += f"""
        <div class="plot-container">
            <img src="{img_str}" alt="Ancestry Distribution">
        </div>
        """
    except Exception as e:
        logger.warning(f"Could not generate ancestry plot: {e}")

    # Generate PCA section
    pca_section = ""
    try:
        if include_pca_panel:
            fig = plot_pca_panel(predictions_df, show_query_as_undefined=False)
        else:
            fig = plot_pca_scatter(predictions_df, show_query_as_undefined=False)
        img_str = encode_figure_to_base64(fig, dpi=150)
        close_figure(fig)
        pca_section = f"""
        <div class="plot-container">
            <img src="{img_str}" alt="PCA Plot">
        </div>
        """

        # Add variance explained
        if result.eigenvalues:
            fig = plot_variance_explained(result.eigenvalues, n_pcs=10)
            img_str = encode_figure_to_base64(fig)
            close_figure(fig)
            pca_section += f"""
            <div class="plot-container">
                <h3>Variance Explained by Principal Components</h3>
                <img src="{img_str}" alt="Variance Explained">
            </div>
            """
    except Exception as e:
        logger.warning(f"Could not generate PCA plots: {e}")
        pca_section = f"<p class='warning'>Could not generate PCA plots: {e}</p>"

    # Generate validation section
    validation_section = _create_validation_section(result)

    # Generate probability distribution section
    probability_section = ""
    try:
        if ANCESTRY_PROB_COL in predictions_df.columns:
            fig = plot_probability_distribution(predictions_df)
            img_str = encode_figure_to_base64(fig)
            close_figure(fig)
            probability_section = f"""
            <div class="plot-container">
                <img src="{img_str}" alt="Probability Distribution">
            </div>
            """
    except Exception as e:
        logger.warning(f"Could not generate probability plot: {e}")

    # Generate predictions table
    predictions_section = _create_predictions_table(predictions_df, max_rows=max_table_rows)

    # Generate config section
    config_section = _create_config_table(config.to_dict())

    # Render template
    html_content = REPORT_TEMPLATE.format(
        title=title,
        timestamp=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        summary_cards=summary_cards,
        ancestry_section=ancestry_section,
        pca_section=pca_section,
        validation_section=validation_section,
        probability_section=probability_section,
        predictions_section=predictions_section,
        config_section=config_section,
    )

    # Write report
    output_path.write_text(html_content)
    logger.info(f"Generated ancestry report: {output_path}")

    return output_path
