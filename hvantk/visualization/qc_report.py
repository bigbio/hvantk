"""
QC HTML Report Generator

This module creates comprehensive HTML reports with embedded plots and summary statistics
for quality control analysis of genomic variant data.

The module provides:
- Self-contained HTML reports with embedded plots
- Professional styling and responsive design
- Quality assessment with color-coded status indicators
- Automated recommendations based on QC thresholds
- Summary statistics and metrics tables
- Support for custom thresholds and parameters

Example:
    >>> from hvantk.hgc import compute_full_qc
    >>> from hvantk.visualization.qc_report import generate_qc_report
    >>> qc_results = compute_full_qc(mt)
    >>> report_path = generate_qc_report(qc_results, 'qc_report.html')
"""

import logging
import base64
from io import BytesIO
from pathlib import Path
from datetime import datetime
from typing import List, Union, Optional

import pandas as pd
import matplotlib

matplotlib.use(
    "Agg"
)  # Use non-interactive backend - must be set before importing pyplot
import matplotlib.pyplot as plt

from .qc_plots import (
    plot_sample_qc_overview,
    plot_variant_qc_overview,
    plot_sample_call_rate_distribution,
    plot_variant_call_rate_distribution,
    plot_allele_frequency_spectrum,
    plot_hwe_pvalues,
    plot_sample_titv_distribution,
    _prepare_sample_qc_data,
    _prepare_variant_qc_data,
)

logger = logging.getLogger(__name__)

# HTML template for the QC report
HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
            line-height: 1.6;
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 30px;
            text-align: center;
        }}
        .header h1 {{
            margin: 0;
            font-size: 2.5em;
            font-weight: 300;
        }}
        .header .subtitle {{
            margin: 10px 0 0 0;
            font-size: 1.1em;
            opacity: 0.9;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
        }}
        .section {{
            background: white;
            padding: 30px;
            margin-bottom: 30px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .section h2 {{
            color: #333;
            margin-top: 0;
            border-bottom: 3px solid #667eea;
            padding-bottom: 10px;
        }}
        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        .summary-card {{
            background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            text-align: center;
        }}
        .summary-card.samples {{
            background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%);
        }}
        .summary-card.variants {{
            background: linear-gradient(135deg, #43e97b 0%, #38f9d7 100%);
        }}
        .summary-card.quality {{
            background: linear-gradient(135deg, #fa709a 0%, #fee140 100%);
        }}
        .summary-card h3 {{
            margin: 0 0 10px 0;
            font-size: 1.1em;
        }}
        .summary-card .value {{
            font-size: 2em;
            font-weight: bold;
            margin: 10px 0;
        }}
        .metrics-table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }}
        .metrics-table th, .metrics-table td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        .metrics-table th {{
            background-color: #f8f9fa;
            font-weight: 600;
            color: #555;
        }}
        .metrics-table tr:hover {{
            background-color: #f5f5f5;
        }}
        .plot-container {{
            margin: 30px 0;
            text-align: center;
        }}
        .plot-container img {{
            max-width: 100%;
            height: auto;
            border-radius: 8px;
            box-shadow: 0 4px 8px rgba(0,0,0,0.1);
        }}
        .plot-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(500px, 1fr));
            gap: 30px;
            margin: 30px 0;
        }}
        .status-badge {{
            padding: 4px 12px;
            border-radius: 20px;
            font-size: 0.9em;
            font-weight: bold;
        }}
        .status-pass {{
            background-color: #d4edda;
            color: #155724;
        }}
        .status-warn {{
            background-color: #fff3cd;
            color: #856404;
        }}
        .status-fail {{
            background-color: #f8d7da;
            color: #721c24;
        }}
        .footer {{
            text-align: center;
            padding: 20px;
            color: #666;
            border-top: 1px solid #ddd;
            margin-top: 40px;
        }}
        .toc {{
            background: #f8f9fa;
            padding: 20px;
            border-radius: 8px;
            margin-bottom: 30px;
        }}
        .toc h3 {{
            margin-top: 0;
        }}
        .toc ul {{
            list-style-type: none;
            padding: 0;
        }}
        .toc li {{
            padding: 5px 0;
        }}
        .toc a {{
            color: #667eea;
            text-decoration: none;
        }}
        .toc a:hover {{
            text-decoration: underline;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>{title}</h1>
            <div class="subtitle">Generated on {timestamp}</div>
        </div>

        <div class="toc">
            <h3>📋 Table of Contents</h3>
            <ul>
                <li><a href="#overview">📊 Overview</a></li>
                <li><a href="#sample-qc">👥 Sample Quality Control</a></li>
                <li><a href="#variant-qc">🧬 Variant Quality Control</a></li>
                <li><a href="#summary">📈 Quality Assessment Summary</a></li>
            </ul>
        </div>

        <div class="section" id="overview">
            <h2>📊 Overview</h2>
            <div class="summary-grid">
                <div class="summary-card samples">
                    <h3>Samples</h3>
                    <div class="value">{n_samples:,}</div>
                    <div>Total samples analyzed</div>
                </div>
                <div class="summary-card variants">
                    <h3>Variants</h3>
                    <div class="value">{n_variants:,}</div>
                    <div>Total variants analyzed</div>
                </div>
                <div class="summary-card quality">
                    <h3>Overall Quality</h3>
                    <div class="value">{overall_quality}</div>
                    <div>Quality assessment</div>
                </div>
            </div>
            
            <h3>📋 Dataset Summary</h3>
            <table class="metrics-table">
                <tr>
                    <th>Metric</th>
                    <th>Value</th>
                    <th>Status</th>
                </tr>
                {overview_metrics}
            </table>
        </div>

        <div class="section" id="sample-qc">
            <h2>👥 Sample Quality Control</h2>
            
            <h3>📈 Sample QC Overview</h3>
            <div class="plot-container">
                <img src="data:image/png;base64,{sample_overview_plot}" alt="Sample QC Overview">
            </div>
            
            <div class="plot-grid">
                <div class="plot-container">
                    <h4>Sample Call Rate Distribution</h4>
                    <img src="data:image/png;base64,{sample_call_rate_plot}" alt="Sample Call Rates">
                </div>
                <div class="plot-container">
                    <h4>Ti/Tv Ratio Distribution</h4>
                    <img src="data:image/png;base64,{sample_titv_plot}" alt="Sample Ti/Tv Ratios">
                </div>
            </div>
            
            <h3>📊 Sample QC Metrics Summary</h3>
            <table class="metrics-table">
                <tr>
                    <th>Metric</th>
                    <th>Mean ± Std</th>
                    <th>Min</th>
                    <th>Max</th>
                    <th>Status</th>
                </tr>
                {sample_metrics_table}
            </table>
        </div>

        <div class="section" id="variant-qc">
            <h2>🧬 Variant Quality Control</h2>
            
            <h3>📈 Variant QC Overview</h3>
            <div class="plot-container">
                <img src="data:image/png;base64,{variant_overview_plot}" alt="Variant QC Overview">
            </div>
            
            <div class="plot-grid">
                <div class="plot-container">
                    <h4>Variant Call Rate Distribution</h4>
                    <img src="data:image/png;base64,{variant_call_rate_plot}" alt="Variant Call Rates">
                </div>
                <div class="plot-container">
                    <h4>Allele Frequency Spectrum</h4>
                    <img src="data:image/png;base64,{allele_freq_plot}" alt="Allele Frequency Spectrum">
                </div>
            </div>
            
            <div class="plot-grid">
                <div class="plot-container">
                    <h4>Hardy-Weinberg Equilibrium p-values</h4>
                    <img src="data:image/png;base64,{hwe_plot}" alt="HWE p-values">
                </div>
            </div>
            
            <h3>📊 Variant QC Metrics Summary</h3>
            <table class="metrics-table">
                <tr>
                    <th>Metric</th>
                    <th>Mean ± Std</th>
                    <th>Min</th>
                    <th>Max</th>
                    <th>Status</th>
                </tr>
                {variant_metrics_table}
            </table>
        </div>

        <div class="section" id="summary">
            <h2>📈 Quality Assessment Summary</h2>
            
            <h3>🎯 Recommendations</h3>
            <div class="recommendations">
                {recommendations}
            </div>
            
            <h3>⚙️ Analysis Parameters</h3>
            <table class="metrics-table">
                <tr>
                    <th>Parameter</th>
                    <th>Value</th>
                </tr>
                {analysis_parameters}
            </table>
        </div>

        <div class="footer">
            <p>Generated by hvantk QC Report Generator | {timestamp}</p>
            <p>For more information, visit the hvantk documentation</p>
        </div>
    </div>
</body>
</html>
"""


def fig_to_base64(fig) -> str:
    """Convert matplotlib figure to base64 string for HTML embedding."""
    buffer = BytesIO()
    fig.savefig(
        buffer,
        format="png",
        dpi=150,
        bbox_inches="tight",
        facecolor="white",
        edgecolor="none",
    )
    buffer.seek(0)
    plot_data = buffer.getvalue()
    buffer.close()

    return base64.b64encode(plot_data).decode("utf-8")


def create_placeholder_plot(title: str, message: str, figsize: tuple = (10, 6)) -> str:
    """Create a placeholder plot with a message and return as base64 string."""
    fig, ax = plt.subplots(figsize=figsize)
    ax.text(
        0.5,
        0.5,
        message,
        ha="center",
        va="center",
        transform=ax.transAxes,
        fontsize=14,
        color="gray",
    )
    ax.set_title(title)
    ax.axis("off")  # Remove axes for cleaner look

    plot_base64 = fig_to_base64(fig)
    plt.close(fig)
    return plot_base64


def get_status_badge(
    value: float,
    good_threshold: float,
    acceptable_threshold: float,
    higher_is_better: bool = True,
) -> str:
    """Generate HTML status badge based on value and thresholds."""
    if higher_is_better:
        if value >= good_threshold:
            return '<span class="status-badge status-pass">✓ Good</span>'
        elif value >= acceptable_threshold:
            return '<span class="status-badge status-warn">⚠ Acceptable</span>'
        else:
            return '<span class="status-badge status-fail">✗ Poor</span>'
    else:
        if value <= good_threshold:
            return '<span class="status-badge status-pass">✓ Good</span>'
        elif value <= acceptable_threshold:
            return '<span class="status-badge status-warn">⚠ Acceptable</span>'
        else:
            return '<span class="status-badge status-fail">✗ Poor</span>'


def generate_sample_metrics_table(sample_df: pd.DataFrame) -> str:
    """Generate HTML table for sample QC metrics."""
    df = _prepare_sample_qc_data(sample_df)

    metrics_to_show = [
        ("call_rate", "Call Rate", 0.95, 0.85, True),
        ("r_ti_tv", "Ti/Tv Ratio", None, None, None),
        ("dp_stats_mean", "Mean Depth", None, None, None),
        ("n_het", "Heterozygous Calls", None, None, None),
        ("n_singleton", "Singleton Variants", None, None, None),
    ]

    rows = []
    for col, label, good_thresh, accept_thresh, higher_better in metrics_to_show:
        if col in df.columns:
            values = df[col].dropna()
            if len(values) > 0:
                mean_val = values.mean()
                std_val = values.std()
                min_val = values.min()
                max_val = values.max()

                if good_thresh is not None:
                    status = get_status_badge(
                        mean_val, good_thresh, accept_thresh, higher_better
                    )
                else:
                    status = '<span class="status-badge status-pass">-</span>'

                if col == "r_ti_tv":
                    mean_str = f"{mean_val:.3f} ± {std_val:.3f}"
                    min_str = f"{min_val:.3f}"
                    max_str = f"{max_val:.3f}"
                elif col in ["call_rate"]:
                    mean_str = f"{mean_val:.3f} ± {std_val:.3f}"
                    min_str = f"{min_val:.3f}"
                    max_str = f"{max_val:.3f}"
                else:
                    mean_str = f"{mean_val:.1f} ± {std_val:.1f}"
                    min_str = f"{min_val:.1f}"
                    max_str = f"{max_val:.1f}"

                rows.append(
                    f"""
                    <tr>
                        <td>{label}</td>
                        <td>{mean_str}</td>
                        <td>{min_str}</td>
                        <td>{max_str}</td>
                        <td>{status}</td>
                    </tr>
                """
                )

    return "".join(rows)


def generate_variant_metrics_table(variant_df: pd.DataFrame) -> str:
    """Generate HTML table for variant QC metrics."""
    df = _prepare_variant_qc_data(variant_df)

    metrics_to_show = [
        ("call_rate", "Call Rate", 0.90, 0.80, True),
        ("AF_alt", "Mean Allele Frequency", None, None, None),
        ("p_value_hwe", "Mean HWE p-value", 0.01, 0.001, True),
        ("n_het", "Mean Heterozygous Calls", None, None, None),
        ("AC_alt", "Mean Allele Count", None, None, None),
    ]

    rows = []
    for col, label, good_thresh, accept_thresh, higher_better in metrics_to_show:
        if col in df.columns:
            values = df[col].dropna()
            if len(values) > 0:
                mean_val = values.mean()
                std_val = values.std()
                min_val = values.min()
                max_val = values.max()

                if good_thresh is not None:
                    status = get_status_badge(
                        mean_val, good_thresh, accept_thresh, higher_better
                    )
                else:
                    status = '<span class="status-badge status-pass">-</span>'

                if col in ["call_rate", "AF_alt"]:
                    mean_str = f"{mean_val:.4f} ± {std_val:.4f}"
                    min_str = f"{min_val:.4f}"
                    max_str = f"{max_val:.4f}"
                elif col == "p_value_hwe":
                    mean_str = f"{mean_val:.2e} ± {std_val:.2e}"
                    min_str = f"{min_val:.2e}"
                    max_str = f"{max_val:.2e}"
                else:
                    mean_str = f"{mean_val:.1f} ± {std_val:.1f}"
                    min_str = f"{min_val:.1f}"
                    max_str = f"{max_val:.1f}"

                rows.append(
                    f"""
                    <tr>
                        <td>{label}</td>
                        <td>{mean_str}</td>
                        <td>{min_str}</td>
                        <td>{max_str}</td>
                        <td>{status}</td>
                    </tr>
                """
                )

    return "".join(rows)


def generate_recommendations(
    sample_df: Optional[pd.DataFrame], variant_df: Optional[pd.DataFrame]
) -> str:
    """Generate QC recommendations based on the data."""
    recommendations = []

    # Sample recommendations - only if sample QC data is available
    if sample_df is not None:
        sample_data = _prepare_sample_qc_data(sample_df)
        if "call_rate" in sample_data.columns:
            sample_call_rates = sample_data["call_rate"].dropna()
            low_call_rate_samples = (sample_call_rates < 0.85).sum()
            if low_call_rate_samples > 0:
                pct = low_call_rate_samples / len(sample_call_rates) * 100
                recommendations.append(
                    f"""
                    <div class="alert alert-warning">
                        <strong>⚠️ Sample Call Rates:</strong> {low_call_rate_samples} samples ({pct:.1f}%) 
                        have call rates below 85%. Consider removing these samples from analysis.
                    </div>
                """
                )
            else:
                recommendations.append(
                    """
                    <div class="alert alert-success">
                        <strong>✅ Sample Call Rates:</strong> All samples have acceptable call rates (≥85%).
                    </div>
                """
                )
    else:
        recommendations.append(
            """
            <div class="alert alert-info">
                <strong>ℹ️ Sample QC:</strong> Sample QC metrics not available - consider running sample QC analysis.
            </div>
        """
        )

    # Variant recommendations - only if variant QC data is available
    if variant_df is not None:
        variant_data = _prepare_variant_qc_data(variant_df)
        if "call_rate" in variant_data.columns:
            variant_call_rates = variant_data["call_rate"].dropna()
            low_call_rate_variants = (variant_call_rates < 0.80).sum()
            if low_call_rate_variants > 0:
                pct = low_call_rate_variants / len(variant_call_rates) * 100
                recommendations.append(
                    f"""
                    <div class="alert alert-warning">
                        <strong>⚠️ Variant Call Rates:</strong> {low_call_rate_variants:,} variants ({pct:.1f}%) 
                        have call rates below 80%. Consider filtering these variants.
                    </div>
                """
                )
            else:
                recommendations.append(
                    """
                    <div class="alert alert-success">
                        <strong>✅ Variant Call Rates:</strong> Most variants have acceptable call rates (≥80%).
                    </div>
                """
                )

        # HWE recommendations - only if variant QC data is available
        if "p_value_hwe" in variant_data.columns:
            hwe_pvals = variant_data["p_value_hwe"].dropna()
            hwe_failing = (hwe_pvals < 1e-6).sum()
            if hwe_failing > 0:
                pct = hwe_failing / len(hwe_pvals) * 100
                recommendations.append(
                    f"""
                    <div class="alert alert-info">
                        <strong>ℹ️ Hardy-Weinberg Equilibrium:</strong> {hwe_failing:,} variants ({pct:.1f}%) 
                        fail HWE test (p < 1e-6). Review these for potential genotyping errors.
                    </div>
                """
                )
    else:
        recommendations.append(
            """
            <div class="alert alert-info">
                <strong>ℹ️ Variant QC:</strong> Variant QC metrics not available - consider running variant QC analysis.
            </div>
        """
        )

    if not recommendations:
        recommendations.append(
            """
            <div class="alert alert-success">
                <strong>✅ Overall Quality:</strong> No major quality issues detected. 
                Dataset appears to be of good quality for downstream analysis.
            </div>
        """
        )

    return "".join(recommendations)


def generate_qc_report(
    qc_results,
    output_path: Union[str, Path],
    title: str = "Quality Control Report",
    include_plots: Optional[List[str]] = None,
) -> Path:
    """
    Generate comprehensive HTML QC report.

    Parameters
    ----------
    qc_results : QCMetrics
        QC results object with computed metrics
    output_path : str or Path
        Output file path for the HTML report
    title : str
        Report title
    include_plots : list, optional
        List of plots to include. Options: 'sample_overview', 'variant_overview',
        'sample_call_rates', 'variant_call_rates', 'allele_frequencies', 'hwe', 'titv'

    Returns
    -------
    Path
        Path to the generated HTML report

    Example
    -------
    >>> qc_results = compute_full_qc(mt)
    >>> report_path = generate_qc_report(qc_results, 'qc_report.html')
    """
    if include_plots is None:
        include_plots = [
            "sample_overview",
            "variant_overview",
            "sample_call_rates",
            "variant_call_rates",
            "allele_frequencies",
            "hwe",
            "titv",
        ]

    logger.info("Generating comprehensive QC HTML report...")

    # Get data - check if QC is available before calling getters
    sample_df = qc_results.get_sample_metrics_df() if qc_results.has_sample_qc else None
    variant_df = (
        qc_results.get_variant_metrics_df() if qc_results.has_variant_qc else None
    )

    # Basic stats with fallback values
    n_samples = len(sample_df) if sample_df is not None else 0
    n_variants = len(variant_df) if variant_df is not None else 0
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Generate plots and convert to base64
    plot_data = {}

    # Initialize all plot keys with placeholder images to avoid KeyError in template formatting
    # These will be overwritten with actual plots if they're in include_plots
    plot_placeholders = {
        "sample_overview_plot": (
            "Sample QC Overview",
            "Plot not included\nin this report",
            (15, 10),
        ),
        "sample_call_rate_plot": (
            "Sample Call Rate Distribution",
            "Plot not included\nin this report",
            (10, 6),
        ),
        "sample_titv_plot": (
            "Sample Ti/Tv Ratio Distribution",
            "Plot not included\nin this report",
            (10, 6),
        ),
        "variant_overview_plot": (
            "Variant QC Overview",
            "Plot not included\nin this report",
            (15, 10),
        ),
        "variant_call_rate_plot": (
            "Variant Call Rate Distribution",
            "Plot not included\nin this report",
            (10, 6),
        ),
        "allele_freq_plot": (
            "Allele Frequency Spectrum",
            "Plot not included\nin this report",
            (10, 6),
        ),
        "hwe_plot": (
            "Hardy-Weinberg Equilibrium P-values",
            "Plot not included\nin this report",
            (10, 6),
        ),
    }

    # Initialize with placeholders
    for key, (plot_title, message, figsize) in plot_placeholders.items():
        plot_data[key] = create_placeholder_plot(plot_title, message, figsize)

    # Sample plots - only generate if sample QC is available
    if qc_results.has_sample_qc and sample_df is not None:
        if "sample_overview" in include_plots:
            logger.info("Generating sample overview plot...")
            fig = plot_sample_qc_overview(sample_df, figsize=(15, 10))
            plot_data["sample_overview_plot"] = fig_to_base64(fig)
            plt.close(fig)

        if "sample_call_rates" in include_plots:
            logger.info("Generating sample call rate plot...")
            fig = plot_sample_call_rate_distribution(sample_df, figsize=(10, 6))
            plot_data["sample_call_rate_plot"] = fig_to_base64(fig)
            plt.close(fig)

        if "titv" in include_plots:
            logger.info("Generating sample Ti/Tv plot...")
            try:
                fig = plot_sample_titv_distribution(sample_df, figsize=(10, 6))
                plot_data["sample_titv_plot"] = fig_to_base64(fig)
                plt.close(fig)
            except (ValueError, KeyError, AttributeError, TypeError) as e:
                # Handle expected recoverable errors (missing data, wrong types, etc.)
                logger.warning(f"Ti/Tv plot failed with {type(e).__name__}: {e}")
                # Create placeholder plot
                plot_data["sample_titv_plot"] = create_placeholder_plot(
                    "Sample Ti/Tv Ratio Distribution",
                    "Ti/Tv data not available\nfor this dataset",
                    (10, 6),
                )
            except Exception:
                # Unexpected errors should not be masked - re-raise them
                logger.exception("Unexpected error generating Ti/Tv plot")
                raise
    else:
        # Update placeholders with "QC not available" message for requested plots when no sample QC
        if "sample_overview" in include_plots:
            plot_data["sample_overview_plot"] = create_placeholder_plot(
                "Sample QC Overview",
                "Sample QC not available\nfor this dataset",
                (15, 10),
            )

        if "sample_call_rates" in include_plots:
            plot_data["sample_call_rate_plot"] = create_placeholder_plot(
                "Sample Call Rate Distribution",
                "Sample call rate data\nnot available",
                (10, 6),
            )

        if "titv" in include_plots:
            plot_data["sample_titv_plot"] = create_placeholder_plot(
                "Sample Ti/Tv Ratio Distribution",
                "Sample Ti/Tv data\nnot available",
                (10, 6),
            )

    # Variant plots - only generate if variant QC is available
    if qc_results.has_variant_qc and variant_df is not None:
        if "variant_overview" in include_plots:
            logger.info("Generating variant overview plot...")
            fig = plot_variant_qc_overview(variant_df, figsize=(15, 10))
            plot_data["variant_overview_plot"] = fig_to_base64(fig)
            plt.close(fig)

        if "variant_call_rates" in include_plots:
            logger.info("Generating variant call rate plot...")
            fig = plot_variant_call_rate_distribution(variant_df, figsize=(10, 6))
            plot_data["variant_call_rate_plot"] = fig_to_base64(fig)
            plt.close(fig)

        if "allele_frequencies" in include_plots:
            logger.info("Generating allele frequency plot...")
            fig = plot_allele_frequency_spectrum(variant_df, figsize=(10, 6))
            plot_data["allele_freq_plot"] = fig_to_base64(fig)
            plt.close(fig)

        if "hwe" in include_plots:
            logger.info("Generating HWE plot...")
            fig = plot_hwe_pvalues(variant_df, figsize=(10, 6))
            plot_data["hwe_plot"] = fig_to_base64(fig)
            plt.close(fig)
    else:
        # Update placeholders with "QC not available" message for requested plots when no variant QC
        if "variant_overview" in include_plots:
            plot_data["variant_overview_plot"] = create_placeholder_plot(
                "Variant QC Overview",
                "Variant QC not available\nfor this dataset",
                (15, 10),
            )

        if "variant_call_rates" in include_plots:
            plot_data["variant_call_rate_plot"] = create_placeholder_plot(
                "Variant Call Rate Distribution",
                "Variant call rate data\nnot available",
                (10, 6),
            )

        if "allele_frequencies" in include_plots:
            plot_data["allele_freq_plot"] = create_placeholder_plot(
                "Allele Frequency Spectrum",
                "Allele frequency data\nnot available",
                (10, 6),
            )

        if "hwe" in include_plots:
            plot_data["hwe_plot"] = create_placeholder_plot(
                "Hardy-Weinberg Equilibrium P-values",
                "Hardy-Weinberg data\nnot available",
                (10, 6),
            )

    # Generate metrics tables with fallback for missing QC data
    if qc_results.has_sample_qc and sample_df is not None:
        sample_metrics_table = generate_sample_metrics_table(sample_df)
    else:
        sample_metrics_table = '<tr><td colspan="6" style="text-align: center; color: #666;">Sample QC metrics not available</td></tr>'

    if qc_results.has_variant_qc and variant_df is not None:
        variant_metrics_table = generate_variant_metrics_table(variant_df)
    else:
        variant_metrics_table = '<tr><td colspan="6" style="text-align: center; color: #666;">Variant QC metrics not available</td></tr>'

    # Generate overview metrics - prepare data with conditional checks
    sample_data = _prepare_sample_qc_data(sample_df) if sample_df is not None else None
    variant_data = (
        _prepare_variant_qc_data(variant_df) if variant_df is not None else None
    )

    overview_metrics = []

    # Sample QC metrics - only if sample QC is available
    if sample_data is not None and "call_rate" in sample_data.columns:
        sample_cr = sample_data["call_rate"].mean()
        status = get_status_badge(sample_cr, 0.95, 0.85, True)
        overview_metrics.append(
            f"""
            <tr>
                <td>Mean Sample Call Rate</td>
                <td>{sample_cr:.3f}</td>
                <td>{status}</td>
            </tr>
        """
        )
    else:
        overview_metrics.append(
            f"""
            <tr>
                <td>Mean Sample Call Rate</td>
                <td>N/A</td>
                <td><span class="badge badge-warning">No QC Data</span></td>
            </tr>
        """
        )

    # Variant QC metrics - only if variant QC is available
    if variant_data is not None and "call_rate" in variant_data.columns:
        variant_cr = variant_data["call_rate"].mean()
        status = get_status_badge(variant_cr, 0.90, 0.80, True)
        overview_metrics.append(
            f"""
            <tr>
                <td>Mean Variant Call Rate</td>
                <td>{variant_cr:.3f}</td>
                <td>{status}</td>
            </tr>
        """
        )
    else:
        overview_metrics.append(
            f"""
            <tr>
                <td>Mean Variant Call Rate</td>
                <td>N/A</td>
                <td><span class="badge badge-warning">No QC Data</span></td>
            </tr>
        """
        )

    # Overall quality assessment with fallback logic
    sample_quality = (
        "Good"
        if (
            sample_data is not None
            and "call_rate" in sample_data.columns
            and sample_data["call_rate"].mean() > 0.90
        )
        else "Review"
    )
    variant_quality = (
        "Good"
        if (
            variant_data is not None
            and "call_rate" in variant_data.columns
            and variant_data["call_rate"].mean() > 0.85
        )
        else "Review"
    )

    # If no QC data is available, mark as "Incomplete"
    if sample_data is None and variant_data is None:
        overall_quality = "Incomplete"
    elif sample_quality == "Good" and variant_quality == "Good":
        overall_quality = "Good"
    else:
        overall_quality = "Review"

    # Generate recommendations - handle None dataframes
    recommendations = generate_recommendations(sample_df, variant_df)

    # Analysis parameters
    analysis_parameters = f"""
        <tr><td>Number of Samples</td><td>{n_samples:,}</td></tr>
        <tr><td>Number of Variants</td><td>{n_variants:,}</td></tr>
        <tr><td>Analysis Date</td><td>{timestamp}</td></tr>
        <tr><td>QC Module Version</td><td>hvantk v1.0</td></tr>
    """

    # Format HTML
    html_content = HTML_TEMPLATE.format(
        title=title,
        timestamp=timestamp,
        n_samples=n_samples,
        n_variants=n_variants,
        overall_quality=overall_quality,
        overview_metrics="".join(overview_metrics),
        sample_metrics_table=sample_metrics_table,
        variant_metrics_table=variant_metrics_table,
        recommendations=recommendations,
        analysis_parameters=analysis_parameters,
        **plot_data,
    )

    # Write HTML file
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html_content)

    logger.info(f"QC HTML report generated: {output_path}")
    return output_path
