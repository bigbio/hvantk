"""Rendering helpers for HGC quality-control summary reports.

Pure presentation logic kept out of the CLI layer so it can be reused from the
Python API and unit-tested in isolation. The input is the ``summary_data`` dict
assembled by ``hvantk hgc qc-summary``:

    {
        "<qc_type>": {
            "file": str,
            "n_samples" | "n_variants": int,
            "summary_stats": {metric: {stat: value, ...}, ...},
            "columns": list[str],
        },
        ...
    }
"""

from typing import Dict

# Order of describe()-style statistics rendered as Markdown table columns.
_SUMMARY_STAT_COLUMNS = ["count", "mean", "std", "min", "25%", "50%", "75%", "max"]


def render_qc_summary_markdown(summary_data: Dict[str, dict]) -> str:
    """Render a QC summary dict into a Markdown report string.

    Parameters
    ----------
    summary_data : dict
        Mapping of ``qc_type`` to a per-type summary dict (see module docstring).

    Returns
    -------
    str
        The Markdown report.
    """
    md_content = "# Quality Control Summary Report\n\n"

    for qc_type, data in summary_data.items():
        md_content += f"## {qc_type.replace('_', ' ').title()}\n\n"
        md_content += f"- **File**: {data['file']}\n"
        md_content += (
            f"- **Count**: {data.get('n_samples', data.get('n_variants', 0)):,}\n"
        )
        md_content += f"- **Metrics**: {', '.join(data['columns'])}\n\n"

        if data["summary_stats"]:
            md_content += "### Summary Statistics\n\n"
            md_content += (
                "| Metric | Count | Mean | Std | Min | 25% | 50% | 75% | Max |\n"
            )
            md_content += "|--------|--------|--------|--------|--------|--------|--------|--------|--------|\n"

            for metric, stats in data["summary_stats"].items():
                if isinstance(stats, dict):
                    row = f"| {metric} |"
                    for stat in _SUMMARY_STAT_COLUMNS:
                        value = stats.get(stat, "N/A")
                        if isinstance(value, (int, float)) and stat != "count":
                            # Use scientific notation for very small magnitudes
                            # (e.g. tiny HWE p-values) so they don't collapse to
                            # "0.000"; keep an exact 0 as "0.000".
                            value = (
                                f"{value:.3f}"
                                if value == 0 or 1e-3 <= abs(value) < 1000
                                else f"{value:.2e}"
                            )
                        row += f" {value} |"
                    md_content += row + "\n"
            md_content += "\n"

    return md_content
