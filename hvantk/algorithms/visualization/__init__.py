"""
Visualization module for hvantk.

This module provides functions and classes for visualizing multiomics data.
"""

from .base import (
    set_default_style,
    save_figure,
    save_figure_to_path,
    encode_figure_to_base64,
    get_colors,
    add_figure_labels,
)

# Lazy facade to avoid importing heavy backends at module import time


def visualize_expression_distribution(*args, **kwargs):
    from .expression.anndata import visualize_expression_distribution as _impl

    return _impl(*args, **kwargs)


def generate_qc_report(*args, **kwargs):
    """Generate the static HTML QC triage report.

    For plotting QC metrics yourself (e.g. for a publication figure), get the
    tables via ``QCMetrics.get_sample_metrics_df()`` /
    ``QCMetrics.get_variant_metrics_df()`` and plot with matplotlib directly.
    """
    from .qc_report import generate_qc_report as _impl

    return _impl(*args, **kwargs)


__all__ = [
    "set_default_style",
    "save_figure",
    "save_figure_to_path",
    "encode_figure_to_base64",
    "get_colors",
    "add_figure_labels",
    "visualize_expression_distribution",
    "generate_qc_report",
]
