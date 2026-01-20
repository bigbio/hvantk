"""
Visualization module for hvantk.

This module provides functions and classes for visualizing multiomics data.
"""

from .base import (
    set_default_style,
    save_figure,
    get_colors,
    add_figure_labels,
)

# Lazy facade to avoid importing heavy backends at module import time


def visualize_expression_distribution(*args, **kwargs):
    from .expression.hail import visualize_expression_distribution as _impl

    return _impl(*args, **kwargs)


# QC plotting functions (lazy imports)
def plot_sample_call_rate_distribution(*args, **kwargs):
    from .qc_plots import plot_sample_call_rate_distribution as _impl

    return _impl(*args, **kwargs)


def plot_sample_titv_distribution(*args, **kwargs):
    from .qc_plots import plot_sample_titv_distribution as _impl

    return _impl(*args, **kwargs)


def plot_sample_depth_distribution(*args, **kwargs):
    from .qc_plots import plot_sample_depth_distribution as _impl

    return _impl(*args, **kwargs)


def plot_sample_qc_overview(*args, **kwargs):
    from .qc_plots import plot_sample_qc_overview as _impl

    return _impl(*args, **kwargs)


def plot_variant_call_rate_distribution(*args, **kwargs):
    from .qc_plots import plot_variant_call_rate_distribution as _impl

    return _impl(*args, **kwargs)


def plot_allele_frequency_spectrum(*args, **kwargs):
    from .qc_plots import plot_allele_frequency_spectrum as _impl

    return _impl(*args, **kwargs)


def plot_hwe_pvalues(*args, **kwargs):
    from .qc_plots import plot_hwe_pvalues as _impl

    return _impl(*args, **kwargs)


def plot_variant_qc_overview(*args, **kwargs):
    from .qc_plots import plot_variant_qc_overview as _impl

    return _impl(*args, **kwargs)


def plot_qc_summary_dashboard(*args, **kwargs):
    from .qc_plots import plot_qc_summary_dashboard as _impl

    return _impl(*args, **kwargs)


def generate_qc_report(*args, **kwargs):
    from .qc_report import generate_qc_report as _impl

    return _impl(*args, **kwargs)


# Interactive QC plotting functions (lazy imports)
def plot_interactive_sample_call_rates(*args, **kwargs):
    from .interactive_qc import plot_interactive_sample_call_rates as _impl

    return _impl(*args, **kwargs)


def plot_interactive_sample_titv(*args, **kwargs):
    from .interactive_qc import plot_interactive_sample_titv as _impl

    return _impl(*args, **kwargs)


def plot_interactive_variant_call_rates(*args, **kwargs):
    from .interactive_qc import plot_interactive_variant_call_rates as _impl

    return _impl(*args, **kwargs)


def plot_interactive_allele_frequencies(*args, **kwargs):
    from .interactive_qc import plot_interactive_allele_frequencies as _impl

    return _impl(*args, **kwargs)


def plot_interactive_hwe_pvalues(*args, **kwargs):
    from .interactive_qc import plot_interactive_hwe_pvalues as _impl

    return _impl(*args, **kwargs)


def plot_interactive_qc_dashboard(*args, **kwargs):
    from .interactive_qc import plot_interactive_qc_dashboard as _impl

    return _impl(*args, **kwargs)


def plot_interactive_sample_scatter(*args, **kwargs):
    from .interactive_qc import plot_interactive_sample_scatter as _impl

    return _impl(*args, **kwargs)


__all__ = [
    "set_default_style",
    "save_figure",
    "get_colors",
    "add_figure_labels",
    "visualize_expression_distribution",
    # QC plotting functions
    "plot_sample_call_rate_distribution",
    "plot_sample_titv_distribution",
    "plot_sample_depth_distribution",
    "plot_sample_qc_overview",
    "plot_variant_call_rate_distribution",
    "plot_allele_frequency_spectrum",
    "plot_hwe_pvalues",
    "plot_variant_qc_overview",
    "plot_qc_summary_dashboard",
    "generate_qc_report",
    # Interactive plotting functions
    "plot_interactive_sample_call_rates",
    "plot_interactive_sample_titv",
    "plot_interactive_variant_call_rates",
    "plot_interactive_allele_frequencies",
    "plot_interactive_hwe_pvalues",
    "plot_interactive_qc_dashboard",
    "plot_interactive_sample_scatter",
]
