"""
Expression visualization utilities.

This module contains plotting functions that operate on Hail MatrixTable
objects and AnnData objects.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

try:
    import hail as hl
except ImportError:  # pragma: no cover
    hl = None  # type: ignore[assignment]


def visualize_expression_distribution(
    mt: hl.MatrixTable,
    n_bins: int = 50,
    log_scale: bool = True,
    title: str = "Expression Value Distribution",
    expr_field: str = "x",
) -> plt.Figure:
    """
    Visualize the distribution of expression values in the MatrixTable using Hail's native histogram.

    Args:
        mt: Hail MatrixTable to visualize
        n_bins: Number of histogram bins
        log_scale: Whether to use log scale for values
        title: Plot title
        expr_field: Name of the entry field containing expression values (default: 'x')

    Returns:
        Matplotlib figure
    """
    # Verify the expression field exists
    if expr_field not in mt.entry.dtype.fields:
        raise ValueError(
            f"Expression field '{expr_field}' not found in entry fields. Available fields: {list(mt.entry.dtype.fields)}"
        )

    # Ensure n_bins is an integer
    n_bins = int(n_bins)

    # Compute min/max for histogram range
    if log_scale:
        min_val = mt.aggregate_entries(
            hl.agg.filter(mt[expr_field] > 0, hl.agg.min(mt[expr_field]))
        )
        max_val = mt.aggregate_entries(
            hl.agg.filter(mt[expr_field] > 0, hl.agg.max(mt[expr_field]))
        )
        if min_val is None or max_val is None or min_val <= 0 or max_val <= 0:
            raise ValueError("No positive values found for log transformation")
        log_min = float(np.floor(np.log10(min_val)))
        log_max = float(np.ceil(np.log10(max_val)))
        if log_min == log_max:
            log_min -= 1
            log_max += 1

        # Use plain Python integer for n_bins
        expr = hl.log10(hl.float64(mt[expr_field]))
        hist = mt.aggregate_entries(hl.agg.hist(expr, log_min, log_max, n_bins))
        x_label = "log10(Expression)"
    else:
        min_val = mt.aggregate_entries(hl.agg.min(mt[expr_field]))
        max_val = mt.aggregate_entries(hl.agg.max(mt[expr_field]))
        if min_val is None or max_val is None:
            raise ValueError("No values found for histogram")
        if min_val == max_val:
            min_val -= 1
            max_val += 1

        # Use plain Python integer for n_bins
        expr = hl.float64(mt[expr_field])
        hist = mt.aggregate_entries(hl.agg.hist(expr, min_val, max_val, n_bins))
        x_label = "Expression"

    # Plot using matplotlib
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(
        hist.bin_edges[:-1],
        hist.bin_freq,
        width=np.diff(hist.bin_edges),
        align="edge",
        edgecolor="black",
    )
    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel("Frequency")
    return fig


def visualize_expression_distribution_ad(
    adata: "ad.AnnData",
    n_bins: int = 50,
    log_scale: bool = True,
    title: str = "Expression Value Distribution",
) -> "plt.Figure":
    """
    Visualize the distribution of expression values in an AnnData object.

    Args:
        adata: AnnData object containing expression data in X.
        n_bins: Number of histogram bins.
        log_scale: Whether to apply log1p transformation before plotting.
        title: Plot title.

    Returns:
        Matplotlib figure with the expression histogram.
    """
    import scipy.sparse

    X = adata.X
    if scipy.sparse.issparse(X):
        X = X.toarray()
    values = np.asarray(X).flatten()

    if log_scale:
        values = np.log1p(values)
        x_label = "log1p(Expression)"
    else:
        x_label = "Expression"

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(values, bins=n_bins, edgecolor="black")
    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel("Frequency")
    return fig
