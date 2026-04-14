"""Expression visualization for AnnData objects."""

from __future__ import annotations

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sp


def visualize_expression_distribution(
    adata: ad.AnnData,
    n_bins: int = 50,
    log_scale: bool = True,
    title: str = "Expression Value Distribution",
) -> plt.Figure:
    """Plot the distribution of expression values from ``adata.X``.

    Parameters
    ----------
    adata : ad.AnnData
        AnnData object with expression values in ``X``.
    n_bins : int
        Number of histogram bins.
    log_scale : bool
        If True, apply ``log1p`` before plotting.
    title : str
        Plot title.

    Returns
    -------
    plt.Figure
        Matplotlib figure containing the histogram.
    """
    X = adata.X
    zero_count = 0
    if sp.issparse(X):
        values = np.asarray(X.data).ravel()
        zero_count = int(X.shape[0] * X.shape[1] - values.size)
    else:
        values = np.asarray(X).ravel()

    if log_scale:
        values = np.log1p(values)
        x_label = "log1p(Expression)"
        zero_value = 0.0
    else:
        x_label = "Expression"
        zero_value = 0.0

    fig, ax = plt.subplots(figsize=(10, 6))
    if zero_count > 0:
        hist_values = np.concatenate(([zero_value], values))
        hist_weights = np.concatenate(([zero_count], np.ones(values.size)))
        ax.hist(hist_values, bins=n_bins, weights=hist_weights, edgecolor="black")
    else:
        ax.hist(values, bins=n_bins, edgecolor="black")
    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel("Frequency")
    return fig
