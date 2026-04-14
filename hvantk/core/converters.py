"""Hail MatrixTable <-> AnnData converters.

Bridges the Hail and scverse ecosystems:
- ``hail_mt_to_anndata``: MT rows (genes) become AnnData vars,
  MT columns (samples) become AnnData obs.
- ``anndata_to_hail_mt``: AnnData back to Hail MatrixTable.
"""

from __future__ import annotations

import logging
from typing import Optional

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

logger = logging.getLogger(__name__)


def hail_mt_to_anndata(
    mt,
    entry_field: str = "x",
    col_key: Optional[str] = None,
    row_key: Optional[str] = None,
    max_dense_entries: int = 5_000_000,
) -> ad.AnnData:
    """Convert a Hail MatrixTable to AnnData.

    MT rows become AnnData vars (genes/features) and MT columns become
    AnnData obs (samples/cells).  This follows the AnnData convention
    of obs x var.

    Parameters
    ----------
    mt : hl.MatrixTable
        Input MatrixTable.
    entry_field : str
        Name of the entry field to use as the expression matrix.
    col_key : str, optional
        Column key field name. If *None*, uses the first key from
        ``mt.col_key``.
    row_key : str, optional
        Row key field name. If *None*, uses the first key from
        ``mt.row_key``.
    max_dense_entries : int
        Maximum allowed ``n_rows * n_cols`` before refusing dense conversion.

    Returns
    -------
    ad.AnnData
        AnnData object with X as float32.
    """
    # Lazy Hail import (via hail_context to apply numpy compat shim first)
    from hvantk.core.hail_context import hl  # noqa: F811

    # Resolve keys
    if row_key is None:
        row_key = list(mt.row_key)[0]
    if col_key is None:
        col_key = list(mt.col_key)[0]

    logger.info(
        "Converting MT to AnnData (row_key=%s, col_key=%s, entry=%s)",
        row_key,
        col_key,
        entry_field,
    )

    n_rows = mt.count_rows()
    n_cols = mt.count_cols()
    if n_rows * n_cols > max_dense_entries:
        raise ValueError(
            "hail_mt_to_anndata builds a dense matrix in memory and is intended "
            f"for small MatrixTables only. Got {n_rows}x{n_cols} "
            f"({n_rows * n_cols:,} entries), limit is {max_dense_entries:,}."
        )

    # Collect row (var) and col (obs) annotations
    row_df = mt.rows().to_pandas()
    row_df = row_df.set_index(row_key)

    col_df = mt.cols().to_pandas()
    col_df = col_df.set_index(col_key)

    # Collect entries and pivot to dense matrix
    entries_df = mt.select_entries(mt[entry_field]).entries().to_pandas()
    pivot = entries_df.pivot(index=col_key, columns=row_key, values=entry_field)

    # Reindex to match collected row/col order
    pivot = pivot.reindex(index=col_df.index, columns=row_df.index)

    X = pivot.values.astype(np.float32)

    adata = ad.AnnData(X=X, obs=col_df, var=row_df)
    logger.info("Created AnnData: %d obs x %d var", adata.n_obs, adata.n_vars)
    return adata


def anndata_to_hail_mt(
    adata: ad.AnnData,
    row_key: str = "gene_id",
    col_key: str = "sample_id",
    entry_field: str = "x",
):
    """Convert AnnData to a Hail MatrixTable.

    Parameters
    ----------
    adata : ad.AnnData
        Input AnnData object.
    row_key : str
        Name for the row key field (taken from ``adata.var.index``).
    col_key : str
        Name for the column key field (taken from ``adata.obs.index``).
    entry_field : str
        Name for the entry field in the resulting MatrixTable.

    Returns
    -------
    hl.MatrixTable
        MatrixTable keyed by *row_key* and *col_key*.
    """
    from hvantk.core.hail_context import hl

    # Handle sparse X
    X = adata.X
    if sp.issparse(X):
        X = X.toarray()

    n_obs, n_var = X.shape
    obs_names = np.asarray(adata.obs.index, dtype=str)
    var_names = np.asarray(adata.var.index, dtype=str)

    logger.info(
        "Converting AnnData (%d obs x %d var) to MT", adata.n_obs, adata.n_vars
    )

    # Vectorized long-format build (rows ordered obs-major, var-minor).
    long_df = pd.DataFrame(
        {
            col_key: np.repeat(obs_names, n_var),
            row_key: np.tile(var_names, n_obs),
            entry_field: np.asarray(X, dtype="float64").ravel(),
        }
    )
    ht = hl.Table.from_pandas(long_df)
    mt = ht.to_matrix_table(
        row_key=[row_key],
        col_key=[col_key],
    )

    # Annotate rows with adata.var columns (if any beyond the index)
    if len(adata.var.columns) > 0:
        var_df = adata.var.copy()
        var_df[row_key] = var_df.index
        var_ht = hl.Table.from_pandas(var_df).key_by(row_key)
        mt = mt.annotate_rows(**{c: var_ht[mt[row_key]][c] for c in adata.var.columns})

    # Annotate cols with adata.obs columns (if any beyond the index)
    if len(adata.obs.columns) > 0:
        obs_df = adata.obs.copy()
        obs_df[col_key] = obs_df.index
        obs_ht = hl.Table.from_pandas(obs_df).key_by(col_key)
        mt = mt.annotate_cols(**{c: obs_ht[mt[col_key]][c] for c in adata.obs.columns})

    logger.info("Created MT: %s", mt.count())
    return mt
