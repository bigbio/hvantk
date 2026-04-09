"""
Utility functions to convert UCSC Cell Browser datasets to AnnData objects.

This module provides functions to load UCSC metadata and expression matrix files
into pandas DataFrames and AnnData objects for downstream analysis.
"""

from __future__ import annotations

import os

from hvantk.core.constants import UCSC_CELL_ID_COLUMN, UCSC_GENE_COLUMN

__all__ = [
    "load_ucsc_metadata",
    "create_anndata_from_ucsc_matrix",
]


import logging

import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

logger = logging.getLogger(__name__)


def load_ucsc_metadata(
    metadata_path: str,
    sep: str = "\t",
    index_col: int = 0,
    index_name: str = UCSC_CELL_ID_COLUMN,
) -> pd.DataFrame:
    """Read UCSC Cell Browser metadata TSV into a DataFrame.

    Parameters
    ----------
    metadata_path : str
        Path to the metadata file.
    sep : str
        Column delimiter (default tab).
    index_col : int
        Which column to use as the index (default 0).
    index_name : str
        Name to assign to the index (default ``"cell_id"``).

    Returns
    -------
    pd.DataFrame
        Metadata with *index_name* as the index and dots replaced by
        underscores in column names.

    Raises
    ------
    FileNotFoundError
        If *metadata_path* does not exist.
    """
    if not os.path.exists(metadata_path):
        raise FileNotFoundError(f"Metadata file not found: {metadata_path}")

    df = pd.read_csv(metadata_path, sep=sep, index_col=index_col)
    df.index.name = index_name
    df.columns = [c.replace(".", "_") for c in df.columns]
    return df


def create_anndata_from_ucsc_matrix(
    expression_matrix_path: str,
    metadata_df: pd.DataFrame = None,
    gene_column: str = UCSC_GENE_COLUMN,
    delimiter: str = "\t",
    split_gene_field: bool = True,
) -> ad.AnnData:
    """Create an AnnData object from a UCSC Cell Browser expression matrix.

    The input file has genes as rows and cells as columns.  The function
    transposes the data to the AnnData convention (obs = cells, var = genes)
    and stores the expression values as a ``scipy.sparse.csr_matrix`` in
    ``float32``.

    Parameters
    ----------
    expression_matrix_path : str
        Path to the expression TSV (genes x cells).
    metadata_df : pd.DataFrame, optional
        Cell metadata to join into ``adata.obs``.  Index must match cell ids.
    gene_column : str
        Name of the first column containing gene identifiers.
    delimiter : str
        Column delimiter (default tab).
    split_gene_field : bool
        If True, split gene names on ``|`` and keep only the first element.

    Returns
    -------
    ad.AnnData
        AnnData with shape (n_cells, n_genes).

    Raises
    ------
    FileNotFoundError
        If *expression_matrix_path* does not exist.
    """
    if not os.path.exists(expression_matrix_path):
        raise FileNotFoundError(
            f"Expression matrix file not found: {expression_matrix_path}"
        )

    logger.info("Reading expression matrix from %s", expression_matrix_path)
    expr_df = pd.read_csv(
        expression_matrix_path, sep=delimiter, index_col=0
    )

    # gene names are in the index after index_col=0
    gene_names = expr_df.index.astype(str)
    if split_gene_field:
        gene_names = gene_names.str.split("|").str[0]
    expr_df.index = gene_names

    # Transpose: genes x cells -> cells x genes
    X = sparse.csr_matrix(expr_df.values.T.astype(np.float32))

    var = pd.DataFrame(index=pd.Index(gene_names, name=gene_column))
    obs = pd.DataFrame(index=pd.Index(expr_df.columns, name=UCSC_CELL_ID_COLUMN))

    adata = ad.AnnData(X=X, obs=obs, var=var)

    # Join metadata into obs if provided
    if metadata_df is not None:
        # Align metadata to obs index
        common = adata.obs.index.intersection(metadata_df.index)
        if len(common) == 0:
            logger.warning(
                "No overlapping cell ids between expression matrix and metadata"
            )
        meta_aligned = metadata_df.reindex(adata.obs.index)
        for col in meta_aligned.columns:
            adata.obs[col] = meta_aligned[col].values

    return adata
