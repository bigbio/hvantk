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
    chunk_size: int = 500,
) -> ad.AnnData:
    """Create an AnnData object from a UCSC Cell Browser expression matrix.

    The input file has genes as rows and cells as columns. This function
    streams the file in chunks of ``chunk_size`` genes, converts each chunk
    to a ``scipy.sparse.csr_matrix``, and vstacks them, so peak RAM stays
    roughly ``chunk_size × n_cells × 4 B`` plus the growing sparse result.
    The stacked matrix is transposed to the AnnData convention
    (obs = cells, var = genes) and stored as ``float32`` CSR.

    Parameters
    ----------
    expression_matrix_path : str
        Path to the expression TSV or `.tsv.gz` (genes x cells).
    metadata_df : pd.DataFrame, optional
        Cell metadata to join into ``adata.obs``. Index must match cell ids.
    gene_column : str
        Name of the first column containing gene identifiers.
    delimiter : str
        Column delimiter (default tab).
    split_gene_field : bool
        If True, split gene names on ``|`` and keep only the first element.
    chunk_size : int
        Gene rows per streaming chunk. Lower ⇒ less peak RAM, slower.

    Returns
    -------
    ad.AnnData
        AnnData with shape (n_cells, n_genes), sparse CSR X.

    Raises
    ------
    FileNotFoundError
        If *expression_matrix_path* does not exist.
    ValueError
        If a chunk fails float32 conversion (names the offending gene row).
    """
    if not os.path.exists(expression_matrix_path):
        raise FileNotFoundError(
            f"Expression matrix file not found: {expression_matrix_path}"
        )

    logger.info(
        "Streaming expression matrix from %s (chunk_size=%d)",
        expression_matrix_path,
        chunk_size,
    )

    gene_name_chunks = []
    sparse_chunks = []
    cell_ids = None
    n_seen = 0

    reader = pd.read_csv(
        expression_matrix_path,
        sep=delimiter,
        header=0,
        index_col=0,
        chunksize=chunk_size,
    )
    for chunk_idx, chunk in enumerate(reader):
        if cell_ids is None:
            cell_ids = list(chunk.columns)

        chunk_genes = chunk.index.astype(str)
        if split_gene_field:
            chunk_genes = chunk_genes.str.split("|").str[0]
        gene_name_chunks.append(np.asarray(chunk_genes))

        try:
            chunk_values = chunk.to_numpy(dtype=np.float32, copy=False)
        except (ValueError, TypeError) as exc:
            bad_gene = chunk_genes[0] if len(chunk_genes) else "<unknown>"
            raise ValueError(
                f"Non-numeric value in chunk starting at gene {bad_gene!r}: {exc}"
            ) from exc

        sparse_chunks.append(sparse.csr_matrix(chunk_values))
        n_seen += len(chunk_genes)
        if (chunk_idx + 1) % 10 == 0:
            logger.info("  streamed %d genes", n_seen)

    if not sparse_chunks:
        raise ValueError(
            f"Expression matrix {expression_matrix_path!r} contained no data rows."
        )

    # Genes-x-cells sparse → cells-x-genes AnnData convention.
    X = sparse.vstack(sparse_chunks, format="csr").T.tocsr()
    gene_names = np.concatenate(gene_name_chunks)

    var = pd.DataFrame(index=pd.Index(gene_names, name=gene_column))
    obs = pd.DataFrame(index=pd.Index(cell_ids, name=UCSC_CELL_ID_COLUMN))

    adata = ad.AnnData(X=X, obs=obs, var=var)

    if metadata_df is not None:
        common = adata.obs.index.intersection(metadata_df.index)
        if len(common) == 0:
            logger.warning(
                "No overlapping cell ids between expression matrix and metadata"
            )
        meta_aligned = metadata_df.reindex(adata.obs.index)
        for col in meta_aligned.columns:
            adata.obs[col] = meta_aligned[col].values

    logger.info("Built AnnData: %d cells × %d genes", adata.shape[0], adata.shape[1])
    return adata
