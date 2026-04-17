"""
Utility functions to convert UCSC Cell Browser datasets to AnnData objects.

This module provides functions to load UCSC metadata and expression matrix files
into pandas DataFrames and AnnData objects for downstream analysis.
"""

from __future__ import annotations

import gzip
import logging
import os
from typing import Iterator

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from anndata.io import sparse_dataset, write_elem
from scipy import sparse

from hvantk.core.constants import UCSC_CELL_ID_COLUMN, UCSC_GENE_COLUMN

__all__ = [
    "load_ucsc_metadata",
    "create_anndata_from_ucsc_matrix",
    "build_ucsc_atlas_backed",
]

logger = logging.getLogger(__name__)


def _open_text(path: str):
    """Open a plain or gzipped text file for reading."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def _iter_ucsc_rows(
    expression_matrix_path: str,
    delimiter: str = "\t",
    split_gene_field: bool = True,
) -> tuple[list[str], Iterator[tuple[str, np.ndarray]]]:
    """Open a UCSC expression TSV (plain or gzipped) and return
    ``(cell_ids, row_iterator)``.

    ``cell_ids`` is the header's cell-column list, extracted eagerly.
    ``row_iterator`` yields ``(gene_name, row_values_float32)`` per data
    line and closes the underlying handle when exhausted.

    The caller must iterate the returned generator to completion or call
    its ``.close()`` method; otherwise the underlying file handle only
    closes at generator finalization, which is GC-timing-dependent.

    Raises
    ------
    FileNotFoundError
        If ``expression_matrix_path`` does not exist.
    ValueError
        If the file is empty or a row has the wrong number of values.
    """
    if not os.path.exists(expression_matrix_path):
        raise FileNotFoundError(
            f"Expression matrix file not found: {expression_matrix_path}"
        )

    fh = _open_text(expression_matrix_path)
    try:
        header = fh.readline().rstrip("\n").rstrip("\r")
        if not header:
            fh.close()
            raise ValueError(
                f"Expression matrix {expression_matrix_path!r} is empty."
            )
        header_fields = header.split(delimiter)
        cell_ids = header_fields[1:]
        n_cells = len(cell_ids)
    except Exception:
        fh.close()
        raise

    def _rows() -> Iterator[tuple[str, np.ndarray]]:
        try:
            for line in fh:
                line = line.rstrip("\n").rstrip("\r")
                if not line:
                    continue
                tab = line.find(delimiter)
                if tab < 0:
                    continue
                gene = line[:tab]
                if split_gene_field:
                    gene = gene.split("|", 1)[0]
                row = np.fromstring(
                    line[tab + 1:], sep=delimiter, dtype=np.float32
                )
                if row.shape[0] != n_cells:
                    raise ValueError(
                        f"Row for gene {gene!r} has {row.shape[0]} parseable "
                        f"float values; expected {n_cells}. This usually "
                        f"means a short row or a non-numeric token in the row."
                    )
                yield gene, row
        finally:
            fh.close()

    return cell_ids, _rows()


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
    streams the file row-by-row via ``gzip``/``open`` and parses each row
    into a ``float32`` numpy array with ``np.fromstring``. Every
    ``chunk_size`` gene rows are flushed into a ``scipy.sparse.csr_matrix``
    chunk; the chunks are then vstacked and transposed to the AnnData
    convention (obs = cells, var = genes). Peak RAM during streaming stays
    roughly ``chunk_size × n_cells × 4 B`` plus the growing sparse result.

    A row-oriented parser is used here rather than ``pandas.read_csv`` or
    ``pyarrow.csv``: single-cell matrices from UCSC can have ~500k columns,
    and both pandas' C engine (on gzip streams) and pyarrow's
    column-oriented batches become orders-of-magnitude slower at that width.

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
        If a row fails float32 conversion (names the offending gene).
    """
    logger.info(
        "Streaming expression matrix from %s (chunk_size=%d)",
        expression_matrix_path,
        chunk_size,
    )

    cell_ids, row_iter = _iter_ucsc_rows(
        expression_matrix_path,
        delimiter=delimiter,
        split_gene_field=split_gene_field,
    )

    gene_name_chunks: list[np.ndarray] = []
    sparse_chunks: list[sparse.csr_matrix] = []
    buf_rows: list[np.ndarray] = []
    buf_genes: list[str] = []
    n_seen = 0
    n_chunks = 0

    def _flush() -> None:
        nonlocal n_chunks
        if not buf_rows:
            return
        sparse_chunks.append(sparse.csr_matrix(np.vstack(buf_rows)))
        gene_name_chunks.append(np.asarray(buf_genes))
        buf_rows.clear()
        buf_genes.clear()
        n_chunks += 1
        if n_chunks % 10 == 0:
            logger.info("  streamed %d genes", n_seen)

    for gene, row in row_iter:
        buf_rows.append(row)
        buf_genes.append(gene)
        n_seen += 1
        if len(buf_rows) >= chunk_size:
            _flush()
    _flush()

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
            raise ValueError(
                f"No overlapping cell ids between expression matrix "
                f"{expression_matrix_path!r} and metadata."
            )
        meta_aligned = metadata_df.reindex(adata.obs.index)
        for col in meta_aligned.columns:
            adata.obs[col] = meta_aligned[col].values

    logger.info("Built AnnData: %d cells × %d genes", adata.shape[0], adata.shape[1])
    return adata


def build_ucsc_atlas_backed(
    expression_matrix_path: str,
    output_path: str,
    metadata_df: pd.DataFrame | None = None,
    gene_column: str = UCSC_GENE_COLUMN,
    delimiter: str = "\t",
    split_gene_field: bool = True,
    column_batch: int = 64,
    overwrite: bool = False,
    uns: dict | None = None,
) -> str:
    """Stream-build an AnnData .h5ad file on disk, appending one batch of
    genes (CSC columns) at a time. Never materializes the full
    ``cells × genes`` matrix.

    Parameters
    ----------
    expression_matrix_path : str
        Path to the UCSC expression TSV (plain or gzipped).
    output_path : str
        Destination ``.h5ad`` path.
    metadata_df : pd.DataFrame, optional
        Cell metadata; reindexed to expression header ``cell_ids``.
    gene_column : str
        Name for the ``var`` index.
    delimiter : str
        Column delimiter in the expression TSV.
    split_gene_field : bool
        If True, split gene ids on ``|`` and keep the first element.
    column_batch : int
        Number of gene columns buffered before a CSC append to disk.
        Peak per-batch RAM ≈ ``column_batch × n_cells × 4 B``
        (≈16 MB at the default for a 520k-cell atlas).
    overwrite : bool
        If False and ``output_path`` exists, raise ``FileExistsError``.
    uns : dict, optional
        Opaque dict written to the h5ad's ``uns`` group (e.g. provenance
        metadata). Callers remain responsible for the dict's structure.

    Returns
    -------
    str
        ``output_path``.
    """
    if os.path.exists(output_path) and not overwrite:
        raise FileExistsError(
            f"{output_path!r} already exists; pass overwrite=True to replace."
        )

    cell_ids, row_iter = _iter_ucsc_rows(
        expression_matrix_path,
        delimiter=delimiter,
        split_gene_field=split_gene_field,
    )
    n_cells = len(cell_ids)

    # Check metadata overlap BEFORE streaming: if zero overlap, raise now
    # rather than after writing gigabytes of X to disk only to leave a
    # half-written h5ad that blocks re-runs.
    if metadata_df is not None:
        if len(pd.Index(cell_ids).intersection(metadata_df.index)) == 0:
            raise ValueError(
                f"No overlapping cell ids between expression matrix "
                f"{expression_matrix_path!r} and metadata."
            )

    logger.info(
        "Backed-write atlas → %s (n_cells=%d, column_batch=%d)",
        output_path, n_cells, column_batch,
    )

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)

    # h5py.File(..., "w") truncates any existing file, so no pre-remove needed.
    f = h5py.File(output_path, "w")
    try:
        # Seed X as an empty (n_cells, 0) CSC; pass indptr dtype via
        # dataset_kwargs so large atlases don't overflow int32.
        seed = sparse.csc_matrix((n_cells, 0), dtype=np.float32)
        write_elem(f, "X", seed, dataset_kwargs={"indptr_dtype": "int64"})
        X = sparse_dataset(f["X"])

        gene_names: list[str] = []
        buf_cols: list[np.ndarray] = []
        buf_genes: list[str] = []
        n_appended = 0
        n_batches = 0

        def _flush() -> None:
            nonlocal n_appended, n_batches
            if not buf_cols:
                return
            # Stack column-vectors into (n_cells, batch_width) CSC block.
            dense_block = np.column_stack(buf_cols)  # (n_cells, batch_width)
            X.append(sparse.csc_matrix(dense_block))
            gene_names.extend(buf_genes)
            n_appended += len(buf_cols)
            n_batches += 1
            buf_cols.clear()
            buf_genes.clear()
            if n_batches % 20 == 0:
                logger.info("  appended %d genes", n_appended)

        for gene, row in row_iter:
            buf_cols.append(row)
            buf_genes.append(gene)
            if len(buf_cols) >= column_batch:
                _flush()
        _flush()

        if n_appended == 0:
            raise ValueError(
                f"Expression matrix {expression_matrix_path!r} contained no data rows."
            )

        # obs: cell_ids + optional metadata join (reject zero-overlap).
        obs = pd.DataFrame(index=pd.Index(cell_ids, name=UCSC_CELL_ID_COLUMN))
        if metadata_df is not None:
            common = obs.index.intersection(metadata_df.index)
            if len(common) == 0:
                raise ValueError(
                    f"No overlapping cell ids between expression matrix "
                    f"{expression_matrix_path!r} and metadata."
                )
            meta_aligned = metadata_df.reindex(obs.index)
            for col in meta_aligned.columns:
                obs[col] = meta_aligned[col].values

        var = pd.DataFrame(index=pd.Index(gene_names, name=gene_column))

        write_elem(f, "obs", obs)
        write_elem(f, "var", var)
        if uns is not None:
            write_elem(f, "uns", uns)
    finally:
        f.close()

    logger.info(
        "Backed atlas written: %d cells × %d genes → %s",
        n_cells, n_appended, output_path,
    )
    return output_path
