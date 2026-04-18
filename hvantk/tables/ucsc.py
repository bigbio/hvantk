"""
Utility functions to convert UCSC Cell Browser datasets to AnnData objects.

This module provides functions to load UCSC metadata and expression matrix files
into pandas DataFrames and AnnData objects for downstream analysis.
"""

from __future__ import annotations

import gzip
import logging
import os
import tempfile
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
    "summarize_ucsc_streaming",
    "finalize_partial_atlas",
    "coerce_obs_for_h5ad",
]

logger = logging.getLogger(__name__)


def _open_text(path: str):
    """Open a plain, gzipped, or bgzipped text file for reading.

    ``.bgz`` (bgzip) is the block-gzip format used by htslib/Hail. It is
    gzip-compatible at the decompressor level, so ``gzip.open`` handles it
    transparently.
    """
    if path.endswith((".gz", ".bgz")):
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


def _iter_ucsc_gene_names_only(
    expression_matrix_path: str,
    delimiter: str = "\t",
    split_gene_field: bool = True,
) -> tuple[list[str], list[str]]:
    """Fast pass over a UCSC expression TSV that reads the header and the
    first column of every data row, skipping the float parse entirely.

    Returns ``(cell_ids, gene_names)``. Useful for recovery paths where the
    per-cell float data is already on disk and only the row labels are
    needed.
    """
    if not os.path.exists(expression_matrix_path):
        raise FileNotFoundError(
            f"Expression matrix file not found: {expression_matrix_path}"
        )
    fh = _open_text(expression_matrix_path)
    try:
        header = fh.readline().rstrip("\n").rstrip("\r")
        if not header:
            raise ValueError(
                f"Expression matrix {expression_matrix_path!r} is empty."
            )
        cell_ids = header.split(delimiter)[1:]
        gene_names: list[str] = []
        for line in fh:
            tab = line.find(delimiter)
            if tab < 0:
                continue
            gene = line[:tab]
            if split_gene_field:
                gene = gene.split("|", 1)[0]
            gene_names.append(gene)
    finally:
        fh.close()
    return cell_ids, gene_names


def coerce_obs_for_h5ad(obs: pd.DataFrame) -> pd.DataFrame:
    """Normalize a ``DataFrame`` so it can be written to h5ad.

    ``anndata``'s vlen-string HDF5 writer raises ``TypeError`` when an
    ``object``-dtype column contains NaN (which pandas stores as ``float``)
    mixed with strings — the implicit ``float → str`` conversion never
    happens. Coerce every ``object`` column through
    ``.where(notna, "").astype(str)`` so NaN becomes the empty string and
    every remaining value is a Python ``str``. Numeric, boolean,
    categorical, and datetime columns are left untouched.

    Returns a copy; does not mutate the input.
    """
    obs = obs.copy()
    for col in obs.select_dtypes(include=["object"]).columns:
        obs[col] = obs[col].where(obs[col].notna(), "").astype(str)
    return obs


def _probe_write_elem(elem, key: str) -> None:
    """Write *elem* to a throwaway HDF5 file under *key* to verify it
    serializes. Raises the underlying error (unchanged) so the caller
    learns *before* committing to a long-running stream that a metadata
    column won't serialize.

    Cheap — writes to a temp file that is removed on exit.
    """
    with tempfile.NamedTemporaryFile(suffix=".h5", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        with h5py.File(tmp_path, "w") as probe:
            write_elem(probe, key, elem)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


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

    # Guard against duplicate cell_id rows (seen in some UCSC/source TSVs with
    # unquoted embedded newlines in string fields, e.g. Asp_2019 celltype).
    # Downstream reindex/join calls raise a cryptic "cannot reindex on an axis
    # with duplicate labels" — dedup here with a clear warning instead.
    n_dup = int(df.index.duplicated().sum())
    if n_dup > 0:
        n_before = len(df)
        df = df.loc[~df.index.duplicated(keep="first")]
        logger.warning(
            "Metadata %s has %d duplicate %s row(s) (%d rows \u2192 %d unique). "
            "Keeping first occurrence of each %s. This often indicates an "
            "upstream parsing issue (e.g. unquoted newlines in string fields); "
            "dedupe the source file if this masks real inconsistencies.",
            metadata_path, n_dup, index_name, n_before, len(df), index_name,
        )
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

    # Build the full obs DataFrame now (before the stream) so we can
    # probe its h5ad-serializability in <1s. This catches cases like
    # object columns with NaN mixed with strings — which anndata's
    # vlen-string writer cannot coerce — before we commit hours of I/O
    # to the X stream only to fail at the closing write_elem.
    obs = pd.DataFrame(index=pd.Index(cell_ids, name=UCSC_CELL_ID_COLUMN))
    if metadata_df is not None:
        meta_aligned = metadata_df.reindex(obs.index)
        for col in meta_aligned.columns:
            obs[col] = meta_aligned[col].values
    obs = coerce_obs_for_h5ad(obs)
    _probe_write_elem(obs, "obs")

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


def finalize_partial_atlas(
    output_path: str,
    expression_matrix_path: str,
    metadata_df: pd.DataFrame | None = None,
    gene_column: str = UCSC_GENE_COLUMN,
    delimiter: str = "\t",
    split_gene_field: bool = True,
    uns: dict | None = None,
) -> str:
    """Finish writing ``/obs``, ``/var``, ``/uns`` into a partially-written
    backed ``.h5ad`` produced by :func:`build_ucsc_atlas_backed`.

    Recovers files whose ``/X`` group was streamed to completion but whose
    metadata write failed or was interrupted (crash, out-of-memory on
    the post-stream serialization, ``write_elem`` type error, etc.). The
    expensive ``X`` data is preserved; we re-parse only the source gzip's
    row labels (first column of each line — minutes, not hours) to
    rebuild ``var``, then write a fresh ``obs`` (coerced via
    :func:`coerce_obs_for_h5ad`) and optional ``uns`` provenance.

    Parameters
    ----------
    output_path : str
        Partial ``.h5ad`` produced by ``build_ucsc_atlas_backed``. Must
        contain a complete ``/X`` group (shape matches source dims).
    expression_matrix_path : str
        The same source TSV/gz that produced the partial file. Used to
        re-derive ``cell_ids`` (header) and ``var`` (first column of each
        row). Source cell count and gene count must match ``/X``'s
        ``shape`` attribute.
    metadata_df : pd.DataFrame, optional
        Cell metadata to rebuild ``/obs`` from. Same shape + coercion as
        ``build_ucsc_atlas_backed``.
    gene_column, delimiter, split_gene_field : see ``build_ucsc_atlas_backed``.
    uns : dict, optional
        Provenance dict to write as ``/uns``.

    Raises
    ------
    ValueError
        If ``/X`` is missing, or source cell/gene counts disagree with
        ``/X``'s shape (the source is not the same one that built this
        partial file — refuse to clobber with mismatched labels).
    """
    with h5py.File(output_path, "r") as f:
        if "X" not in f:
            raise ValueError(
                f"{output_path!r} has no /X group; nothing to recover — "
                f"re-run build_ucsc_atlas_backed from scratch."
            )
        x_shape = tuple(int(s) for s in f["X"].attrs["shape"])
    n_cells_x, n_genes_x = x_shape
    logger.info(
        "Finalizing partial atlas at %s (/X shape=%d × %d)",
        output_path, n_cells_x, n_genes_x,
    )

    cell_ids, gene_names = _iter_ucsc_gene_names_only(
        expression_matrix_path,
        delimiter=delimiter,
        split_gene_field=split_gene_field,
    )
    if len(cell_ids) != n_cells_x:
        raise ValueError(
            f"Source cell count ({len(cell_ids)}) does not match the "
            f"partial atlas /X rows ({n_cells_x}). The partial file was "
            f"likely built from a different source — refusing to finalize."
        )
    if len(gene_names) != n_genes_x:
        raise ValueError(
            f"Source gene count ({len(gene_names)}) does not match the "
            f"partial atlas /X columns ({n_genes_x}). The source and "
            f"partial file disagree — refusing to finalize."
        )

    obs = pd.DataFrame(index=pd.Index(cell_ids, name=UCSC_CELL_ID_COLUMN))
    if metadata_df is not None:
        if len(obs.index.intersection(metadata_df.index)) == 0:
            raise ValueError(
                f"No overlapping cell ids between {expression_matrix_path!r} "
                f"header and the provided metadata_df."
            )
        meta_aligned = metadata_df.reindex(obs.index)
        for col in meta_aligned.columns:
            obs[col] = meta_aligned[col].values
    obs = coerce_obs_for_h5ad(obs)
    var = pd.DataFrame(index=pd.Index(gene_names, name=gene_column))

    with h5py.File(output_path, "a") as f:
        for key in ("obs", "var", "uns"):
            if key in f:
                del f[key]
        write_elem(f, "obs", obs)
        write_elem(f, "var", var)
        if uns is not None:
            write_elem(f, "uns", uns)

    logger.info(
        "Finalized partial atlas at %s (%d cells × %d genes)",
        output_path, n_cells_x, n_genes_x,
    )
    return output_path


def summarize_ucsc_streaming(
    expression_matrix_path: str,
    metadata_df: pd.DataFrame,
    group_by: str | list[str],
    filter_by: dict[str, str | list[str]] | None = None,
    min_cells_per_group: int = 10,
    gene_column: str = UCSC_GENE_COLUMN,
    delimiter: str = "\t",
    split_gene_field: bool = True,
) -> ad.AnnData:
    """Stream a UCSC expression matrix and aggregate per-group per-gene
    statistics in one pass.

    Returns an AnnData of shape ``(n_groups, n_genes)`` with ``X`` set to
    the per-group mean (float32) and ``layers`` containing ``sum``,
    ``count_nonzero``, ``fraction_expressed``. ``obs`` includes the
    original group-by columns plus ``n_cells``.

    Memory: one ``float32[n_cells]`` row buffer + accumulators of size
    ``n_groups × n_genes × 16 B`` (sum float64 + count_nz int64). Independent
    of ``n_cells`` past the single-row buffer.
    """
    by = [group_by] if isinstance(group_by, str) else list(group_by)

    # Apply filter on the DataFrame before factorizing groups.
    work = metadata_df.copy()
    if filter_by:
        for field, value in filter_by.items():
            if field not in work.columns:
                raise ValueError(
                    f"filter_by field {field!r} not in metadata columns: "
                    f"{sorted(work.columns)}"
                )
            if isinstance(value, (list, tuple, set)):
                work = work[work[field].isin(list(value))]
            else:
                work = work[work[field] == value]
        if work.empty:
            raise ValueError("filter_by produced zero cells.")

    for col in by:
        if col not in work.columns:
            raise ValueError(
                f"group_by column {col!r} not in metadata columns: "
                f"{sorted(work.columns)}"
            )

    cell_ids, row_iter = _iter_ucsc_rows(
        expression_matrix_path,
        delimiter=delimiter,
        split_gene_field=split_gene_field,
    )
    n_cells = len(cell_ids)
    logger.info(
        "Fused aggregate stream → %s (n_cells=%d, group_by=%s)",
        expression_matrix_path, n_cells, by,
    )

    # Align metadata to expression header order, keep only cells present
    # in both, and drop NaN-in-group rows.
    aligned = work.reindex(cell_ids)
    keep_mask = aligned[by].notna().all(axis=1).to_numpy()
    if not keep_mask.any():
        raise ValueError(
            "No cells overlap between metadata (post-filter) and expression header."
        )
    # Count cells present in post-filter metadata that were dropped due to
    # NaN in any group-by column (excludes cells absent from metadata entirely).
    n_dropped_nan = int((~keep_mask & aligned.index.isin(work.index)).sum())
    if n_dropped_nan:
        example = cell_ids[int(np.where(~keep_mask)[0][0])]
        logger.info(
            "Dropped %d cells with NaN in group-by columns (example: %s).",
            n_dropped_nan, example,
        )

    labels = (
        aligned.loc[keep_mask, by].astype(str).agg("_".join, axis=1)
        if len(by) > 1
        else aligned.loc[keep_mask, by[0]].astype(str)
    )
    codes, uniques = pd.factorize(labels, sort=True)
    n_groups = len(uniques)

    # group_idx per cell in expression-header order, -1 for dropped cells.
    group_idx = np.full(n_cells, -1, dtype=np.int64)
    group_idx[keep_mask] = codes
    valid = group_idx >= 0
    n_cells_per_group = np.bincount(group_idx[valid], minlength=n_groups).astype(np.int64)

    sum_matrix = np.zeros((n_groups, 0), dtype=np.float64)
    count_nz_matrix = np.zeros((n_groups, 0), dtype=np.int64)
    gene_names: list[str] = []

    # Pre-size buffers in chunks to avoid per-gene resize; column-extend
    # accumulators lazily.
    BLOCK = 1024
    sum_block = np.zeros((n_groups, BLOCK), dtype=np.float64)
    count_block = np.zeros((n_groups, BLOCK), dtype=np.int64)
    block_fill = 0

    def _flush_block():
        nonlocal sum_matrix, count_nz_matrix, block_fill
        if block_fill == 0:
            return
        sum_matrix = np.concatenate([sum_matrix, sum_block[:, :block_fill]], axis=1)
        count_nz_matrix = np.concatenate(
            [count_nz_matrix, count_block[:, :block_fill]], axis=1
        )
        block_fill = 0

    # Heartbeat every PROGRESS_EVERY genes so long-running streams show
    # progress; matches the cadence of the backed builder.
    PROGRESS_EVERY = 5000
    n_streamed = 0

    for gene, row in row_iter:
        row_valid = row[valid]
        sum_block[:, block_fill] = np.bincount(
            group_idx[valid], weights=row_valid, minlength=n_groups
        )
        count_block[:, block_fill] = np.bincount(
            group_idx[valid],
            weights=(row_valid != 0).astype(np.float64),
            minlength=n_groups,
        ).astype(np.int64)
        gene_names.append(gene)
        block_fill += 1
        n_streamed += 1
        if block_fill == BLOCK:
            _flush_block()
        if n_streamed % PROGRESS_EVERY == 0:
            logger.info("  streamed %d genes", n_streamed)
    _flush_block()

    if len(gene_names) == 0:
        raise ValueError(
            f"Expression matrix {expression_matrix_path!r} contained no data rows."
        )

    n_cells_col = n_cells_per_group[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        mean = np.where(n_cells_col > 0, sum_matrix / n_cells_col, 0.0).astype(np.float32)
        fraction_expressed = np.where(
            n_cells_col > 0, count_nz_matrix / n_cells_col, 0.0
        ).astype(np.float32)

    # Build per-field obs columns from the source DataFrame rather than
    # splitting the composite label — values containing "_" would corrupt
    # the split. Take one representative row per code from the pre-aggregate
    # metadata and reindex to the code range (0..n_groups-1), which matches
    # ``uniques`` order because pd.factorize(sort=True) returns codes that
    # index into the sorted uniques.
    by_src = aligned.loc[keep_mask, by].copy()
    by_src["_code"] = codes
    per_group = (
        by_src.groupby("_code", sort=True).first().reindex(np.arange(n_groups))
    )

    obs = pd.DataFrame({"n_cells": n_cells_per_group}, index=pd.Index(uniques))
    for col in by:
        obs[col] = per_group[col].values

    var = pd.DataFrame(index=pd.Index(gene_names, name=gene_column))

    keep_groups = obs["n_cells"].to_numpy() >= min_cells_per_group
    if not keep_groups.any():
        # Build a short report for the error message.
        report = "; ".join(
            f"{u}={int(n)}" for u, n in zip(uniques, n_cells_per_group)
        )
        raise ValueError(
            f"All groups fell below min_cells_per_group={min_cells_per_group}. "
            f"Observed: {report}"
        )

    adata = ad.AnnData(
        X=mean[keep_groups],
        obs=obs[keep_groups].copy(),
        var=var,
        layers={
            "sum": sum_matrix[keep_groups].astype(np.float32),
            "count_nonzero": count_nz_matrix[keep_groups].astype(np.int64),
            "fraction_expressed": fraction_expressed[keep_groups],
            "mean": mean[keep_groups],
        },
    )
    logger.info(
        "Aggregated AnnData: %d groups × %d genes (dropped %d groups below min_cells)",
        adata.n_obs, adata.n_vars, int((~keep_groups).sum()),
    )
    return adata
