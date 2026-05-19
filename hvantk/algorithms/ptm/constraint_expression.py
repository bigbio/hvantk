"""Expression-source adapter for the PTM constraint pipeline.

Unifies three backends — Hail MatrixTable, AnnData ``.h5ad``, and pre-computed
tabular files — behind a single contract: a genes x groups wide pandas
``DataFrame`` of aggregate expression.

Downstream analysis (:mod:`hvantk.ptm.constraint`) consumes this single shape
regardless of the origin modality.
"""

from __future__ import annotations

import logging
import os
from typing import Dict, Literal, Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

__all__ = ["load_gene_by_group_matrix"]

Source = Literal["hail-mt", "anndata", "tabular"]
AggFunc = Literal["median", "mean", "median_nonzero"]

_MAX_GROUPS = 200


def load_gene_by_group_matrix(
    source: Source,
    path: str,
    grouping: str,
    aggfunc: AggFunc = "median",
    min_cells_per_group: int = 50,
    entry_field: str = "x",
    gene_key: Optional[str] = None,
) -> pd.DataFrame:
    """Load an expression source into a genes x groups wide ``pd.DataFrame``.

    Parameters
    ----------
    source
        Backend identifier. One of ``"hail-mt"``, ``"anndata"``, ``"tabular"``.
    path
        Filesystem path to the data (``.mt`` directory, ``.h5ad`` file, or a
        parquet/pickle/TSV for tabular).
    grouping
        For Hail MT / AnnData: name of the column metadata field to aggregate
        over. For tabular: ignored (the file is assumed to already be wide).
    aggfunc
        Aggregation function per group. ``median`` (default), ``mean``, or
        ``median_nonzero`` (median of strictly-positive entries; useful for
        sparse single-cell sources).
    min_cells_per_group
        Groups smaller than this are dropped. Only used by Hail MT and AnnData
        backends.
    entry_field
        Hail MT entry field holding the expression value. Default ``"x"``.
    gene_key
        Optional override for the row key used as the gene identifier. If
        ``None`` the backend's default is used (``var_names`` for AnnData,
        first row key for Hail MT).

    Returns
    -------
    pd.DataFrame
        Rows = genes (index named ``gene_id``), columns = groups, values =
        aggregate expression. Guaranteed numeric with no NaN rows.
    """
    if source == "tabular":
        df = _load_from_tabular(path)
    elif source == "anndata":
        df = _load_from_anndata(
            path,
            grouping=grouping,
            aggfunc=aggfunc,
            min_cells_per_group=min_cells_per_group,
        )
    elif source == "hail-mt":
        df = _load_from_hail_mt(
            path,
            grouping=grouping,
            aggfunc=aggfunc,
            min_cells_per_group=min_cells_per_group,
            entry_field=entry_field,
            gene_key=gene_key,
        )
    else:
        raise ValueError(
            f"Unknown source '{source}'. Expected one of "
            "'hail-mt', 'anndata', 'tabular'."
        )

    return _validate_matrix(df)


def _validate_matrix(df: pd.DataFrame) -> pd.DataFrame:
    """Common checks + numeric coercion applied to every backend output."""
    if df.empty:
        raise ValueError("Expression matrix is empty after loading.")

    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.dropna(how="all", axis=0)

    n_groups = df.shape[1]
    if n_groups > _MAX_GROUPS:
        raise ValueError(
            f"Grouping produced {n_groups} groups (> {_MAX_GROUPS}); refusing "
            "to continue. Check that the field is categorical, not numeric."
        )
    if n_groups < 2:
        raise ValueError(
            f"Need >= 2 groups to run constraint analysis; got {n_groups}."
        )

    df = df.fillna(0.0)
    min_val = float(df.values.min()) if df.size else 0.0
    if min_val < 0:
        logger.warning(
            "Expression matrix contains negative values (min=%.3f); clipping to 0.",
            min_val,
        )
        df = df.clip(lower=0.0)

    if df.index.name is None:
        df.index.name = "gene_id"

    logger.info("Expression matrix ready: %d genes x %d groups", *df.shape)
    return df


def _load_from_tabular(path: str) -> pd.DataFrame:
    """Load a pre-aggregated gene x group matrix from disk."""
    ext = os.path.splitext(path)[1].lower()
    logger.info("Loading tabular expression matrix from %s (ext=%s)", path, ext)

    if ext in {".parquet", ".pq"}:
        df = pd.read_parquet(path)
    elif ext in {".pkl", ".pickle"}:
        df = pd.read_pickle(path)
    elif ext in {".tsv", ".txt"}:
        df = pd.read_csv(path, sep="\t", index_col=0)
    elif ext in {".csv"}:
        df = pd.read_csv(path, sep=",", index_col=0)
    else:
        raise ValueError(
            f"Unsupported tabular extension '{ext}' at {path}. "
            "Use .parquet/.pkl/.tsv/.csv."
        )

    if df.index.name is None or df.index.name == 0:
        df.index.name = "gene_id"

    object_cols = df.select_dtypes(include="object").columns
    for col in object_cols:
        df[col] = (
            df[col]
            .astype(str)
            .str.replace(",", ".", regex=False)
            .replace({"": np.nan, "nan": np.nan})
        )
    return df


def _load_from_anndata(
    path: str,
    grouping: str,
    aggfunc: AggFunc,
    min_cells_per_group: int,
) -> pd.DataFrame:
    """Load an AnnData ``.h5ad`` and aggregate to a gene x group matrix.

    ``summarize_expression_ad`` only emits per-group means. For the
    ``median`` / ``median_nonzero`` paths we fall back to a direct numpy
    aggregation over ``adata.X`` so the CLI's ``--expression-metric`` flag
    is honoured instead of silently collapsing to the mean.
    """
    from hvantk.core.anndata_utils import load_anndata
    from hvantk.algorithms.expression.matrix_utils import summarize_expression_ad

    adata = load_anndata(path)
    if grouping not in adata.obs.columns:
        raise KeyError(
            f"Grouping field '{grouping}' not found in adata.obs. "
            f"Available: {list(adata.obs.columns)[:20]}"
        )

    if aggfunc == "mean":
        long_df = summarize_expression_ad(
            adata,
            group_by=grouping,
            min_cells_per_group=min_cells_per_group,
        )
        if long_df.empty:
            raise ValueError(
                f"No groups survived min_cells_per_group={min_cells_per_group}."
            )
        wide = long_df.pivot(index="gene_id", columns="group", values="mean")
        wide.columns.name = None
        return wide

    if aggfunc in {"median", "median_nonzero"}:
        return _aggregate_anndata_direct(adata, grouping, aggfunc, min_cells_per_group)

    raise ValueError(f"Unknown aggfunc '{aggfunc}'.")


def _aggregate_anndata_direct(
    adata,
    grouping: str,
    aggfunc: AggFunc,
    min_cells_per_group: int,
) -> pd.DataFrame:
    """Aggregate AnnData cells → groups with median or median_nonzero per gene."""
    import scipy.sparse as sp

    labels = adata.obs[grouping].astype(str).to_numpy()
    unique_groups = [g for g in sorted(set(labels))]

    gene_ids = list(adata.var_names)
    X = adata.X
    is_sparse = sp.issparse(X)

    columns: Dict[str, np.ndarray] = {}
    for grp in unique_groups:
        mask = labels == grp
        n = int(mask.sum())
        if n < min_cells_per_group:
            logger.info(
                "Dropping group '%s' (n=%d < min_cells_per_group=%d)",
                grp,
                n,
                min_cells_per_group,
            )
            continue
        sub = X[mask, :]
        if is_sparse:
            dense = sub.toarray()
        else:
            dense = np.asarray(sub)

        if aggfunc == "median_nonzero":
            with np.errstate(invalid="ignore"):
                dense = dense.astype(float)
                dense[dense <= 0] = np.nan
                col_values = np.nanmedian(dense, axis=0)
                col_values = np.nan_to_num(col_values, nan=0.0)
        else:
            col_values = np.median(dense, axis=0)

        columns[grp] = col_values

    if not columns:
        raise ValueError(
            f"No groups survived min_cells_per_group={min_cells_per_group}."
        )

    df = pd.DataFrame(columns, index=pd.Index(gene_ids, name="gene_id"))
    return df


def _load_from_hail_mt(
    path: str,
    grouping: str,
    aggfunc: AggFunc,
    min_cells_per_group: int,
    entry_field: str,
    gene_key: Optional[str],
) -> pd.DataFrame:
    """Aggregate a Hail MatrixTable into a gene x group matrix via TSV export.

    Follows the Notebook E/G pattern: ``group_cols_by(...) → aggregate(...) →
    entries() → export(TSV) → pandas`` to avoid OOM from ``mt.to_pandas()``.
    """
    from hvantk.core.hail_context import hl, init_hail

    init_hail()

    logger.info("Reading Hail MatrixTable from %s", path)
    mt = hl.read_matrix_table(path)

    col_fields = set(mt.col)
    if grouping not in col_fields:
        raise KeyError(
            f"Grouping field '{grouping}' not found in mt.col. "
            f"Available: {list(col_fields)[:20]}"
        )

    if gene_key is None:
        row_keys = list(mt.row_key)
        if not row_keys:
            raise ValueError("Hail MT has no row_key; cannot infer gene identifier.")
        gene_key = row_keys[0]

    if entry_field not in set(mt.entry):
        raise KeyError(
            f"Entry field '{entry_field}' not found. "
            f"Available: {list(mt.entry)}"
        )

    entry = mt[entry_field]
    if aggfunc == "mean":
        agg_expr = hl.agg.mean(entry)
    else:
        raise ValueError(
            "Hail MT backend currently supports only aggfunc='mean'. "
            "Use AnnData or tabular backends for median-based aggregation."
        )

    grouped = mt.group_cols_by(mt[grouping]).aggregate(
        _agg=agg_expr,
        _n=hl.agg.count(),
    )

    et = grouped.entries()
    et = et.filter(et._n >= min_cells_per_group)

    tmp_tsv = hl.utils.new_temp_file(extension="tsv")
    logger.info("Exporting aggregated entries to %s", tmp_tsv)
    et = et.select(et[grouping], et._agg)
    et.export(tmp_tsv)

    df_long = pd.read_csv(tmp_tsv, sep="\t")

    grouping_col = grouping
    if grouping_col not in df_long.columns:
        candidate = [c for c in df_long.columns if c.endswith(grouping)]
        if candidate:
            grouping_col = candidate[0]
        else:
            raise RuntimeError(
                f"Could not locate grouping column '{grouping}' in exported TSV; "
                f"columns were {list(df_long.columns)}"
            )

    if gene_key not in df_long.columns:
        candidate = [c for c in df_long.columns if c.endswith(gene_key)]
        if not candidate:
            raise RuntimeError(
                f"Could not locate gene key '{gene_key}' in exported TSV; "
                f"columns were {list(df_long.columns)}"
            )
        gene_col = candidate[0]
    else:
        gene_col = gene_key

    df_long["_agg"] = (
        df_long["_agg"]
        .astype(str)
        .str.replace(",", ".", regex=False)
    )
    df_long["_agg"] = pd.to_numeric(df_long["_agg"], errors="coerce")

    wide = df_long.pivot(index=gene_col, columns=grouping_col, values="_agg")
    wide.index.name = "gene_id"
    wide.columns.name = None
    return wide
