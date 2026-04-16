"""AnnData-based expression analysis utilities.

Provides describe, filter, and per-group summarization functions that operate
directly on AnnData objects. Hail-MatrixTable equivalents were removed when
expression I/O migrated to AnnData (.h5ad).
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional, Union

import anndata as ad
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

__all__ = [
    "describe_expression_ad",
    "filter_by_metadata_ad",
    "summarize_expression_ad",
]


def describe_expression_ad(adata: ad.AnnData) -> Dict[str, Any]:
    """Return a summary dict describing an AnnData expression object.

    Prefers the cached ``adata.uns["column_summary"]`` (set by
    :func:`hvantk.core.anndata_utils.annotate_column_summary_ad`); otherwise
    computes a lightweight summary from ``adata.obs`` on the fly.

    Returns
    -------
    dict
        ``{"n_obs", "n_vars", "fields"}`` where each field entry has
        ``{"name", "dtype"}`` plus ``min``/``max`` for numeric or
        ``n_unique`` for categorical.
    """
    fields_info: List[Dict[str, Any]] = []

    if "column_summary" in adata.uns:
        for field_name, info in sorted(adata.uns["column_summary"].items()):
            entry: Dict[str, Any] = {"name": field_name, "dtype": info["dtype"]}
            if info["dtype"] == "categorical":
                entry["n_unique"] = info.get("n_unique")
            else:
                entry["min"] = info.get("min")
                entry["max"] = info.get("max")
            fields_info.append(entry)
    else:
        for col in adata.obs.columns:
            series = adata.obs[col]
            if pd.api.types.is_numeric_dtype(series):
                fields_info.append(
                    {
                        "name": col,
                        "dtype": "numeric",
                        "min": float(np.nanmin(series.values)),
                        "max": float(np.nanmax(series.values)),
                    }
                )
            else:
                fields_info.append(
                    {
                        "name": col,
                        "dtype": "categorical",
                        "n_unique": int(series.nunique()),
                    }
                )

    return {"n_obs": adata.n_obs, "n_vars": adata.n_vars, "fields": fields_info}


def filter_by_metadata_ad(
    adata: ad.AnnData,
    filters: Dict[str, Union[str, List[str]]],
) -> ad.AnnData:
    """Return a copy of *adata* keeping only obs matching all *filters*.

    Parameters
    ----------
    filters
        Mapping ``obs_column -> value`` or ``obs_column -> [values]``.
    """
    mask = np.ones(adata.n_obs, dtype=bool)
    for col, values in filters.items():
        if not isinstance(values, list):
            values = [values]
        mask &= adata.obs[col].isin(values).values
    return adata[mask].copy()


def summarize_expression_ad(
    adata: ad.AnnData,
    group_by: Union[str, List[str]],
    filter_by: Optional[Dict[str, Union[str, List[str]]]] = None,
    min_cells_per_group: int = 10,
) -> ad.AnnData:
    """Collapse an AnnData expression matrix into a per-group, per-gene AnnData.

    Thin wrapper around :func:`scanpy.get.aggregate` that also derives
    ``fraction_expressed`` and attaches per-group cell counts.

    Parameters
    ----------
    adata
        Expression AnnData (obs = cells/samples, var = genes).
    group_by
        One or more obs columns to group by. Multi-column groupings produce
        ``obs_names`` joined with ``_`` (matches ``scanpy.get.aggregate``).
    filter_by
        Optional pre-filter passed to :func:`filter_by_metadata_ad`.
    min_cells_per_group
        Drop groups with fewer than this many cells.

    Returns
    -------
    ad.AnnData
        Shape ``(n_groups, n_genes)`` with layers ``mean``, ``sum``,
        ``count_nonzero``, ``fraction_expressed``. ``obs["n_cells"]`` stores
        the per-group cell count.
    """
    import scanpy as sc

    if filter_by:
        adata = filter_by_metadata_ad(adata, filter_by)

    by = [group_by] if isinstance(group_by, str) else list(group_by)

    agg = sc.get.aggregate(
        adata,
        by=by,
        func=["mean", "sum", "count_nonzero"],
    )

    # sc.get.aggregate does not stash per-group cell counts; compute from the
    # pre-aggregate adata, joining group_by columns with "_" to match
    # agg.obs_names.
    group_labels = adata.obs[by].astype(str).agg("_".join, axis=1)
    n_cells = group_labels.value_counts().reindex(agg.obs_names).astype(int)

    agg.obs["n_cells"] = n_cells.values
    agg.layers["fraction_expressed"] = (
        np.asarray(agg.layers["count_nonzero"]) / n_cells.values[:, None]
    )

    keep = agg.obs["n_cells"].values >= min_cells_per_group
    return agg[keep].copy()
