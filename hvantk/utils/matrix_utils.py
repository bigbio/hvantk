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
    min_cells_per_group: int = 1,
) -> pd.DataFrame:
    """Collapse an AnnData expression matrix into a per-group, per-gene summary.

    Parameters
    ----------
    adata
        Expression AnnData (obs = cells/samples, var = genes).
    group_by
        One or more obs columns. Multiple columns are joined with ``_``.
    filter_by
        Optional pre-filter passed to :func:`filter_by_metadata_ad`.
    min_cells_per_group
        Drop groups with fewer cells than this threshold.

    Returns
    -------
    pd.DataFrame
        Long-form with columns ``gene_id, group, mean, fraction_expressed,
        n_cells``.

    Notes
    -----
    Does not mutate *adata*. Uses ``groupby(...).indices`` for O(n_cells)
    group splitting (no per-row index lookups).
    """
    if filter_by:
        adata = filter_by_metadata_ad(adata, filter_by)

    if isinstance(group_by, str):
        group_by = [group_by]

    if len(group_by) == 1:
        labels = adata.obs[group_by[0]].astype(str)
    else:
        labels = adata.obs[group_by[0]].astype(str)
        for col in group_by[1:]:
            labels = labels + "_" + adata.obs[col].astype(str)

    # Local Series — does not touch adata.obs.
    label_series = pd.Series(labels.values, index=np.arange(adata.n_obs))

    gene_ids = adata.var_names.tolist()
    X = adata.X

    records: List[Dict[str, Any]] = []
    for group_name, positions in label_series.groupby(label_series).groups.items():
        n_cells = len(positions)
        if n_cells < min_cells_per_group:
            continue
        sub = X[np.asarray(positions), :]
        means = np.asarray(sub.mean(axis=0)).ravel()
        frac_expr = np.asarray((sub > 0).mean(axis=0)).ravel()
        for j, gid in enumerate(gene_ids):
            records.append(
                {
                    "gene_id": gid,
                    "group": group_name,
                    "mean": float(means[j]),
                    "fraction_expressed": float(frac_expr[j]),
                    "n_cells": n_cells,
                }
            )

    return pd.DataFrame(records)
