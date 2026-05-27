"""AnnData utility functions for column summaries.

I/O functions (save_anndata, load_anndata) live in
:mod:`hvantk.core.io.anndata_io` to honor the intra-core rule that
core/models must not perform disk I/O.

The legacy ``build_anndata_metadata`` writer (and its Hail Table sibling
``build_table_metadata``) was retired alongside the rest of the
``hvantk_metadata`` globals mechanism. Provenance now lives exclusively
on the ``Artifact.provenance`` field stamped by ``run_builder_for_spec``.
"""

import logging
from typing import Any, Dict

import anndata as ad
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def annotate_column_summary_ad(
    adata: ad.AnnData,
    max_levels: int = 100,
    top_n_levels: int = 10,
) -> None:
    """Compute obs column metadata summary and store in ``adata.uns``.

    For each column in ``adata.obs``:
    - Numeric columns: dtype, min, max, mean, n_missing.
    - Categorical columns with <= *max_levels* unique values: dtype, n_unique,
      levels (sorted list), n_missing.
    - High-cardinality categorical (> *max_levels*): dtype, n_unique,
      top_levels (top N by frequency), n_missing.

    Parameters
    ----------
    adata : ad.AnnData
        Annotated data object. Modified in place.
    max_levels : int, optional
        Threshold to distinguish low- from high-cardinality categoricals.
    top_n_levels : int, optional
        Number of top levels to keep for high-cardinality categoricals.
    """
    summary: Dict[str, Dict[str, Any]] = {}

    for col in adata.obs.columns:
        series = adata.obs[col]
        n_missing = int(series.isna().sum())

        if pd.api.types.is_numeric_dtype(series):
            summary[col] = {
                "dtype": "numeric",
                "min": float(np.nanmin(series.values)),
                "max": float(np.nanmax(series.values)),
                "mean": float(np.nanmean(series.values)),
                "n_missing": n_missing,
            }
        else:
            # Treat as categorical (includes pd.Categorical, object, string)
            n_unique = int(series.nunique())
            if n_unique <= max_levels:
                levels = sorted(series.dropna().unique().tolist())
                summary[col] = {
                    "dtype": "categorical",
                    "n_unique": n_unique,
                    "levels": levels,
                    "n_missing": n_missing,
                }
            else:
                top_levels = (
                    series.value_counts().head(top_n_levels).index.tolist()
                )
                summary[col] = {
                    "dtype": "categorical",
                    "n_unique": n_unique,
                    "top_levels": top_levels,
                    "n_missing": n_missing,
                }

    adata.uns["column_summary"] = summary
