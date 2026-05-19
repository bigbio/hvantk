"""AnnData utility functions for provenance, column summaries, and I/O."""

import logging
import os
from datetime import datetime
from typing import Any, Dict

import anndata as ad
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

_ANNDATA_SOURCE_DESCRIPTIONS: Dict[str, str] = {
    "ucsc": "Single-cell expression data from UCSC Cell Browser",
    "expressionatlas": "Bulk/single-cell RNA-seq from EMBL-EBI Expression Atlas",
    "cptac": "Proteomics expression data from the Clinical Proteomic Tumor Analysis Consortium",
}


def _normalize_source_name(source_name: str) -> str:
    """Normalize source name by lowercasing and removing non-alphanumeric chars."""
    return "".join(ch for ch in source_name.lower() if ch.isalnum())


def _get_hvantk_version() -> str:
    """Get hvantk version from package metadata."""
    try:
        from importlib.metadata import version

        return version("hvantk")
    except Exception:
        return "unknown"


def build_anndata_metadata(source_name: str, input_path: str) -> Dict[str, Any]:
    """Build provenance metadata dict for ``adata.uns["hvantk_metadata"]``.

    Parameters
    ----------
    source_name : str
        Human-readable name of the data source.
    input_path : str
        Path to the raw input file.

    Returns
    -------
    dict
        Metadata dict with keys: hvantk_version, source_name,
        source_description, raw_input_path, build_date.
    """
    normalized = _normalize_source_name(source_name)
    return {
        "hvantk_version": _get_hvantk_version(),
        "source_name": source_name,
        "source_description": _ANNDATA_SOURCE_DESCRIPTIONS.get(normalized, ""),
        "raw_input_path": input_path,
        "build_date": datetime.now().isoformat(),
    }


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


def save_anndata(
    adata: ad.AnnData,
    path: str,
    overwrite: bool = True,
) -> None:
    """Save AnnData to ``.h5ad`` file.

    Parameters
    ----------
    adata : ad.AnnData
        Annotated data object to save.
    path : str
        Output file path (should end in ``.h5ad``).
    overwrite : bool, optional
        If False and *path* already exists, raise :class:`FileExistsError`.

    Raises
    ------
    FileExistsError
        If *overwrite* is False and *path* exists.
    """
    if not overwrite and os.path.exists(path):
        raise FileExistsError(f"File already exists: {path}")

    logger.info("Saving AnnData (%d obs x %d var) to %s", adata.n_obs, adata.n_vars, path)
    adata.write_h5ad(path)


def load_anndata(path: str) -> ad.AnnData:
    """Load AnnData from ``.h5ad`` file.

    Parameters
    ----------
    path : str
        Path to the ``.h5ad`` file.

    Returns
    -------
    ad.AnnData
        The loaded annotated data object.
    """
    logger.info("Loading AnnData from %s", path)
    return ad.read_h5ad(path)
