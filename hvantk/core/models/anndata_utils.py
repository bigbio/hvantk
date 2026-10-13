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


class TranscriptKeyedError(ValueError):
    """Raised by :func:`require_gene_keyed` for a transcript-keyed AnnData.

    A distinct type (rather than a bare ``ValueError``) so a CLI can catch
    exactly this usage error and report it cleanly, without also swallowing an
    unrelated ``ValueError`` raised deeper in the same call (e.g. a numpy
    broadcasting failure) as if it were the same kind of problem.
    """


def require_gene_keyed(adata: ad.AnnData, consumer: str) -> None:
    """Raise if ``adata.var`` is keyed by transcript id rather than gene id.

    A transcript-level Expression Atlas build indexes ``var`` by the transcript
    id (``create_anndata_from_expression_atlas``, issue #349) and keeps the
    gene id as a regular ``"Gene ID"`` column instead -- ``var.index.name ==
    "transcript_id"``. Generic AnnData consumers that read ``var_names``/
    ``adata.var.index`` as the gene identifier (marker extraction, per-group
    summarization, the PTM constraint expression loader, the annotation
    matrix reducer) would otherwise silently emit one row per transcript,
    labelled as a gene (issue #383).

    Call this at the top of any such consumer's AnnData entry point, before it
    treats ``var_names``/``adata.var.index`` as gene identifiers.

    Parameters
    ----------
    adata : ad.AnnData
        The AnnData to check. Not modified.
    consumer : str
        Human-readable name of the calling command/function, used only to
        make the error message actionable (e.g. ``"hvantk expression
        markers"``).

    Raises
    ------
    TranscriptKeyedError
        If ``adata.var.index.name == "transcript_id"``.
    """
    if adata.var.index.name == "transcript_id":
        raise TranscriptKeyedError(
            f"{consumer}: this AnnData is keyed by transcript id "
            "(var.index.name == 'transcript_id'), not gene id, so its var_names "
            "cannot be read as genes here. Rebuild it with `--plugin-arg "
            "level=gene` (e.g. `hvantk reprocess expression-atlas:dataset "
            "--raw-dir <dir> --output <out>.h5ad --plugin-arg "
            "accession=<ACCESSION> --plugin-arg level=gene`, or `--skip-download` "
            "if <dir> already holds the raw files) to get a gene-keyed matrix, "
            "or use a gene-level export."
        )


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
                top_levels = series.value_counts().head(top_n_levels).index.tolist()
                summary[col] = {
                    "dtype": "categorical",
                    "n_unique": n_unique,
                    "top_levels": top_levels,
                    "n_missing": n_missing,
                }

    adata.uns["column_summary"] = summary
