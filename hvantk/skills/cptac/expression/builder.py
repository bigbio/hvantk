"""AnnData builder for the CPTAC long-format protein expression resource.

This module owns ``build_cptac_ad``, the canonical builder that turns a
long-format CPTAC expression TSV plus its sample metadata into an
``anndata.AnnData`` object (samples x genes). It was migrated out of
:mod:`hvantk.tables.matrix_builders` so that everything CPTAC-specific
(builders, helpers, downloader, dataset class, tests, fixtures, SKILLs)
lives under the plugin folder at :mod:`hvantk.skills.cptac`.

The shared AnnData helpers (``build_anndata_metadata``,
``annotate_column_summary_ad``, ``save_anndata``) intentionally stay in
``hvantk/core/anndata_utils.py`` because they are reused by every anndata
builder.
"""

from __future__ import annotations

import logging
from typing import Optional

import anndata as ad

logger = logging.getLogger(__name__)


def build_cptac_ad(
    expression_path: str,
    metadata_path: str,
    output_path: Optional[str] = None,
    gene_id_col: str = "GeneID",
    gene_name_col: Optional[str] = "Gene Name",
    sample_id_col: str = "SampleID",
    expression_col: str = "Expression",
    overwrite: bool = False,
) -> "ad.AnnData":
    """Build an AnnData object from long-format CPTAC expression + metadata.

    Parameters
    ----------
    expression_path : str
        Path to long-format expression TSV/CSV.
    metadata_path : str
        Path to sample metadata TSV/CSV.
    output_path : str, optional
        If provided, save the AnnData as ``.h5ad``.
    gene_id_col : str
        Column containing gene identifiers.
    gene_name_col : str, optional
        Column containing gene names.
    sample_id_col : str
        Column containing sample identifiers.
    expression_col : str
        Column containing expression values.
    overwrite : bool
        Allow overwriting *output_path* if it exists.

    Returns
    -------
    ad.AnnData
        Expression AnnData with metadata in ``obs`` and provenance in ``uns``.
    """
    import pandas as pd

    from hvantk.skills.cptac.shared.cptac import create_anndata_from_cptac_long
    from hvantk.core.anndata_utils import (
        build_anndata_metadata,
        annotate_column_summary_ad,
        save_anndata,
    )

    logger.info("Reading CPTAC expression from %s", expression_path)
    expr_df = pd.read_csv(expression_path, sep=None, engine="python")

    logger.info("Reading CPTAC metadata from %s", metadata_path)
    meta_df = pd.read_csv(metadata_path, sep=None, engine="python")
    meta_df = meta_df.set_index(sample_id_col)

    logger.info("Building CPTAC AnnData")
    adata = create_anndata_from_cptac_long(
        expression_df=expr_df,
        metadata_df=meta_df,
        gene_id_col=gene_id_col,
        gene_name_col=gene_name_col,
        sample_id_col=sample_id_col,
        expression_col=expression_col,
    )

    adata.uns["hvantk_metadata"] = build_anndata_metadata("CPTAC", expression_path)
    annotate_column_summary_ad(adata)

    if output_path:
        save_anndata(adata, output_path, overwrite=overwrite)

    return adata
