"""AnnData builder for the EBI Expression Atlas resource.

This module owns ``build_expression_atlas_ad``, the canonical builder that
turns an Expression Atlas baseline bulk-RNA-seq expression TSV plus its SDRF
metadata file into an ``anndata.AnnData`` object (samples x genes). It was
migrated out of :mod:`hvantk.tables.matrix_builders` so that everything
Expression Atlas-specific (builder, SDRF helpers, downloader, dataset class,
tests, fixtures, SKILL) lives under the plugin folder at
:mod:`hvantk.skills.expression_atlas`.

The shared AnnData helpers (``build_anndata_metadata``,
``annotate_column_summary_ad``, ``save_anndata``) intentionally stay in
``hvantk/core/anndata_utils.py`` because they are reused by every anndata
builder (UCSC, CPTAC, ...).
"""

from __future__ import annotations

import logging
from typing import Optional

import anndata as ad

logger = logging.getLogger(__name__)


def build_expression_atlas_ad(
    expression_matrix_path: str,
    sdrf_file: str,
    output_path: Optional[str] = None,
    gene_column: str = "Gene ID",
    gene_name_column: str = "Gene Name",
    delimiter: str = "\t",
    overwrite: bool = False,
) -> "ad.AnnData":
    """Build an AnnData object from Expression Atlas expression + SDRF metadata.

    Parameters
    ----------
    expression_matrix_path : str
        Path to Expression Atlas expression TSV (genes x samples).
    sdrf_file : str
        Path to SDRF metadata file.
    output_path : str, optional
        If provided, save the AnnData as ``.h5ad``.
    gene_column : str
        Name of the gene identifier column (default ``"Gene ID"``).
    gene_name_column : str
        Name of the gene name column (default ``"Gene Name"``).
    delimiter : str
        Column delimiter (default tab).
    overwrite : bool
        Allow overwriting *output_path* if it exists.

    Returns
    -------
    ad.AnnData
        Expression AnnData with SDRF metadata in ``obs`` and provenance in
        ``uns``.
    """
    from hvantk.skills.expression_atlas.shared.expression_atlas import (
        convert_sdrf_to_dataframe,
        create_anndata_from_expression_atlas,
    )
    from hvantk.core.models.anndata_utils import (
        build_anndata_metadata,
        annotate_column_summary_ad,
        save_anndata,
    )

    logger.info("Loading SDRF metadata from %s", sdrf_file)
    metadata_df = convert_sdrf_to_dataframe(sdrf_file)

    logger.info("Creating AnnData from Expression Atlas matrix")
    adata = create_anndata_from_expression_atlas(
        expression_matrix_path=expression_matrix_path,
        metadata_df=metadata_df,
        gene_id_column=gene_column,
        gene_name_column=gene_name_column,
        delimiter=delimiter,
    )

    adata.uns["hvantk_metadata"] = build_anndata_metadata(
        "ExpressionAtlas", expression_matrix_path
    )
    annotate_column_summary_ad(adata)

    if output_path:
        save_anndata(adata, output_path, overwrite=overwrite)

    return adata
