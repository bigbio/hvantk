"""
AnnData builders for converting raw expression sources into AnnData (.h5ad) objects.

This module provides unified builders that wrap existing per-source functions
and normalize options and naming.
"""

from __future__ import annotations

import anndata as ad
import logging
import os
from typing import Optional

logger = logging.getLogger(__name__)

# Auto-select the backed builder for UCSC inputs larger than this threshold.
BACKED_BUILDER_THRESHOLD_BYTES = 1 * 1024 * 1024 * 1024  # 1 GiB

__all__ = [
    "build_ucsc_ad",
    "build_expression_atlas_ad",
    "build_cptac_ad",
    "build_cptac_phospho_ad",
]


def build_ucsc_ad(
    expression_matrix_path: str,
    metadata_path: str,
    output_path: Optional[str] = None,
    gene_column: str = "gene",
    delimiter: str = "\t",
    split_gene_field: bool = True,
    overwrite: bool = False,
    chunk_size: int = 500,
    backed: bool | None = None,
    column_batch: int = 64,
) -> "ad.AnnData":
    """Build an AnnData object from UCSC Cell Browser expression + metadata.

    Parameters
    ----------
    expression_matrix_path : str
        Path to expression TSV (genes x cells).
    metadata_path : str
        Path to metadata TSV.
    output_path : str, optional
        If provided, save the AnnData as ``.h5ad``.
    gene_column : str
        Name of the gene identifier column (default ``"gene"``).
    delimiter : str
        Column delimiter (default tab).
    split_gene_field : bool
        Split pipe-separated gene names, keeping the first element.
    overwrite : bool
        Allow overwriting *output_path* if it exists.
    chunk_size : int
        Gene rows per streaming chunk (default 500).
    backed : bool, optional
        Force the backed-write builder. When ``None`` (default), auto-select
        based on the expression matrix file size (``> BACKED_BUILDER_THRESHOLD_BYTES``
        triggers backed mode). When ``True``, ``output_path`` is required.
    column_batch : int
        Cell-column batch size used by the backed builder (default 64).

    Returns
    -------
    ad.AnnData
        Expression AnnData with metadata in ``obs`` and provenance in ``uns``.
        In backed mode, returns a read-backed AnnData handle (``backed='r'``).
    """
    import anndata as ad

    from hvantk.tables.ucsc import (
        load_ucsc_metadata,
        create_anndata_from_ucsc_matrix,
        build_ucsc_atlas_backed,
    )
    from hvantk.core.anndata_utils import (
        build_anndata_metadata,
        annotate_column_summary_ad,
        save_anndata,
    )

    logger.info("Loading UCSC metadata from %s", metadata_path)
    metadata_df = load_ucsc_metadata(metadata_path, sep=delimiter)

    # Auto-select backed mode by input size unless caller forces.
    if backed is None:
        size = os.path.getsize(expression_matrix_path)
        backed = size > BACKED_BUILDER_THRESHOLD_BYTES
        logger.info(
            "build_ucsc_ad: auto-selected backed=%s (input size %.2f GiB, threshold %.2f GiB)",
            backed, size / (1024**3), BACKED_BUILDER_THRESHOLD_BYTES / (1024**3),
        )

    if backed:
        if output_path is None:
            raise ValueError("backed=True requires output_path (writes directly to disk).")
        provenance = {"hvantk_metadata": build_anndata_metadata("UCSC", expression_matrix_path)}
        build_ucsc_atlas_backed(
            expression_matrix_path=expression_matrix_path,
            output_path=output_path,
            metadata_df=metadata_df,
            gene_column=gene_column,
            delimiter=delimiter,
            split_gene_field=split_gene_field,
            column_batch=column_batch,
            overwrite=overwrite,
            uns=provenance,
        )
        # Return a backed-mode handle — shape + obs/var without materializing X.
        # annotate_column_summary_ad would need to scan the full X matrix and
        # is intentionally skipped for backed atlases (v1 trade-off).
        logger.info(
            "Backed atlas built at %s; skipping annotate_column_summary_ad "
            "(would materialize X). Returning a backed AnnData handle.",
            output_path,
        )
        return ad.read_h5ad(output_path, backed="r")

    logger.info("Creating AnnData from UCSC expression matrix (in-memory)")
    adata = create_anndata_from_ucsc_matrix(
        expression_matrix_path=expression_matrix_path,
        metadata_df=metadata_df,
        gene_column=gene_column,
        delimiter=delimiter,
        split_gene_field=split_gene_field,
        chunk_size=chunk_size,
    )
    adata.uns["hvantk_metadata"] = build_anndata_metadata(
        "UCSC", expression_matrix_path
    )
    annotate_column_summary_ad(adata)
    if output_path:
        save_anndata(adata, output_path, overwrite=overwrite)
    return adata


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
    from hvantk.tables.expression_atlas import (
        convert_sdrf_to_dataframe,
        create_anndata_from_expression_atlas,
    )
    from hvantk.core.anndata_utils import (
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

    from hvantk.tables.cptac import create_anndata_from_cptac_long
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


def build_cptac_phospho_ad(
    expression_path: str,
    metadata_path: str,
    output_path: Optional[str] = None,
    site_id_col: str = "SiteID",
    sample_id_col: str = "SampleID",
    overwrite: bool = False,
) -> "ad.AnnData":
    """Build an AnnData from CPTAC phosphoproteomics wide-format matrix + metadata.

    Parameters
    ----------
    expression_path : str
        Path to wide-format phospho expression TSV/CSV (sites x samples).
    metadata_path : str
        Path to sample metadata TSV/CSV.
    output_path : str, optional
        If provided, save the AnnData as ``.h5ad``.
    site_id_col : str
        Column containing site identifiers (e.g. ``TP53_S315``).
    sample_id_col : str
        Column containing sample identifiers in metadata.
    overwrite : bool
        Allow overwriting *output_path* if it exists.

    Returns
    -------
    ad.AnnData
        Expression AnnData (samples x sites) with parsed site annotations
        in ``var`` and metadata in ``obs``.
    """
    import pandas as pd

    from hvantk.tables.cptac import create_anndata_from_cptac_phospho
    from hvantk.core.anndata_utils import (
        build_anndata_metadata,
        annotate_column_summary_ad,
        save_anndata,
    )

    logger.info("Reading CPTAC phospho expression from %s", expression_path)
    expr_df = pd.read_csv(expression_path, sep=None, engine="python")

    logger.info("Reading CPTAC phospho metadata from %s", metadata_path)
    meta_df = pd.read_csv(metadata_path, sep=None, engine="python")
    meta_df = meta_df.set_index(sample_id_col)

    logger.info("Building CPTAC phospho AnnData")
    adata = create_anndata_from_cptac_phospho(
        expression_df=expr_df,
        metadata_df=meta_df,
        site_id_col=site_id_col,
        sample_id_col=sample_id_col,
    )

    adata.uns["hvantk_metadata"] = build_anndata_metadata("CPTAC", expression_path)
    annotate_column_summary_ad(adata)

    if output_path:
        save_anndata(adata, output_path, overwrite=overwrite)

    return adata
