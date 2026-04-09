"""
MatrixTable builders for converting raw sources into Hail MatrixTables (MT).

This module provides unified builders that wrap existing per-source functions
and normalize options and naming.
"""

from __future__ import annotations

import logging
from typing import Optional, Dict, List

logger = logging.getLogger(__name__)

__all__ = [
    "build_ucsc_mt",
    "build_ucsc_ad",
    "build_expression_atlas_mt",
    "build_expression_atlas_ad",
    "build_cptac_mt",
    "build_cptac_ad",
    "build_cptac_phospho_mt",
    "build_cptac_phospho_ad",
]


def build_ucsc_mt(
    expression_matrix_path: str,
    metadata_path: str,
    output_mt: Optional[str] = None,
    gene_column: str = "gene",
    metadata_index_col: int = 0,
    delimiter: str = "\t",
    min_partitions: int = 50,
    force_bgz: bool = True,
    split_gene_field: bool = True,
    overwrite: bool = False,
    auto_convert_bgz: bool = False,
) -> hl.MatrixTable:
    """
    Build a MatrixTable from UCSC Cell Browser expression + metadata files.

    Returns a MatrixTable with columns keyed by 'cell_id' and rows keyed by the
    specified gene_column (default: 'gene'). If output_mt is provided, checkpoint
    is written there.
    """
    import hail as hl
    from hvantk.core.metadata import build_matrix_metadata
    from hvantk.tables.ucsc import (
        convert_ucsc_metadata_to_hail_table,
        create_mt_from_ucsc_expression_matrix,
    )

    logger.info("Converting UCSC metadata to Hail Table")
    metadata_ht = convert_ucsc_metadata_to_hail_table(
        metadata_path=metadata_path,
        sep=delimiter,
        index_col=metadata_index_col,
    )

    logger.info("Creating MatrixTable from UCSC expression matrix")
    mt = create_mt_from_ucsc_expression_matrix(
        expression_matrix_path=expression_matrix_path,
        output_path=None,  # checkpoint after column summary annotation
        delimiter=delimiter,
        row_fields={gene_column: hl.tstr},
        row_key=gene_column,
        split_gene_field=split_gene_field,
        min_partitions=min_partitions,
        force_bgz=force_bgz,
        overwrite=overwrite,
        metadata_ht=metadata_ht,
        auto_convert_bgz=auto_convert_bgz,
    )

    from hvantk.utils.matrix_utils import annotate_column_summary

    mt = annotate_column_summary(mt)
    mt = mt.annotate_globals(
        hvantk_metadata=build_matrix_metadata("UCSC", expression_matrix_path, mt)
    )

    if output_mt:
        mt = mt.checkpoint(output_mt, overwrite=overwrite)

    return mt


def build_ucsc_ad(
    expression_matrix_path: str,
    metadata_path: str,
    output_path: Optional[str] = None,
    gene_column: str = "gene",
    delimiter: str = "\t",
    split_gene_field: bool = True,
    overwrite: bool = False,
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

    Returns
    -------
    ad.AnnData
        Expression AnnData with metadata in ``obs`` and provenance in ``uns``.
    """
    import anndata as ad

    from hvantk.tables.ucsc import load_ucsc_metadata, create_anndata_from_ucsc_matrix
    from hvantk.core.anndata_utils import (
        build_anndata_metadata,
        annotate_column_summary_ad,
        save_anndata,
    )

    logger.info("Loading UCSC metadata from %s", metadata_path)
    metadata_df = load_ucsc_metadata(metadata_path)

    logger.info("Creating AnnData from UCSC expression matrix")
    adata = create_anndata_from_ucsc_matrix(
        expression_matrix_path=expression_matrix_path,
        metadata_df=metadata_df,
        gene_column=gene_column,
        delimiter=delimiter,
        split_gene_field=split_gene_field,
    )

    adata.uns["hvantk_metadata"] = build_anndata_metadata(
        "UCSC", expression_matrix_path
    )
    annotate_column_summary_ad(adata)

    if output_path:
        save_anndata(adata, output_path, overwrite=overwrite)

    return adata


def build_expression_atlas_mt(
    expression_matrix_path: str,
    sdrf_file: str,
    output_mt: Optional[str] = None,
    gene_column: str = "Gene ID",
    sample_id_column: str = "sample_id",
    delimiter: str = "\t",
    min_partitions: int = 50,
    force_bgz: bool = False,
    overwrite: bool = False,
    auto_convert_bgz: bool = False,
) -> hl.MatrixTable:
    """
    Build a MatrixTable from Expression Atlas matrix + SDRF metadata.

    Returns a MatrixTable with columns keyed by 'sample_id' and rows keyed by
    the specified gene_column (default: 'Gene ID'). If output_mt is provided,
    checkpoint is written there.
    """
    import hail as hl
    from hvantk.core.metadata import build_matrix_metadata
    from hvantk.tables.expression_atlas import (
        convert_sdrf_to_hail_table,
        create_mt_from_expression_atlas_matrix,
    )

    logger.info("Converting SDRF metadata to Hail Table")
    metadata_ht = convert_sdrf_to_hail_table(
        sdrf_file=sdrf_file,
        output_file=None,
        keys=[sample_id_column],
        repartition=min_partitions,
        overwrite=overwrite,
    )

    logger.info("Creating MatrixTable from Expression Atlas matrix")
    mt = create_mt_from_expression_atlas_matrix(
        expression_matrix_path=expression_matrix_path,
        output_path=None,  # checkpoint after column summary annotation
        delimiter=delimiter,
        row_fields={gene_column: hl.tstr},
        row_key=gene_column,
        min_partitions=min_partitions,
        force_bgz=force_bgz,
        overwrite=overwrite,
        metadata_ht=metadata_ht,
        auto_convert_bgz=auto_convert_bgz,
    )

    from hvantk.utils.matrix_utils import annotate_column_summary

    mt = annotate_column_summary(mt)
    mt = mt.annotate_globals(
        hvantk_metadata=build_matrix_metadata(
            "ExpressionAtlas", expression_matrix_path, mt
        )
    )

    if output_mt:
        mt = mt.checkpoint(output_mt, overwrite=overwrite)

    return mt


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


def build_cptac_mt(
    expression_path: str,
    metadata_path: str,
    output_mt: Optional[str] = None,
    gene_id_col: str = "GeneID",
    gene_name_col: Optional[str] = "Gene Name",
    sample_id_col: str = "SampleID",
    expression_col: str = "Expression",
    categorical_cols: Optional[List[str]] = None,
    numeric_cols: Optional[List[str]] = None,
    overwrite: bool = False,
) -> hl.MatrixTable:
    """
    Build a CPTAC MatrixTable from expression and metadata TSV/CSV inputs.

    Reads the files into pandas DataFrames, constructs a MatrixTable using
    hvantk.tables.cptac helpers, and optionally checkpoints to output_mt.
    """
    import pandas as pd
    from hvantk.core.metadata import build_matrix_metadata
    from hvantk.tables.cptac import create_cptac_matrix_table

    logger.info("Reading CPTAC expression table from %s", expression_path)
    expr_df = pd.read_csv(
        expression_path,
        sep="\t" if expression_path.endswith((".tsv", ".tsv.bgz", ".tsv.gz")) else ",",
    )

    logger.info("Reading CPTAC metadata table from %s", metadata_path)
    meta_df = pd.read_csv(
        metadata_path,
        sep="\t" if metadata_path.endswith((".tsv", ".tsv.bgz", ".tsv.gz")) else ",",
    )

    logger.info("Building CPTAC MatrixTable")
    mt = create_cptac_matrix_table(
        expression_df=expr_df,
        metadata_df=meta_df,
        gene_id_col=gene_id_col,
        gene_name_col=gene_name_col,
        sample_id_col=sample_id_col,
        expression_col=expression_col,
        categorical_cols=categorical_cols,
        numeric_cols=numeric_cols,
    )

    from hvantk.utils.matrix_utils import annotate_column_summary

    mt = annotate_column_summary(mt)
    mt = mt.annotate_globals(
        hvantk_metadata=build_matrix_metadata("CPTAC", expression_path, mt)
    )

    if output_mt:
        logger.info("Checkpointing CPTAC MatrixTable to %s", output_mt)
        mt = mt.checkpoint(output_mt, overwrite=overwrite)

    return mt


def build_cptac_phospho_mt(
    expression_path: str,
    metadata_path: str,
    output_mt: Optional[str] = None,
    site_id_col: str = "SiteID",
    sample_id_col: str = "SampleID",
    categorical_cols: Optional[List[str]] = None,
    numeric_cols: Optional[List[str]] = None,
    overwrite: bool = False,
) -> "hl.MatrixTable":
    """Build a CPTAC phospho MatrixTable from a sites-by-samples matrix.

    Parameters
    ----------
    expression_path : str
        Path to matrix CSV (sites x samples, from cptac_phospho_datasets).
    metadata_path : str
        Path to metadata CSV (sample clinical info).
    output_mt : str, optional
        Path to checkpoint the MatrixTable.
    site_id_col : str
        Column name for site identifiers (default: "SiteID").
    sample_id_col : str
        Column name for sample identifiers in metadata (default: "SampleID").
    categorical_cols : list of str, optional
        Metadata columns to cast as string.
    numeric_cols : list of str, optional
        Metadata columns to cast as float.
    overwrite : bool
        If True, overwrite existing output.

    Returns
    -------
    hl.MatrixTable
    """
    import hail as hl
    import pandas as pd

    from hvantk.core.metadata import build_matrix_metadata
    from hvantk.utils.matrix_utils import annotate_column_summary

    logger.info("Building CPTAC phospho MatrixTable")

    # Read expression matrix (sites x samples)
    expr_df = pd.read_csv(expression_path, index_col=0)
    logger.info("Expression matrix: %d sites x %d samples", *expr_df.shape)

    # Melt to coordinate format: SiteID, SampleID, Intensity
    expr_long = expr_df.stack(dropna=False).reset_index()
    expr_long.columns = [site_id_col, sample_id_col, "Intensity"]
    expr_long[site_id_col] = expr_long[site_id_col].astype(str)
    expr_long[sample_id_col] = expr_long[sample_id_col].astype(str)

    # Parse site ID into gene_symbol, amino_acid, residue_pos
    # Expected format: Gene_AminoAcidPosition (e.g. TP53_S315, MAPK1_T185)
    import re
    _site_id_re = re.compile(r"^(.+)_([STY])(\d+)$")

    def _parse_site_id(sid):
        m = _site_id_re.match(sid)
        if m:
            return m.group(1), m.group(2), int(m.group(3))
        return sid, "", 0

    site_info = {sid: _parse_site_id(sid) for sid in expr_long[site_id_col].unique()}
    expr_long["gene_symbol"] = expr_long[site_id_col].map(lambda s: site_info[s][0])
    expr_long["amino_acid"] = expr_long[site_id_col].map(lambda s: site_info[s][1])
    expr_long["residue_pos"] = expr_long[site_id_col].map(lambda s: site_info[s][2])

    # Create Hail Table from long-format expression
    ht_expr = hl.Table.from_pandas(expr_long)
    ht_expr = ht_expr.key_by(site_id_col, sample_id_col)

    # Convert to MatrixTable
    mt = ht_expr.to_matrix_table(
        row_key=[site_id_col],
        col_key=[sample_id_col],
        row_fields=["gene_symbol", "amino_acid", "residue_pos"],
    )

    # Read and join metadata
    meta_df = pd.read_csv(metadata_path, index_col=0)
    meta_df.index = meta_df.index.astype(str)
    meta_df.index.name = sample_id_col

    if categorical_cols:
        for c in categorical_cols:
            if c in meta_df.columns:
                meta_df[c] = meta_df[c].astype(str)
    if numeric_cols:
        for c in numeric_cols:
            if c in meta_df.columns:
                meta_df[c] = pd.to_numeric(meta_df[c], errors="coerce")

    ht_meta = hl.Table.from_pandas(meta_df.reset_index())
    ht_meta = ht_meta.key_by(sample_id_col)
    mt = mt.annotate_cols(**ht_meta[mt.col_key])

    mt = annotate_column_summary(mt)
    mt = mt.annotate_globals(
        hvantk_metadata=build_matrix_metadata("CPTAC-phospho", expression_path, mt)
    )

    if output_mt:
        logger.info("Checkpointing CPTAC phospho MatrixTable to %s", output_mt)
        mt = mt.checkpoint(output_mt, overwrite=overwrite)

    return mt


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
    expr_df = pd.read_csv(expression_path, sep="\t")

    logger.info("Reading CPTAC metadata from %s", metadata_path)
    meta_df = pd.read_csv(metadata_path, sep="\t")
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
    expr_df = pd.read_csv(expression_path, sep="\t")

    logger.info("Reading CPTAC phospho metadata from %s", metadata_path)
    meta_df = pd.read_csv(metadata_path, sep="\t")
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
