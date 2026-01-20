"""
MatrixTable builders for converting raw sources into Hail MatrixTables (MT).

This module provides unified builders that wrap existing per-source functions
and normalize options and naming.
"""

from __future__ import annotations

import logging
from typing import Optional, Dict, List

import hail as hl

logger = logging.getLogger(__name__)

__all__ = [
    "build_ucsc_mt",
    "build_expression_atlas_mt",
    "build_cptac_mt",
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
) -> hl.MatrixTable:
    """
    Build a MatrixTable from UCSC Cell Browser expression + metadata files.

    Returns a MatrixTable with columns keyed by 'cell_id' and rows keyed by the
    specified gene_column (default: 'gene'). If output_mt is provided, checkpoint
    is written there.
    """
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
        output_path=output_mt,
        delimiter=delimiter,
        row_fields={gene_column: hl.tstr},
        row_key=gene_column,
        split_gene_field=split_gene_field,
        min_partitions=min_partitions,
        force_bgz=force_bgz,
        overwrite=overwrite,
        metadata_ht=metadata_ht,
    )
    return mt


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
) -> hl.MatrixTable:
    """
    Build a MatrixTable from Expression Atlas matrix + SDRF metadata.

    Returns a MatrixTable with columns keyed by 'sample_id' and rows keyed by
    the specified gene_column (default: 'Gene ID'). If output_mt is provided,
    checkpoint is written there.
    """
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
        output_path=output_mt,
        delimiter=delimiter,
        row_fields={gene_column: hl.tstr},
        row_key=gene_column,
        min_partitions=min_partitions,
        force_bgz=force_bgz,
        overwrite=overwrite,
        metadata_ht=metadata_ht,
    )
    return mt


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

    if output_mt:
        logger.info("Checkpointing CPTAC MatrixTable to %s", output_mt)
        mt = mt.checkpoint(output_mt, overwrite=overwrite)

    return mt
