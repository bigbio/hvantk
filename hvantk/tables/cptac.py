"""
Module for processing CPTAC (Clinical Proteomic Tumor Analysis Consortium) data.

This module provides functions to convert CPTAC expression and metadata data
into Hail MatrixTables and Tables for downstream analysis.
"""

import logging
from typing import Optional, List

import hail as hl
import pandas as pd

logger = logging.getLogger(__name__)

__all__ = [
    "convert_cptac_expression_to_matrix_table",
    "convert_cptac_metadata_to_table",
    "create_cptac_matrix_table",
    "save_cptac_matrix_table",
]


def convert_cptac_expression_to_matrix_table(
    expression_df: pd.DataFrame,
    gene_id_col: str = 'GeneID',
    gene_name_col: Optional[str] = 'Gene Name',
    sample_id_col: str = 'SampleID',
    expression_col: str = 'Expression'
) -> hl.MatrixTable:
    """
    Converts a CPTAC expression dataframe to a Hail MatrixTable.

    Args:
        expression_df: DataFrame containing gene expression data
        gene_id_col: Column name containing gene IDs
        gene_name_col: Column name containing gene names (optional)
        sample_id_col: Column name containing sample IDs
        expression_col: Column name containing expression values

    Returns:
        A Hail MatrixTable with genes as rows and samples as columns

    Raises:
        ValueError: If required columns are missing from the input DataFrame
    """
    logger.info("Converting CPTAC expression data to MatrixTable")

    # Validate required columns
    required_cols = [gene_id_col, sample_id_col, expression_col]
    missing_cols = [col for col in required_cols if col not in expression_df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    # Create a copy of the DataFrame to avoid modifying the original
    df = expression_df.copy()

    # First, we'll create a dictionary mapping gene IDs to gene names
    # This way we can add gene_name as a row annotation after creating the MatrixTable
    gene_name_dict = {}
    if gene_name_col and gene_name_col in df.columns:
        # Create a dictionary mapping gene IDs to their names
        gene_name_dict = df.drop_duplicates(subset=[gene_id_col]).set_index(gene_id_col)[gene_name_col].to_dict()

    # Work directly with the coordinate format - no need to pivot
    # Select only the required columns for the MatrixTable
    coord_df = df[[gene_id_col, sample_id_col, expression_col]].copy()

    # Ensure key columns are strings for consistent joins
    coord_df[gene_id_col] = coord_df[gene_id_col].astype(str)
    coord_df[sample_id_col] = coord_df[sample_id_col].astype(str)

    # Convert to Hail Table
    coord_ht = hl.Table.from_pandas(coord_df)

    # Convert to MatrixTable using coordinate representation
    mt = coord_ht.to_matrix_table(
        row_key=[gene_id_col],
        col_key=[sample_id_col],
        row_fields=[],
        col_fields=[]
    )

    # If we have gene names, add them as row annotations
    if gene_name_dict:
        # Create a literal dictionary in Hail
        gene_name_dict_expr = hl.literal(gene_name_dict)
        # Use the dictionary to annotate rows with gene names
        mt = mt.annotate_rows(gene_name=gene_name_dict_expr.get(mt[gene_id_col]))

    logger.info(f"Created MatrixTable with {mt.count_rows()} genes and {mt.count_cols()} samples")
    return mt


def convert_cptac_metadata_to_table(
    metadata_df: pd.DataFrame,
    sample_id_col: str = 'SampleID',
    categorical_cols: Optional[List[str]] = None,
    numeric_cols: Optional[List[str]] = None
) -> hl.Table:
    """
    Converts CPTAC metadata to a Hail Table.

    Args:
        metadata_df: DataFrame containing sample metadata
        sample_id_col: Column name containing sample IDs
        categorical_cols: List of column names to be treated as categorical variables
        numeric_cols: List of column names to be treated as numeric variables

    Returns:
        A Hail Table with samples as rows and metadata as columns

    Raises:
        ValueError: If sample_id_col is missing from the input DataFrame or duplicate sample IDs exist
    """
    logger.info("Converting CPTAC metadata to Table")

    # Check if sample_id_col exists in metadata_df
    if sample_id_col not in metadata_df.columns:
        raise ValueError(f"Sample ID column '{sample_id_col}' not found in metadata")
    metadata_df = metadata_df.copy()
    metadata_df[sample_id_col] = metadata_df[sample_id_col].astype(str)

    # Detect duplicate sample IDs (fast fail to avoid incorrect joins / explode)
    if metadata_df[sample_id_col].duplicated().any():
        dup_counts = metadata_df[sample_id_col][metadata_df[sample_id_col].duplicated(keep=False)].value_counts()
        # Limit list length in message if extremely large
        duplicate_list = dup_counts.index.tolist()
        if len(duplicate_list) > 50:
            shown = duplicate_list[:50]
            more = len(duplicate_list) - 50
            display_ids = f"{shown} (+{more} more)"
        else:
            display_ids = str(duplicate_list)
        raise ValueError(
            "Duplicate sample IDs found in metadata ({} duplicates across {} unique IDs). Example duplicates: {}. Counts: {}".format(
                dup_counts.sum() - dup_counts.shape[0],  # total duplicate entries beyond first occurrences
                len(dup_counts),
                display_ids,
                dup_counts.to_dict()
            )
        )

    # Convert to Hail Table
    ht = hl.Table.from_pandas(metadata_df)

    # Set sample ID as key
    ht = ht.key_by(sample_id_col)

    # Convert specified columns to appropriate types
    if categorical_cols:
        for col in categorical_cols:
            if col in ht.row:
                ht = ht.annotate(**{col: hl.str(ht[col])})

    if numeric_cols:
        for col in numeric_cols:
            if col in ht.row:
                ht = ht.annotate(**{col: hl.float64(ht[col])})

    logger.info(f"Created metadata Table with {ht.count()} samples")
    return ht


def create_cptac_matrix_table(
    expression_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    gene_id_col: str = 'GeneID',
    gene_name_col: Optional[str] = 'Gene Name',
    sample_id_col: str = 'SampleID',
    expression_col: str = 'Expression',
    categorical_cols: Optional[List[str]] = None,
    numeric_cols: Optional[List[str]] = None
) -> hl.MatrixTable:
    """
    Creates a complete CPTAC MatrixTable with expression data and metadata.

    Args:
        expression_df: DataFrame containing gene expression data
        metadata_df: DataFrame containing sample metadata
        gene_id_col: Column name containing gene IDs
        gene_name_col: Column name containing gene names (optional)
        sample_id_col: Column name containing sample IDs
        expression_col: Column name containing expression values
        categorical_cols: List of metadata columns to be treated as categorical
        numeric_cols: List of metadata columns to be treated as numeric

    Returns:
        A Hail MatrixTable with genes as rows, samples as columns, and metadata annotated

    Raises:
        ValueError: If there are mismatches between expression and metadata sample IDs
    """
    logger.info("Creating complete CPTAC MatrixTable")

    # Create expression MatrixTable
    mt = convert_cptac_expression_to_matrix_table(
        expression_df,
        gene_id_col=gene_id_col,
        gene_name_col=gene_name_col,
        sample_id_col=sample_id_col,
        expression_col=expression_col
    )

    # Create metadata Table
    metadata_ht = convert_cptac_metadata_to_table(
        metadata_df,
        sample_id_col=sample_id_col,
        categorical_cols=categorical_cols,
        numeric_cols=numeric_cols
    )

    # Check for sample ID mismatches - get the sample IDs from the MatrixTable columns
    # Use mt.col_key to access the column keys (sample IDs) properly
    # Extract the sample ID values from the struct since col_key returns Struct objects
    meta_samples = set(metadata_ht[sample_id_col].collect())

    expr_samples = set(mt.aggregate_cols(hl.agg.collect(mt.col_key[sample_id_col])))

    if expr_samples != meta_samples:
        missing_in_meta = expr_samples - meta_samples
        missing_in_expr = meta_samples - expr_samples
        error_messages = []
        if missing_in_meta:
            error_messages.append(f"Samples in expression data but not in metadata: {missing_in_meta}")
        if missing_in_expr:
            error_messages.append(f"Samples in metadata but not in expression data: {missing_in_expr}")
        raise ValueError("Sample ID mismatches found. " + "; ".join(error_messages))

    # Annotate MatrixTable with metadata
    # Use the sample_id_col to join with metadata
    mt = mt.annotate_cols(**metadata_ht[mt.col_key])

    logger.info("Successfully created annotated MatrixTable")
    return mt


def save_cptac_matrix_table(
    mt: hl.MatrixTable,
    output_path: str,
    overwrite: bool = False
) -> None:
    """
    Saves a CPTAC MatrixTable to disk.

    Args:
        mt: The MatrixTable to save
        output_path: Path where the MatrixTable will be saved
        overwrite: Whether to overwrite existing files

    Raises:
        ValueError: If the output path already exists and overwrite is False
    """
    logger.info(f"Saving MatrixTable to {output_path}")

    try:
        mt.write(output_path, overwrite=overwrite)
        logger.info("Successfully saved MatrixTable")
    except Exception as e:
        logger.error(f"Failed to save MatrixTable: {str(e)}")
        raise
