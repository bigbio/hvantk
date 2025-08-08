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
        for _, row in df.iterrows():
            gene_id = row[gene_id_col]
            gene_name = row[gene_name_col]
            gene_name_dict[gene_id] = gene_name

    # Work directly with the coordinate format - no need to pivot
    # Select only the required columns for the MatrixTable
    coord_df = df[[gene_id_col, sample_id_col, expression_col]].copy()
    
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
        ValueError: If sample_id_col is missing from the input DataFrame
    """
    logger.info("Converting CPTAC metadata to Table")
    
    if sample_id_col not in metadata_df.columns:
        raise ValueError(f"Sample ID column '{sample_id_col}' not found in metadata")
    
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
    expr_samples = set(mt.col_key.collect())
    meta_samples = set(metadata_ht[sample_id_col].collect())

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
