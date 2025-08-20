"""
Utility functions to convert UCSC Cell Browser datasets to Hail Matrix Tables.

This module provides functions to convert UCSC metadata and expression matrix files
into Hail Tables and MatrixTables for downstream genetic and single-cell analysis.
"""

import os
import hail as hl
import pandas as pd

from hvantk.utils.expressions import split_field_expr
from hvantk.core.constants import UCSC_CELL_ID_COLUMN, UCSC_GENE_COLUMN

__all__ = [
    "convert_ucsc_metadata_to_hail_table",
    "create_mt_from_ucsc_expression_matrix",
]


def _replace_dots_in_column_names(df: pd.DataFrame) -> pd.DataFrame:
    """
    Replaces dots in the column names of a DataFrame with underscores.

    :param df: Input DataFrame with column names containing dots.
    :type df: pd.DataFrame
    :return: DataFrame with dots replaced by underscores in column names.
    :rtype: pd.DataFrame
    """
    df.columns = df.columns.str.replace(".", "_", regex=False)
    return df


def convert_ucsc_metadata_to_hail_table(
    metadata_path: str,
    sep: str = "\t",
    index_col: int = 0,
    index_name: str = UCSC_CELL_ID_COLUMN,
) -> hl.Table:
    """
    Converts a tabular metadata file from UCSC expression matrix into a Hail Table.

    :param metadata_path: Path to the metadata file to be converted.
    :type metadata_path: str
    :param sep: Delimiter character used to separate columns in the metadata file.
    :type sep: str
    :param index_col: Column index to be used as the index of the DataFrame.
    :type index_col: int
    :param index_name: Name of the index column when added as a DataFrame column.
    :type index_name: str
    :return: Converted Hail Table with the specified key.
    :rtype: hail.Table
    :raises FileNotFoundError: If the specified metadata file does not exist.
    """

    # Check if the metadata file exists
    if not os.path.exists(metadata_path):
        raise FileNotFoundError(f"Metadata file not found: {metadata_path}")

    # Load the metadata file into a Pandas DataFrame
    df = pd.read_csv(metadata_path, sep=sep, index_col=index_col)

    # Annotate the index as a column
    df.reset_index(inplace=True)

    # Rename the index column
    df.rename(columns={"index": index_name}, inplace=True)

    # Replace dots in column names with underscores
    df = _replace_dots_in_column_names(df)

    # Convert Pandas DataFrame to Hail Table with the specified key
    ht = hl.Table.from_pandas(df, key=index_name)

    return ht


def create_mt_from_ucsc_expression_matrix(
    expression_matrix_path: str,
    output_path: str = None,
    delimiter: str = "\t",
    row_fields: dict = None,
    row_key: str = UCSC_GENE_COLUMN,
    split_gene_field: bool = True,
    min_partitions: int = 50,
    force_bgz: bool = True,
    overwrite: bool = True,
    metadata_ht: hl.Table = None,
) -> hl.MatrixTable:
    """
    Creates a Hail MatrixTable from a given UCSC expression matrix file.
    Annotates the MatrixTable with metadata if provided.

    :param expression_matrix_path: Path to the UCSC expression matrix file.
    :param output_path: Optional path to save the resulting MatrixTable as a checkpoint.
    :param delimiter: Delimiter used in the expression matrix file.
    :param row_fields: Dictionary defining row fields and their types for the matrix. Mapping str to hl.type.
    :param row_key: Key for rows, which must match a field in row_fields.
    :param split_gene_field: Whether to split the gene field and use the first element.
    :param min_partitions: Minimum number of partitions for the imported MatrixTable.
    :param force_bgz: Whether to force BGZF compression for the input file.
    :param overwrite: Whether to allow overwriting existing files at the output path.
    :param metadata_ht: Optional Hail Table containing metadata for annotating columns. Keys must match.
    :return: The Hail MatrixTable generated from the expression matrix file.
    :rtype: Hail MatrixTable
    :raises FileNotFoundError: If the specified expression matrix file does not exist.
    :raises FileExistsError: If the output path exists and overwrite is set to False.
    :raises ValueError: If metadata_ht is provided but not a Hail Table,
        or if its key does not match the matrix table's column key.
    """

    if not os.path.exists(expression_matrix_path):
        raise FileNotFoundError(
            f"Expression matrix file not found: {expression_matrix_path}"
        )

    # Set default row_fields if None
    if row_fields is None:
        row_fields = {UCSC_GENE_COLUMN: hl.tstr}


    # Check if the output path exists and overwrite is set to False
    if output_path and os.path.exists(output_path) and not overwrite:
        raise FileExistsError(
            f"Output path already exists: {output_path}. Set overwrite=True to overwrite."
        )

    # Import the matrix table with the specified parameters
    mt = hl.import_matrix_table(
        expression_matrix_path,
        delimiter=delimiter,
        row_fields=row_fields,
        row_key=row_key,
        min_partitions=min_partitions,
        force_bgz=force_bgz,
    )

    # rename col_id to cell_id
    mt = mt.rename({"col_id": UCSC_CELL_ID_COLUMN})

    # Optionally split the gene field if specified (A|B -> A)
    if split_gene_field and row_key in mt.row:
        # Unkey the matrix table to split the gene field
        mt = mt.key_rows_by()
        # Split the <gene> field and re-annotate the matrix table
        mt = mt.annotate_rows(**{row_key: split_field_expr(mt, field_name=row_key)})
        # Re-key the matrix table with the split gene field
        mt = mt.key_rows_by(mt[row_key])

    # Optionally add metadata to the matrix table if provided
    if metadata_ht is not None:
        # Check if the metadata table is a Hail Table
        if not isinstance(metadata_ht, hl.Table):
            raise ValueError("metadata_ht must be a Hail Table.")

        # Annotate the matrix table with metadata
        # Ensure the metadata table has the same key as the matrix table
        if str(metadata_ht.key) != str(mt.col_key):
            raise ValueError("metadata_ht must have the same key as the matrix table.")

        mt = mt.annotate_cols(metadata=metadata_ht[mt.col_key])

    # If a checkpoint path is specified, save the matrix table with a checkpoint
    if output_path:
        mt = mt.checkpoint(output_path, overwrite=overwrite)

    return mt
