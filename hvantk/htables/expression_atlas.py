import pandas as pd
import hail as hl
import os


def _import_sdrf(sdrf_file: str, **kwargs) -> pd.DataFrame:
    """
    Import SDRF (Sample and Data Relationship Format) file and return as pandas DataFrame.
    SDRF files are tab-separated, no header, and contain metadata about samples and their characteristics.
    The second column is usually NaN and can be ignored.

    Parameters:
        sdrf_file (str): Path to the SDRF file
        **kwargs: Additional arguments to pass to pandas.read_csv. Default parameters can be overridden:
                 - sep: "\t"
                 - comment: "#"
                 - header: None
                 - column names: ['accession', 'unused', 'sample_id', 'column_type', 'column_name', 'column_value']

    Returns:
        pd.DataFrame: Processed SDRF data

    Raises:
        FileNotFoundError: If the SDRF file does not exist
        pd.errors.EmptyDataError: If the SDRF file is empty
    """

    try:
        # Define default parameters as a configuration dict
        default_params = {
            "sep": "\t",
            "comment": "#",
            "header": None,
            "names": ['accession', 'unused', 'sample_id', 'column_type', 'column_name', 'column_value']
        }

        # Update defaults with any provided kwargs
        read_params = {**default_params, **kwargs}

        # Infer number of columns if names not provided in kwargs
        if 'names' not in kwargs:
            with open(sdrf_file) as f:
                first_line = f.readline().strip()
                num_cols = len(first_line.split(read_params['sep']))
                read_params['usecols'] = range(min(6, num_cols))

        # Read file
        sdrf_df = pd.read_csv(sdrf_file, **read_params)

        # Clean column names
        sdrf_df.columns = sdrf_df.columns.str.strip()

        # Drop unused columns
        sdrf_df.drop(columns=['unused'], inplace=True, errors='ignore')

        return sdrf_df

    except FileNotFoundError as e:
        raise FileNotFoundError(f"SDRF file not found: {sdrf_file}") from e
    except pd.errors.EmptyDataError as e:
        raise pd.errors.EmptyDataError(f"SDRF file is empty: {sdrf_file}") from e


def _reshape_sdrf_long_to_wide_format(df_sdrf: pd.DataFrame,
                                          include_only_factors: bool = False,
                                          include_only_characteristic: bool = False) -> pd.DataFrame:
        """
        Reshapes an input DataFrame from long format to wide format based on specific
        filtering and pivoting rules.

        Args:
            df_sdrf (pd.DataFrame): Input DataFrame containing columns: 'sample_name',
                'column_type', 'column_name', and 'column_value'.
            include_only_factors (bool, optional): If True, include only factor rows. Defaults to False.
            include_only_characteristic (bool, optional): If True, include only characteristic rows.
                Defaults to False.

        Returns:
            pd.DataFrame: Reshaped DataFrame in wide format with:
                - Columns representing distinct column_name values
                - Rows representing sample_name entries
                - Values from column_value
                - Clean column names (spaces replaced with underscores)

        Raises:
            ValueError: If both include_only_factors and include_only_characteristic are True
        """
        if include_only_factors and include_only_characteristic:
            raise ValueError("Cannot set both include_only_factors and include_only_characteristic to True")

        df = df_sdrf.copy()

        # Filter rows based on column_type
        if include_only_characteristic:
            df = df[df['column_type'] == 'characteristic']
        elif include_only_factors:
            df = df[df['column_type'] == 'factor']
        # else: use all column types

        # Handle duplicate sample_id/column_name combinations by taking last value
        pivot_df = df[['sample_id', 'column_name', 'column_value']].drop_duplicates(
            subset=['sample_id', 'column_name'],
            keep='last'
        )

        # Reshape from long to wide format
        wide_df = pivot_df.pivot(
            index='sample_id',
            columns='column_name',
            values='column_value'
        )

        # Reset index to make sample_id a column
        wide_df.reset_index(inplace=True)

        # Clean column names (replace spaces and special chars with underscores)
        wide_df.columns = [str(col).strip().replace(' ', '_').replace('(', '').replace(')', '')
                           for col in wide_df.columns]

        return wide_df


def convert_sdrf_to_hail_table(
        sdrf_file: str,
        output_file: str=None,
        keys=None,
        repartition: int = 50,
        overwrite: bool = False,
        **kwargs
) -> hl.Table:
    """
    Convert an SDRF (Sample and Data Relationship Format) file to a Hail Table.

    This function imports an SDRF file, reshapes it from long to wide format,
    and converts it to a Hail Table that is persisted to disk.

    Parameters:
        sdrf_file (str): Path to the SDRF file.
        output_file (str): Path to save the Hail Table.
        keys (list): List of column names to use as keys for the Hail Table. Default is ['sample_id'].
        repartition (int): Number of partitions for the Hail Table. Default is 50.
        overwrite (bool): Whether to overwrite the output file if it exists. Default is False.
        **kwargs: Additional arguments to pass to pandas.read_csv when importing the SDRF file.

    Returns:
        hl.Table: Hail Table created from the SDRF file.

    Examples:
        >>> ht = convert_sdrf_to_hail_table('path/to/sdrf.tsv', 'output/path/ht')
        >>> ht = convert_sdrf_to_hail_table('path/to/sdrf.tsv', 'output/path/ht',
                                           keys=['subject_id', 'sample_id'],
                                           overwrite=True)
    """
    # Import SDRF file
    if keys is None:
        keys = ['sample_id']
    df_sdrf = _import_sdrf(sdrf_file, **kwargs)

    # Reshape SDRF DataFrame from long to wide format
    df_wide = _reshape_sdrf_long_to_wide_format(df_sdrf)

    # Convert the DataFrame to a Hail Table
    ht = (hl.Table.from_pandas(df_wide)
          .key_by(*keys)
          .repartition(repartition)
          .persist())

    # Write the Hail Table to a file
    if output_file is not None:
        ht = ht.checkpoint(output_file, overwrite=overwrite)

    return ht


def create_mt_from_expression_atlas_matrix(
    expression_matrix_path: str,
    output_path: str | None = None,
    delimiter: str = "\t",
    row_fields: dict | None = None,
    row_key: str = "GeneID",
    min_partitions: int = 50,
    force_bgz: bool = True,
    overwrite: bool = True,
    metadata_ht: hl.Table = None,
) -> hl.MatrixTable:
    """
    Creates a Hail MatrixTable from an expression atlas matrix file with optional metadata.

    Parameters
    ----------
    expression_matrix_path : str
        Path to the expression matrix file
    output_path : str, optional
        Path to save the resulting MatrixTable as a checkpoint
    delimiter : str, default="\t"
        Delimiter used in the expression matrix file
    row_fields : dict, default={"Gene ID": hl.str, "Gene Name": hl.str, "GeneID": hl.tstr}
        Dictionary mapping field names to Hail types for row fields
    row_key : str, default='GeneID'
        Key for rows, must match a field in row_fields
    min_partitions : int, default=50
        Minimum number of partitions for the imported MatrixTable
    force_bgz : bool, default=True
        Whether to force BGZF compression for the input file
    overwrite : bool, default=True
        Whether to allow overwriting existing files at the output path
    metadata_ht : hl.Table, optional
        Hail Table containing metadata for annotating columns

    Returns
    -------
    hl.MatrixTable
        MatrixTable generated from the expression matrix file

    Raises
    ------
    FileNotFoundError
        If the specified expression matrix file does not exist
    FileExistsError
        If the output path exists and overwrite is set to False
    ValueError
        If metadata_ht is provided but invalid or incompatible
    """
    # Validate input file existence
    if not os.path.exists(expression_matrix_path):
        raise FileNotFoundError(
            f"Expression matrix file not found: {expression_matrix_path}"
        )

    # Validate row_fields
    if row_fields is None:
        row_fields = {"Gene ID": hl.tstr, "Gene Name": hl.tstr, "GeneID": hl.tstr}

    # Check if output path exists when overwrite is False
    if output_path and os.path.exists(output_path) and not overwrite:
        raise FileExistsError(
            f"Output path already exists: {output_path}. Set overwrite=True to overwrite."
        )

    # Import the matrix table
    mt = hl.import_matrix_table(
        expression_matrix_path,
        delimiter=delimiter,
        row_fields=row_fields,
        row_key=row_key,
        entry_type=hl.tfloat32,
        min_partitions=min_partitions,
        force_bgz=force_bgz,
    )

    # Rename column identifier to more descriptive name
    mt = mt.rename({"col_id": "sample_id"})

    # Add metadata if provided
    if metadata_ht is not None:
        _validate_metadata(metadata_ht, mt)
        mt = mt.annotate_cols(metadata=metadata_ht[mt.col_key])

    # Save checkpoint if output path is specified
    if output_path:
        mt = mt.checkpoint(output_path, overwrite=overwrite)

    return mt


def _validate_metadata(metadata_ht, mt):
    """
    Validates that the metadata table is compatible with the matrix table.

    Parameters
    ----------
    metadata_ht : hl.Table
        Metadata table to validate
    mt : hl.MatrixTable
        Matrix table to validate against

    Raises
    ------
    ValueError
        If metadata_ht is not a Hail Table or has incompatible keys
    """
    if not isinstance(metadata_ht, hl.Table):
        raise ValueError("metadata_ht must be a Hail Table.")

    if str(metadata_ht.key) != str(mt.col_key):
        raise ValueError("metadata_ht must have the same key as the matrix table.")
