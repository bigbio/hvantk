from __future__ import annotations

import anndata as ad
import pandas as pd
import os

__all__ = [
    "convert_sdrf_to_dataframe",
    "create_anndata_from_expression_atlas",
]


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
            "names": [
                "accession",
                "unused",
                "sample_id",
                "column_type",
                "column_name",
                "column_value",
            ],
        }

        # Update defaults with any provided kwargs
        read_params = {**default_params, **kwargs}

        # Infer number of columns if names not provided in kwargs
        if "names" not in kwargs:
            with open(sdrf_file) as f:
                first_line = f.readline().strip()
                num_cols = len(first_line.split(read_params["sep"]))
                read_params["usecols"] = range(min(6, num_cols))

        # Read file
        sdrf_df = pd.read_csv(sdrf_file, **read_params)

        # Clean column names
        sdrf_df.columns = sdrf_df.columns.str.strip()

        # Drop unused columns
        sdrf_df.drop(columns=["unused"], inplace=True, errors="ignore")

        return sdrf_df

    except FileNotFoundError as e:
        raise FileNotFoundError(f"SDRF file not found: {sdrf_file}") from e
    except pd.errors.EmptyDataError as e:
        raise pd.errors.EmptyDataError(f"SDRF file is empty: {sdrf_file}") from e


def _reshape_sdrf_long_to_wide_format(
    df_sdrf: pd.DataFrame,
    include_only_factors: bool = False,
    include_only_characteristic: bool = False,
) -> pd.DataFrame:
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
        raise ValueError(
            "Cannot set both include_only_factors and include_only_characteristic to True"
        )

    df = df_sdrf.copy()

    # Filter rows based on column_type
    if include_only_characteristic:
        df = df[df["column_type"] == "characteristic"]
    elif include_only_factors:
        df = df[df["column_type"] == "factor"]
    # else: use all column types

    # Handle duplicate sample_id/column_name combinations by taking last value
    pivot_df = df[["sample_id", "column_name", "column_value"]].drop_duplicates(
        subset=["sample_id", "column_name"], keep="last"
    )

    # Reshape from long to wide format
    wide_df = pivot_df.pivot(
        index="sample_id", columns="column_name", values="column_value"
    )

    # Reset index to make sample_id a column
    wide_df.reset_index(inplace=True)

    # Clean column names (replace spaces and special chars with underscores)
    wide_df.columns = [
        str(col).strip().replace(" ", "_").replace("(", "").replace(")", "")
        for col in wide_df.columns
    ]

    return wide_df


def convert_sdrf_to_dataframe(sdrf_file: str, **kwargs) -> pd.DataFrame:
    """Parse an SDRF file into a wide-format DataFrame indexed by sample_id.

    Uses the existing ``_import_sdrf`` and ``_reshape_sdrf_long_to_wide_format``
    helpers to read the raw SDRF and pivot it so that each row represents one
    sample and columns represent its characteristics/factors.

    Parameters
    ----------
    sdrf_file : str
        Path to the SDRF file.
    **kwargs
        Extra keyword arguments forwarded to ``_import_sdrf`` (and ultimately
        to ``pandas.read_csv``).

    Returns
    -------
    pd.DataFrame
        Wide-format DataFrame with one row per sample, indexed by
        ``sample_id``.
    """
    df_long = _import_sdrf(sdrf_file, **kwargs)
    df_wide = _reshape_sdrf_long_to_wide_format(df_long)
    df_wide = df_wide.set_index("sample_id")
    return df_wide


def create_anndata_from_expression_atlas(
    expression_matrix_path: str,
    metadata_df: pd.DataFrame = None,
    gene_id_column: str = "Gene ID",
    gene_name_column: str = "Gene Name",
    delimiter: str = "\t",
) -> "ad.AnnData":
    """Create an AnnData object from an Expression Atlas expression matrix.

    The input TSV is gene-centric (rows = genes, columns = samples) with
    ``gene_id_column`` and ``gene_name_column`` as the first two columns
    followed by one column per sample.  The matrix is transposed so that
    observations are samples and variables are genes, following the AnnData
    convention.

    Parameters
    ----------
    expression_matrix_path : str
        Path to the expression matrix TSV file.
    metadata_df : pd.DataFrame, optional
        Sample metadata indexed by sample_id.  Columns are joined into
        ``obs``.
    gene_id_column : str
        Name of the gene identifier column (default ``"Gene ID"``).
    gene_name_column : str
        Name of the gene name column (default ``"Gene Name"``).
    delimiter : str
        Column delimiter (default tab).

    Returns
    -------
    ad.AnnData
        Expression AnnData with genes in ``var`` and samples in ``obs``.
    """
    import anndata as ad
    import numpy as np

    df = pd.read_csv(expression_matrix_path, sep=delimiter)

    # Extract gene annotations
    gene_ids = df[gene_id_column].values
    gene_names = df[gene_name_column].values if gene_name_column in df.columns else None

    # Sample columns are everything except the gene id/name columns
    non_sample_cols = {gene_id_column, gene_name_column} & set(df.columns)
    sample_cols = [c for c in df.columns if c not in non_sample_cols]

    # Build expression matrix (samples x genes)
    X = df[sample_cols].values.T.astype(np.float32)

    # var DataFrame (genes)
    var = pd.DataFrame(index=pd.Index(gene_ids, name="gene_id"))
    if gene_names is not None:
        var[gene_name_column] = gene_names

    # obs DataFrame (samples)
    obs = pd.DataFrame(index=pd.Index(sample_cols, name="sample_id"))

    adata = ad.AnnData(X=X, obs=obs, var=var)

    # Join metadata into obs if provided
    if metadata_df is not None:
        common = adata.obs.index.intersection(metadata_df.index)
        for col in metadata_df.columns:
            adata.obs[col] = metadata_df.loc[common, col].reindex(adata.obs.index)

    return adata


