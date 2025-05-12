import logging

import click
import hail as hl

from hvantk.htables.expression_atlas import (
    convert_sdrf_to_hail_table,
    create_mt_from_expression_atlas_matrix,
)
from hvantk.settings import CONTEXT_SETTINGS

logger = logging.getLogger(__name__)


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """Tools for working with Expression Atlas datasets."""
    pass


@click.command(
    "expression-atlas-matrix",
    short_help="Create Hail matrix table from Expression Atlas matrix and SDRF metadata",
)
@click.option(
    "--expression_matrix",
    "-e",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    default=None,
    help="Path to the expression matrix file (TSV or TSV.GZ file)",
)
@click.option(
    "--sdrf_file",
    "-s",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    default=None,
    help="Path to the SDRF metadata file (TSV format)",
)
@click.option(
    "--output_mt",
    "-o",
    required=True,
    type=click.Path(file_okay=False, dir_okay=True),
    default="data/expression_atlas_matrix.mt",
    help="Output path for the Hail matrix table",
)
@click.option(
    "--gene_column",
    default="Gene ID",
    help="Column name for gene identifiers in the expression matrix (default: 'Gene ID')",
)
@click.option(
    "--sample_id_column",
    default="sample_id",
    help="Column name for sample IDs in the metadata (default: 'sample_id')",
)
@click.option(
    "--delimiter",
    "-d",
    default="\t",
    type=str,
    help="Delimiter for the expression matrix file (default: tab)",
)
@click.option(
    "--min_partitions",
    "-p",
    default=50,
    type=int,
    help="Minimum number of partitions for the Matrix Table",
)
@click.option(
    "--force_bgz",
    is_flag=True,
    help="Force bgz compression for the input matrix expression file",
)
@click.option(
    "--overwrite",
    is_flag=True,
    help="Overwrite existing files at output path",
)
def make_expression_atlas_matrix_table(
    sdrf_file,
    expression_matrix,
    output_mt,
    gene_column,
    sample_id_column,
    delimiter,
    min_partitions,
    force_bgz,
    overwrite,
):
    """
    Create a Hail Matrix Table from Expression Atlas SDRF metadata and expression matrix.

    :param sdrf_file: Path to the Expression Atlas SDRF metadata file in TSV format.
    :param expression_matrix: Path to the Expression Atlas gene expression matrix file (TSV format).
    :param output_mt: Output directory path for the generated Hail Matrix Table.
    :param gene_column: Column name in the expression matrix file corresponding to genes.
    :param sample_id_column: Column name for sample IDs in the metadata.
    :param delimiter: Delimiter used in the expression matrix file (default is tab).
    :param min_partitions: The desired minimum number of partitions to use for the Matrix Table.
    :param force_bgz: Specifies whether to force BGZ compression on the input matrix file.
    :param overwrite: Indicates if existing files at the output path should be overwritten.

    :return: Hail Matrix Table created from the given Expression Atlas metadata and expression matrix.
    :rtype: hl.MatrixTable
    """

    metadata_output = f"{output_mt}.metadata.ht"
    logger.info(f"Converting SDRF metadata file to Hail Table: {sdrf_file}")

    # Convert SDRF to Hail Table using existing function
    metadata_ht = convert_sdrf_to_hail_table(
        sdrf_file=sdrf_file,
        output_file=metadata_output,
        keys=[sample_id_column],
        repartition=min_partitions,
        overwrite=overwrite,
    )

    logger.info(f"Successfully created metadata table at: {metadata_output}")
    logger.info(f"Creating matrix table from expression matrix: {expression_matrix}")

    # Use the create_mt_from_expression_atlas_matrix function for matrix table creation
    mt = create_mt_from_expression_atlas_matrix(
        expression_matrix_path=expression_matrix,
        output_path=output_mt,
        delimiter=delimiter,
        row_fields={gene_column: hl.tstr},
        row_key=gene_column,
        min_partitions=min_partitions,
        force_bgz=force_bgz,
        overwrite=overwrite,
        metadata_ht=metadata_ht,
    )

    # Log information about the created matrix table
    num_rows = mt.count_rows()
    num_cols = mt.count_cols()
    logger.info(f"Successfully created matrix table at: {output_mt}")
    logger.info(
        f"Matrix table dimensions: {num_rows} rows (genes) × {num_cols} columns (samples)"
    )

    # Print nicely formatted mt.describe() to the console
    logger.info("Matrix Table description:")
    logger.info(mt.describe())

    return mt


# Add the command to the CLI group
cli.add_command(make_expression_atlas_matrix_table)

if __name__ == "__main__":
    cli()
