import logging

import click
import hail as hl

from hvantk.settings import CONTEXT_SETTINGS
from hvantk.utils.constants import (
    UCSC_CELL_ID_COLUMN,
    UCSC_GENE_COLUMN
)
from hvantk.utils.ucsc import (
    convert_ucsc_metadata_to_hail_table,
    create_mt_from_ucsc_expression_matrix
)

logger = logging.getLogger(__name__)


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """A package for gene and variant annotation."""
    pass


@click.command("ucsc-matrix", short_help="Create Hail matrix table from UCSC metadata and expression matrix")
@click.option('--expression_matrix',
              '-e',
              required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=False),
              default=None,
              help='Path to the expression matrix file (Block-compressed TSV file)')
@click.option('--metadata',
              '-m',
              required=True,
              type=click.Path(exists=True, file_okay=True, dir_okay=False),
              default=None,
              help='Path to the metadata file (TSV format)'
)
@click.option('--output_mt',
              '-o',
              required=True,
              type=click.Path(file_okay=False, dir_okay=True),
              default='ucsc_expression_matrix.mt',
              help='Output path for the Hail matrix table'
)
@click.option('--gene_column', 
              default=UCSC_GENE_COLUMN,
              help='Column name for gene in the expression matrix (default: UCSC_GENE_COLUMN)')
@click.option('--split_gene_field', 
              default=True, 
              help='Split gene field and use first element(e.g., A|B -> A)')
@click.option('--metadata_index_col', 
              default=0,
              help='Index column in metadata file (default: 0)')
@click.option('--delimiter',
                '-d',
                default='\t',
                type=str,
                help='Delimiter for the expression matrix file (default: tab)')
@click.option('--min_partitions', 
              '-p',
              default=50, 
              type=int,
              help='Minimum number of partitions for the Matrix Table')
@click.option('--force_bgz', 
              default=True, 
              help='Force bgz compression for the input matrix expression file')
@click.option('--overwrite', 
              default=True, 
              help='Overwrite existing files at output path')
def make_ucsc_matrix_table(metadata,
                           expression_matrix,
                           output_mt,
                           gene_column,
                           metadata_index_col,
                           delimiter,
                           min_partitions,
                           force_bgz,
                           split_gene_field,
                           overwrite):
    """
    Create a Hail Matrix Table from UCSC metadata and expression matrix.

    :param ctx: The click context object, passed automatically by the Click CLI framework.
    :param metadata: Path to the UCSC metadata file in TSV format.
    :param expression_matrix: Path to the UCSC gene expression matrix file (block-compressed TSV format).
    :param output_mt: Output directory path for the generated Hail Matrix Table.
    :param gene_column: Column name in the expression matrix file corresponding to genes.
    :param metadata_index_col: Index column in the metadata file to be used as the primary key.
    :param delimiter: Delimiter used in the expression matrix file (default is tab).
    :param min_partitions: The desired minimum number of partitions to use for the Matrix Table.
    :param force_bgz: Specifies whether to force BGZ compression on the input matrix file.
    :param split_gene_field: Determines whether to split the gene field in the matrix
                             and use the first element (e.g., from "A|B", use "A").
    :param overwrite: Indicates if existing files at the output path should be overwritten.

    :return: Hail Matrix Table created from the given UCSC metadata and expression matrix.
    :rtype: hl.MatrixTable
    """
    
    logger.info(f"Converting metadata file: {metadata}")
    # Convert metadata to Hail Table
    metadata_ht = convert_ucsc_metadata_to_hail_table(
        metadata_path=metadata,
        sep=delimiter,
        index_col=metadata_index_col,
        index_name=UCSC_CELL_ID_COLUMN
    )
    
    logger.info(f"Creating matrix table from expression matrix: {expression_matrix}")
    # Create matrix table
    mt = create_mt_from_ucsc_expression_matrix(
        expression_matrix_path=expression_matrix,
        output_path=output_mt,
        delimiter=delimiter,
        row_fields={gene_column: hl.tstr},
        row_key=gene_column,
        split_gene_field=split_gene_field,
        min_partitions=min_partitions,
        force_bgz=force_bgz,
        overwrite=overwrite,
        metadata_ht=metadata_ht
    )
    
    num_rows = mt.count_rows()
    num_cols = mt.count_cols()
    logger.info(f"Successfully created matrix table at: {output_mt}")
    logger.info(f"Matrix table dimensions: {num_rows} rows (genes) × {num_cols} columns (cells)")

    # print nicely formatted mt.describe() to the console
    logger.info("Matrix Table description:")
    logger.info(mt.describe())
    
    return mt


if __name__ == "__main__":
    cli()