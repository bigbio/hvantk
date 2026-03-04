"""
Flexible CLI to build individual Hail MatrixTables from raw inputs.

Examples:
  hvantk mkmatrix ucsc -e /path/to/expr.tsv.bgz -m /path/to/meta.tsv -o /out/mt
  hvantk mkmatrix expression-atlas -e /path/to/matrix.tsv -s /path/to/atlas.sdrf.tsv -o /out/atlas.mt
  hvantk mkmatrix cptac --expression /path/expr.tsv --metadata /path/meta.tsv -o /out/cptac.mt
"""

import logging

import click

from hvantk.core.config import CONTEXT_SETTINGS

logger = logging.getLogger(__name__)


def _build_ucsc_mt(**kwargs):
    from hvantk.tables.matrix_builders import build_ucsc_mt

    return build_ucsc_mt(**kwargs)


def _build_expression_atlas_mt(**kwargs):
    from hvantk.tables.matrix_builders import build_expression_atlas_mt

    return build_expression_atlas_mt(**kwargs)


def _build_cptac_mt(**kwargs):
    from hvantk.tables.matrix_builders import build_cptac_mt

    return build_cptac_mt(**kwargs)


@click.group("mkmatrix", context_settings=CONTEXT_SETTINGS)
def mkmatrix_group():
    """Create a single Hail MatrixTable from a raw input file."""
    pass


@mkmatrix_group.command("ucsc")
@click.option(
    "-e",
    "--expression-matrix",
    "expression_matrix",
    required=True,
    type=click.Path(exists=True),
)
@click.option("-m", "--metadata", required=True, type=click.Path(exists=True))
@click.option("-o", "--output-mt", "output_mt", required=True, type=click.Path())
@click.option("-g", "--gene-column", default="gene", show_default=True)
@click.option(
    "--ix",
    "--metadata-index-col",
    "metadata_index_col",
    default=0,
    type=int,
    show_default=True,
)
@click.option("-d", "--delimiter", default="\t", show_default=True)
@click.option("-p", "--min-partitions", default=50, type=int, show_default=True)
@click.option("--force-bgz/--no-force-bgz", default=True, show_default=True)
@click.option(
    "--split-gene-field/--no-split-gene-field", default=True, show_default=True
)
@click.option("-w", "--overwrite", is_flag=True)
@click.option(
    "--auto-convert-bgz",
    is_flag=True,
    help="Automatically convert plain gzip files to BGZF before import",
)
def mkmatrix_ucsc(
    expression_matrix,
    metadata,
    output_mt,
    gene_column,
    metadata_index_col,
    delimiter,
    min_partitions,
    force_bgz,
    split_gene_field,
    overwrite,
    auto_convert_bgz,
):
    """Build a MatrixTable from UCSC Cell Browser expression + metadata files."""
    logger.info("Building UCSC MatrixTable")
    mt = _build_ucsc_mt(
        expression_matrix_path=expression_matrix,
        metadata_path=metadata,
        output_mt=output_mt,
        gene_column=gene_column,
        metadata_index_col=metadata_index_col,
        delimiter=delimiter,
        min_partitions=min_partitions,
        force_bgz=force_bgz,
        split_gene_field=split_gene_field,
        overwrite=overwrite,
        auto_convert_bgz=auto_convert_bgz,
    )
    click.echo(f"MatrixTable created at {output_mt}")
    mt.describe()


@mkmatrix_group.command("expression-atlas")
@click.option(
    "-e",
    "--expression-matrix",
    "expression_matrix",
    required=True,
    type=click.Path(exists=True),
)
@click.option("-s", "--sdrf", required=True, type=click.Path(exists=True))
@click.option("-o", "--output-mt", "output_mt", required=True, type=click.Path())
@click.option("-g", "--gene-column", default="Gene ID", show_default=True)
@click.option(
    "--sid",
    "--sample-id-column",
    "sample_id_column",
    default="sample_id",
    show_default=True,
)
@click.option("-d", "--delimiter", default="\t", show_default=True)
@click.option("-p", "--min-partitions", default=50, type=int, show_default=True)
@click.option("--force-bgz/--no-force-bgz", default=False, show_default=True)
@click.option("-w", "--overwrite", is_flag=True)
@click.option(
    "--auto-convert-bgz",
    is_flag=True,
    help="Automatically convert plain gzip files to BGZF before import",
)
def mkmatrix_expression_atlas(
    expression_matrix,
    sdrf,
    output_mt,
    gene_column,
    sample_id_column,
    delimiter,
    min_partitions,
    force_bgz,
    overwrite,
    auto_convert_bgz,
):
    """Build a MatrixTable from Expression Atlas matrix + SDRF metadata."""
    logger.info("Building Expression Atlas MatrixTable")
    mt = _build_expression_atlas_mt(
        expression_matrix_path=expression_matrix,
        sdrf_file=sdrf,
        output_mt=output_mt,
        gene_column=gene_column,
        sample_id_column=sample_id_column,
        delimiter=delimiter,
        min_partitions=min_partitions,
        force_bgz=force_bgz,
        overwrite=overwrite,
        auto_convert_bgz=auto_convert_bgz,
    )
    click.echo(f"MatrixTable created at {output_mt}")
    mt.describe()


@mkmatrix_group.command("cptac")
@click.option(
    "-e",
    "--expression",
    required=True,
    type=click.Path(exists=True),
    help="Path to CPTAC expression TSV/CSV",
)
@click.option(
    "-m",
    "--metadata",
    required=True,
    type=click.Path(exists=True),
    help="Path to CPTAC metadata TSV/CSV",
)
@click.option("-o", "--output-mt", "output_mt", required=True, type=click.Path())
@click.option("-g", "--gene-id-col", default="GeneID", show_default=True)
@click.option("-n", "--gene-name-col", default="Gene Name", show_default=True)
@click.option(
    "--sid",
    "--sample-id-col",
    "sample_id_col",
    default="SampleID",
    show_default=True,
)
@click.option("-x", "--expression-col", default="Expression", show_default=True)
@click.option(
    "-c",
    "--categorical-cols",
    default=None,
    help="Comma-separated categorical metadata columns",
)
@click.option(
    "-u",
    "--numeric-cols", default=None, help="Comma-separated numeric metadata columns"
)
@click.option("-w", "--overwrite", is_flag=True)
def mkmatrix_cptac(
    expression,
    metadata,
    output_mt,
    gene_id_col,
    gene_name_col,
    sample_id_col,
    expression_col,
    categorical_cols,
    numeric_cols,
    overwrite,
):
    """Build a CPTAC MatrixTable from expression and metadata tables."""
    logger.info("Building CPTAC MatrixTable")

    def _split(val):
        if not val:
            return None
        parts = [x.strip() for x in val.split(",") if x.strip()]
        return parts or None

    mt = _build_cptac_mt(
        expression_path=expression,
        metadata_path=metadata,
        output_mt=output_mt,
        gene_id_col=gene_id_col,
        gene_name_col=gene_name_col,
        sample_id_col=sample_id_col,
        expression_col=expression_col,
        categorical_cols=_split(categorical_cols),
        numeric_cols=_split(numeric_cols),
        overwrite=overwrite,
    )
    click.echo(f"MatrixTable created at {output_mt}")
    mt.describe()
