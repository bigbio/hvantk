"""
Flexible CLI to build individual AnnData (.h5ad) expression objects from raw inputs.

Examples:
  hvantk mkmatrix ucsc -e /path/to/expr.tsv -m /path/to/meta.tsv -o /out/expr.h5ad
  hvantk mkmatrix expression-atlas -e /path/to/matrix.tsv -s /path/to/atlas.sdrf.tsv -o /out/atlas.h5ad
  hvantk mkmatrix cptac --expression /path/expr.tsv --metadata /path/meta.tsv -o /out/cptac.h5ad
"""

import logging

import click

from hvantk.core.config import CONTEXT_SETTINGS

logger = logging.getLogger(__name__)


def _build_ucsc_ad(**kwargs):
    from hvantk.skills.ucsc_cellbrowser.builder import build_ucsc_ad

    return build_ucsc_ad(**kwargs)


def _build_expression_atlas_ad(**kwargs):
    from hvantk.skills.expression_atlas.builder import build_expression_atlas_ad

    return build_expression_atlas_ad(**kwargs)


def _build_cptac_ad(**kwargs):
    from hvantk.skills.cptac.expression.builder import build_cptac_ad

    return build_cptac_ad(**kwargs)


def _build_cptac_phospho_ad(**kwargs):
    from hvantk.skills.cptac.phospho.builder import build_cptac_phospho_ad

    return build_cptac_phospho_ad(**kwargs)


@click.group("mkmatrix", context_settings=CONTEXT_SETTINGS)
def mkmatrix_group():
    """Create a single AnnData expression object (.h5ad) from a raw input file."""
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
@click.option(
    "-o",
    "--output",
    "output_path",
    required=True,
    type=click.Path(),
    help="Output path for AnnData (.h5ad).",
)
@click.option("-g", "--gene-column", default="gene", show_default=True)
@click.option("-d", "--delimiter", default="\t", show_default=True)
@click.option(
    "--split-gene-field/--no-split-gene-field", default=True, show_default=True
)
@click.option(
    "--chunk-size",
    type=int,
    default=500,
    show_default=True,
    help="Gene rows per streaming chunk. Lower ⇒ less peak RAM, slower.",
)
@click.option(
    "--backed/--no-backed",
    default=None,
    help=(
        "Use the incremental backed-write builder (CSC columns appended to "
        ".h5ad on disk; no RAM materialization). Default: auto — backed if "
        "the input is larger than ~1 GiB."
    ),
)
@click.option(
    "--column-batch",
    type=int,
    default=64,
    show_default=True,
    help="Gene columns per backed-append block (backed mode only).",
)
@click.option("-w", "--overwrite", is_flag=True)
def mkmatrix_ucsc(
    expression_matrix,
    metadata,
    output_path,
    gene_column,
    delimiter,
    split_gene_field,
    chunk_size,
    backed,
    column_batch,
    overwrite,
):
    """Build an AnnData object from UCSC Cell Browser expression + metadata files."""
    logger.info("Building UCSC AnnData")
    adata = _build_ucsc_ad(
        expression_matrix_path=expression_matrix,
        metadata_path=metadata,
        output_path=output_path,
        gene_column=gene_column,
        delimiter=delimiter,
        split_gene_field=split_gene_field,
        chunk_size=chunk_size,
        backed=backed,
        column_batch=column_batch,
        overwrite=overwrite,
    )
    click.echo(f"AnnData created at {output_path}")
    click.echo(f"  Shape: {adata.shape[0]} obs x {adata.shape[1]} vars")


@mkmatrix_group.command("expression-atlas")
@click.option(
    "-e",
    "--expression-matrix",
    "expression_matrix",
    required=True,
    type=click.Path(exists=True),
)
@click.option("-s", "--sdrf", required=True, type=click.Path(exists=True))
@click.option(
    "-o",
    "--output",
    "output_path",
    required=True,
    type=click.Path(),
    help="Output path for AnnData (.h5ad).",
)
@click.option("-g", "--gene-column", default="Gene ID", show_default=True)
@click.option("-d", "--delimiter", default="\t", show_default=True)
@click.option("-w", "--overwrite", is_flag=True)
def mkmatrix_expression_atlas(
    expression_matrix,
    sdrf,
    output_path,
    gene_column,
    delimiter,
    overwrite,
):
    """Build an AnnData object from Expression Atlas matrix + SDRF metadata."""
    logger.info("Building Expression Atlas AnnData")
    adata = _build_expression_atlas_ad(
        expression_matrix_path=expression_matrix,
        sdrf_file=sdrf,
        output_path=output_path,
        gene_column=gene_column,
        delimiter=delimiter,
        overwrite=overwrite,
    )
    click.echo(f"AnnData created at {output_path}")
    click.echo(f"  Shape: {adata.shape[0]} obs x {adata.shape[1]} vars")


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
@click.option(
    "-o",
    "--output",
    "output_path",
    required=True,
    type=click.Path(),
    help="Output path for AnnData (.h5ad).",
)
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
@click.option("-w", "--overwrite", is_flag=True)
def mkmatrix_cptac(
    expression,
    metadata,
    output_path,
    gene_id_col,
    gene_name_col,
    sample_id_col,
    expression_col,
    overwrite,
):
    """Build a CPTAC AnnData object from expression and metadata tables."""
    logger.info("Building CPTAC AnnData")

    adata = _build_cptac_ad(
        expression_path=expression,
        metadata_path=metadata,
        output_path=output_path,
        gene_id_col=gene_id_col,
        gene_name_col=gene_name_col,
        sample_id_col=sample_id_col,
        expression_col=expression_col,
        overwrite=overwrite,
    )
    click.echo(f"AnnData created at {output_path}")
    click.echo(f"  Shape: {adata.shape[0]} obs x {adata.shape[1]} vars")


@mkmatrix_group.command("cptac-phospho")
@click.option(
    "-e",
    "--expression",
    required=True,
    type=click.Path(exists=True),
    help="Path to CPTAC phospho matrix CSV (sites x samples)",
)
@click.option(
    "-m",
    "--metadata",
    required=True,
    type=click.Path(exists=True),
    help="Path to CPTAC phospho metadata CSV",
)
@click.option(
    "-o",
    "--output",
    "output_path",
    required=True,
    type=click.Path(),
    help="Output path for AnnData (.h5ad).",
)
@click.option("--sid", "--sample-id-col", "sample_id_col", default="SampleID", show_default=True)
@click.option("-w", "--overwrite", is_flag=True)
def mkmatrix_cptac_phospho(
    expression,
    metadata,
    output_path,
    sample_id_col,
    overwrite,
):
    """Build a CPTAC phospho AnnData object from site intensities and metadata."""
    logger.info("Building CPTAC phospho AnnData")

    adata = _build_cptac_phospho_ad(
        expression_path=expression,
        metadata_path=metadata,
        output_path=output_path,
        sample_id_col=sample_id_col,
        overwrite=overwrite,
    )
    click.echo(f"AnnData created at {output_path}")
    click.echo(f"  Shape: {adata.shape[0]} obs x {adata.shape[1]} vars")
