"""
Generate annotation tables from multiple (raw) sources

"""

import click

import logging

logger = logging.getLogger(__name__)

from hvantk.core.config import CONTEXT_SETTINGS, RAW_DATA_PATH, RAW_DATA_PATHS, set_raw_data_path

output_dir_default = f"{RAW_DATA_PATH}/annotation_tables"


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """A package for gene and variant annotation."""
    pass


def make_annotation_tables_from_raw_sources(
    raw_data_path: str,
    ccr: bool = False,  # TODO: implement CCR table
    interactome: bool = False,
    temporal_rnaseq: bool = False,
    clinvar: bool = False,
    gevir: bool = False,
    scell_heart_deg: bool = False,
    hca_rnaseq: bool = False,
    gene_ensembl: bool = False,
    gnomad_metrics: bool = False,
    output_dir: str = output_dir_default,
    default_ref_genome: str = "GRCh38",
):
    # Import table creators lazily to avoid heavy Hail import at CLI load time
    from hvantk.tables.creators import (
        create_interactome_tb,
        create_clinvar_tb,
        create_gevir_tb,
        create_gnomad_constraint_gene_metrics_tb,
    )

    # set the raw data path
    """
    Generates selected gene and variant annotation tables from raw data sources.

    For each enabled annotation type, creates the corresponding Hail table from raw data
    and saves it to the specified output directory. Only interactome, ClinVar, GEVIR, and
    gnomAD metrics tables are currently implemented; other options are accepted but not used.

    Args:
        raw_data_path: Path to the directory containing raw annotation data.
        output_dir: Directory where generated annotation tables will be saved.
        default_ref_genome: Reference genome build identifier (e.g., "GRCh38").

    Note:
        Only the interactome, ClinVar, GEVIR, and gnomAD metrics tables are generated.
        Other flags are reserved for future implementation.
    """
    logger.info(f"Setting raw data path to {raw_data_path}")
    set_raw_data_path(raw_data_path)

    if interactome:
        logger.info("Creating interactome table")
        input_path = RAW_DATA_PATHS["interactome_path"]
        output_path = f"{output_dir}/interactome.{default_ref_genome}.ht"
        create_interactome_tb(
            input_path=input_path,
            output_path=output_path,
            overwrite=True,
        )
        logger.info(
            f"Interactome table created at {output_path}"
        )

    if clinvar:
        logger.info("Creating ClinVar table")
        input_path = RAW_DATA_PATHS["clinvar_path"]
        output_path = f"{output_dir}/clinvar.{default_ref_genome}.ht"
        create_clinvar_tb(
            input_path=input_path,
            output_path=output_path,
            overwrite=True,
            reference_genome=default_ref_genome,
        )
        logger.info(
            f"ClinVar table created at {output_path}"
        )

    if gevir:
        logger.info("Creating GEVIR table")
        input_path = RAW_DATA_PATHS["gevir_path"]
        output_path = f"{output_dir}/gevir.metrics.ht"
        create_gevir_tb(
            input_path=input_path,
            output_path=output_path,
            overwrite=True,
        )
        logger.info(f"GEVIR table created at {output_path}")

    if gnomad_metrics:
        logger.info("Creating gnomAD metrics table")
        input_path = RAW_DATA_PATHS["gnomad_metrics_path"]
        output_path = f"{output_dir}/gnomad.metrics.ht"
        create_gnomad_constraint_gene_metrics_tb(
            input_path=input_path,
            output_path=output_path,
            overwrite=True,
        )
        logger.info(f"gnomAD metrics table created at {output_path}")


@click.command("mktables", short_help="Create annotation tables from raw sources.")
@click.option(
    "--raw_data_path",
    default=RAW_DATA_PATH,
    type=str,
    required=True,
    help="Path to raw data directory.",
)
@click.option(
    "--output_dir",
    default=output_dir_default,
    type=str,
    required=True,
    help="Output directory to copy created Hail tables",
)
@click.option("--ccr", is_flag=True, help="Create/update CCR table from source.")
@click.option(
    "--interactome", is_flag=True, help="Create/update CCR table from source."
)
@click.option(
    "--clinvar", is_flag=True, help="Create/update Clinvar table from source."
)
@click.option(
    "--gevir", is_flag=True, help="Create/update GeVIR score table from raw source."
)
@click.option(
    "--scell_heart_deg",
    is_flag=True,
    help="Create/update table with DEGs from cardiac-specific cell clusters",
)
@click.option(
    "--hca_rnaseq",
    is_flag=True,
    help="Create/update table with gene/cell expression levels from HCA dataset (UCSC)",
)
@click.option(
    "--gene_ensembl",
    is_flag=True,
    help="Create/update gene annotation table from Ensembl.",
)
@click.option(
    "--gnomad_metrics",
    is_flag=True,
    help="Create/update transcript-specific constraint metrics from gnomad database",
)
@click.option(
    "--default_ref_genome",
    default="GRCh38",
    type=str,
    help="Default reference genome to start Hail. Only GRCh38 is supported for now.",
)
@click.pass_context
def make_annotation_tables_cli(
    ctx,
    raw_data_path,
    output_dir,
    ccr,
    interactome,
    clinvar,
    gevir,
    scell_heart_deg,
    hca_rnaseq,
    gene_ensembl,
    gnomad_metrics,
    default_ref_genome,
):

    # exit if no flat parameter is set
    """
    Handles the CLI command for generating annotation tables from raw data sources.

    Validates that at least one annotation table flag is set, then initiates the creation
    of the selected tables using the provided raw data path, output directory, and
    reference genome. Aborts execution if no table creation flags are specified.
    """
    if not any(
        [
            ccr,
            interactome,
            clinvar,
            gevir,
            scell_heart_deg,
            hca_rnaseq,
            gene_ensembl,
            gnomad_metrics,
        ]
    ):
        click.echo(
            "No flag set. Please set at least one flag to create/update a table."
        )
        ctx.abort()

    logger.info("Starting make_annotation_tables_from_raw_sources")
    make_annotation_tables_from_raw_sources(
        raw_data_path,
        ccr,
        interactome,
        clinvar,
        gevir,
        scell_heart_deg,
        hca_rnaseq,
        gene_ensembl,
        gnomad_metrics,
        output_dir,
        default_ref_genome,
    )
    logger.info("make_annotation_tables_from_raw_sources completed")


if __name__ == "__main__":
    logger.info("Starting make_annotation_tables_cli")
    cli()
    logger.info("make_annotation_tables_cli completed")
