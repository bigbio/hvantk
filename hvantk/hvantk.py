import click
import logging

logger = logging.getLogger(__name__)

from hvantk.core.config import CONTEXT_SETTINGS
from hvantk.commands.ucsc_downloader import ucsc_downloader
from hvantk.commands.expression_atlas_downloader import download_experiments
from hvantk.commands.clingen_downloader import clingen_downloader
from hvantk.commands.hgnc_downloader import hgnc_downloader
from hvantk.commands.clinvar_downloader import clinvar_downloader
from hvantk.commands.make_table_cli import mktable_group
from hvantk.commands.make_table_batch_cli import mktable_batch_cli
from hvantk.commands.make_matrix_cli import mkmatrix_group
from hvantk.commands.make_matrix_batch_cli import mkmatrix_batch_cli
from hvantk.commands.catalog_cli import catalog
from hvantk.commands.hgc import hgc_group
from hvantk.commands.psroc_cli import psroc_cmd
from hvantk.commands.enrichex_cli import enrichex_group
from hvantk.commands.ancestry_cli import ancestry_inference_cmd
from hvantk.commands.build_1k_genome_cli import build_1k_genome_cmd

# Main CLI entry point for the package (hvantk)


def setup_logging(verbosity: int = 0, log_file: str | None = None):
    """Configure centralized logging for the hvantk CLI.

    Parameters
    ----------
    verbosity : int
        Verbosity level: 0 = WARNING (default), 1 = INFO, 2+ = DEBUG.
    log_file : str or None
        Optional path to a file where log output will be written.
    """
    level = {0: logging.WARNING, 1: logging.INFO, 2: logging.DEBUG}.get(
        verbosity, logging.DEBUG
    )

    handlers: list[logging.Handler] = [logging.StreamHandler()]
    if log_file:
        handlers.append(logging.FileHandler(log_file))

    logging.basicConfig(
        level=level,
        format="%(asctime)s %(name)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=handlers,
        force=True,
    )


@click.group(
    "hvantk",
    help="A python package for gene and variant annotation with joint genotyping capabilities.",
    context_settings=CONTEXT_SETTINGS,
)
@click.option(
    "-v",
    "--verbose",
    count=True,
    help="Increase verbosity (-v: INFO, -vv: DEBUG)",
)
@click.option(
    "--log-file",
    type=click.Path(),
    default=None,
    help="Write logs to file",
)
def cli(verbose, log_file):
    """
    Entry point for the hvantk command-line interface.

    Serves as the root CLI group for gene and variant annotation commands
    with integrated joint genotyping workflows.
    """
    setup_logging(verbose, log_file)
    logger.info("Starting hvantk CLI")


cli.add_command(ucsc_downloader)
cli.add_command(download_experiments)
cli.add_command(clingen_downloader)
cli.add_command(hgnc_downloader)
cli.add_command(clinvar_downloader)
cli.add_command(mktable_group)  # per-table builder
cli.add_command(mktable_batch_cli)  # batch builder from recipe
cli.add_command(mkmatrix_group)  # per-matrix builder
cli.add_command(mkmatrix_batch_cli)  # batch matrix builder from recipe
cli.add_command(catalog)  # catalog operations
cli.add_command(hgc_group)  # HGC joint genotyping commands
cli.add_command(psroc_cmd)  # PSROC prediction score ROC analysis
cli.add_command(enrichex_group)  # EnrichEx gene set enrichment analysis
cli.add_command(ancestry_inference_cmd)  # Ancestry inference pipeline
cli.add_command(build_1k_genome_cmd)  # 1000 Genomes MatrixTable builder

# Lazy import: convert-bgz is a new command that may not be installed yet
# in existing editable installs until `poetry install` is re-run.
try:
    from hvantk.commands.convert_bgz_cli import convert_bgz_cmd

    cli.add_command(convert_bgz_cmd)  # GZ → BGZF conversion utility
except ImportError as e:
    logger.debug("Could not load convert-bgz command (%s): %s", "hvantk.commands.convert_bgz_cli", e)


def main():
    """
    Runs the main entry point for the hvantk CLI application.

    Invokes the command-line interface to process user commands and logs the start and completion of the main function.
    """
    logger.info("Running main function")
    cli()
    logger.info("Main function completed")


if __name__ == "__main__":
    main()
