import click
import logging

logger = logging.getLogger(__name__)

from hvantk.core.config import CONTEXT_SETTINGS
from hvantk.commands.ucsc_downloader import ucsc_downloader
from hvantk.commands.expression_atlas_downloader import download_experiments
from hvantk.commands.make_table_cli import mktable_group
from hvantk.commands.make_table_batch_cli import mktable_batch_cli
from hvantk.commands.make_matrix_cli import mkmatrix_group
from hvantk.commands.make_matrix_batch_cli import mkmatrix_batch_cli
from hvantk.commands.catalog_cli import catalog
from hvantk.commands.hgc_cli import hgc_group


# Main CLI entry point for the package (hvantk)


@click.group(
    "hvantk",
    help="A python package for gene and variant annotation with joint genotyping capabilities.",
    context_settings=CONTEXT_SETTINGS,
)
def cli():
    """
    Entry point for the hvantk command-line interface.

    Serves as the root CLI group for gene and variant annotation commands
    with integrated joint genotyping workflows.
    """
    logger.info("Starting hvantk CLI")
    pass


cli.add_command(ucsc_downloader)
cli.add_command(download_experiments)
cli.add_command(mktable_group)  # per-table builder
cli.add_command(mktable_batch_cli)  # batch builder from recipe
cli.add_command(mkmatrix_group)  # per-matrix builder
cli.add_command(mkmatrix_batch_cli)  # batch matrix builder from recipe
cli.add_command(catalog)  # catalog operations
cli.add_command(hgc_group)  # HGC joint genotyping commands


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
