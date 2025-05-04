import click
import logging

logger = logging.getLogger(__name__)

from hvantk.settings import CONTEXT_SETTINGS
from hvantk.commands.make_annotation_tables_cli import make_annotation_tables_cli
from hvantk.commands.ucsc_downloader import ucsc_downloader


# Main CLI entry point for the package (hvantk)


@click.group(
    "hvantk",
    help="A python package for gene and variant annotation.",
    context_settings=CONTEXT_SETTINGS,
)
def cli():
    """A python package for gene and variant annotation."""
    logger.info("Starting hvantk CLI")
    pass


cli.add_command(ucsc_downloader)
cli.add_command(make_annotation_tables_cli)


def main():
    logger.info("Running main function")
    cli()
    logger.info("Main function completed")


if __name__ == "__main__":
    main()
