import click
import logging

logger = logging.getLogger(__name__)

from hvantk.core.config import CONTEXT_SETTINGS
from hvantk.tools.plugins.download_cli import download_group
from hvantk.tools.infra.utils_cli import utils_group
from hvantk.tools.infra.catalog_cli import catalog as catalog_group
from hvantk.tools.genesets.genesets_cli import genesets_group
from hvantk.tools.hgc import hgc_group
from hvantk.tools.ptm.psroc_cli import psroc_cmd
from hvantk.tools.enrichex import enrichex_group
from hvantk.tools.ancestry.ancestry_cli import ancestry_inference_cmd
from hvantk.tools.expression.summarize_expression_cli import expression_group
from hvantk.tools.ptm.ptm_cli import ptm_group
from hvantk.tools.qtl.qtlcascade_cli import qtlcascade_group
from hvantk.tools.plugins.plugins_cli import plugins_group
from hvantk.tools.plugins.tools_cli import tools_group
from hvantk.tools.plugins.drift_cli import drift_cmd
from hvantk.tools.plugins.reprocess_cli import reprocess_cmd
from hvantk.tools.rerank import rerank_cmd
from hvantk.tools.annotation.annotate_cli import annotate_group
from hvantk.tools.cohort.cohort_cli import cohort_group

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


cli.add_command(download_group)
cli.add_command(utils_group)
cli.add_command(catalog_group)
cli.add_command(genesets_group)
cli.add_command(hgc_group)
cli.add_command(psroc_cmd)
cli.add_command(enrichex_group)
cli.add_command(ancestry_inference_cmd)
cli.add_command(expression_group)
cli.add_command(ptm_group)
cli.add_command(qtlcascade_group)
cli.add_command(plugins_group)
cli.add_command(tools_group)
cli.add_command(drift_cmd)
cli.add_command(reprocess_cmd)
cli.add_command(rerank_cmd)
cli.add_command(annotate_group)
cli.add_command(cohort_group)


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
