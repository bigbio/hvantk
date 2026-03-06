"""
EnrichEx CLI Commands - Gene Set Enrichment Analysis

This module provides CLI commands for gene set enrichment analysis:
- overlap: Test gene list overlap enrichment (Fisher's exact)
- burden: Case-control burden analysis (Hail-native regression)

Note: Gene set extraction has been moved to `hvantk utils extract genesets`.
"""

import logging

import click

from hvantk.core.config import CONTEXT_SETTINGS

logger = logging.getLogger(__name__)


@click.group(
    name="enrichex",
    help="Gene set enrichment analysis commands.",
    context_settings=CONTEXT_SETTINGS,
)
@click.option(
    "--log-level",
    default="INFO",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR"]),
    help="Set logging level",
)
@click.pass_context
def enrichex_group(ctx, log_level):
    """EnrichEx command group for gene set enrichment analysis.

    This command group provides tools for testing gene set enrichment
    using overlap analysis (Fisher's exact test) and burden testing
    (case-control association tests).

    Examples
    --------
    Test gene list enrichment:
        hvantk enrichex overlap -g genes.txt -s gene_sets.json -o results.tsv

    Run burden analysis:
        hvantk enrichex burden -m cohort.mt -p pheno.ht -s gene_sets.json -o results.tsv
    """
    ctx.ensure_object(dict)
    ctx.obj["log_level"] = log_level

    # Apply the selected log level to the root logger
    numeric_level = getattr(logging, log_level.upper(), logging.INFO)
    logging.getLogger().setLevel(numeric_level)


# Import and register subcommands
from hvantk.commands.enrichex_cli.burden_cli import register_burden_commands
from hvantk.commands.enrichex_cli.overlap_cli import register_overlap_commands

register_overlap_commands(enrichex_group)
register_burden_commands(enrichex_group)
