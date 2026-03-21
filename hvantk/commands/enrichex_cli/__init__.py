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
@click.pass_context
def enrichex_group(ctx):
    """EnrichEx command group for gene set enrichment analysis.

    Logging is configured centrally via ``hvantk -v`` / ``hvantk --log-file``.

    Examples
    --------
    Test gene list enrichment:
        hvantk enrichex overlap -g genes.txt -s gene_sets.json -o results.tsv

    Run burden analysis:
        hvantk enrichex burden -m cohort.mt -p pheno.ht -s gene_sets.json -o results.tsv
    """
    ctx.ensure_object(dict)


# Import and register subcommands
from hvantk.commands.enrichex_cli.burden_cli import register_burden_commands
from hvantk.commands.enrichex_cli.overlap_cli import register_overlap_commands

register_overlap_commands(enrichex_group)
register_burden_commands(enrichex_group)
