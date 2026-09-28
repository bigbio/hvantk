"""HGC CLI Commands Module - Main Entry Point.

Exposes a Click group whose subcommand modules (``combine_cli``, ``convert_cli``,
``pipeline_cli``, ``qc_cli``) are themselves cheap to import: none of them imports
``hvantk.algorithms.hgc`` (and therefore Hail) at module scope any more. Each command
imports the Hail-dependent implementation inside its own function body instead, so the
subcommands can be registered here at import time -- ``hvantk hgc --help`` and unrelated
subcommands no longer pay for Hail; the heavy imports only fire when the user actually
runs ``hvantk hgc ...``.
"""

import logging

import click

from .combine_cli import register_combine_commands
from .convert_cli import register_convert_commands
from .pipeline_cli import register_pipeline_command
from .qc_cli import register_qc_commands

logger = logging.getLogger(__name__)


@click.group(
    name="hgc",
    help="HGC (Hail-based Genotype Combiner) commands for joint genotyping workflows",
)
@click.pass_context
def hgc_group(ctx):
    """HGC command group for joint genotyping operations.

    Logging is configured centrally via ``hvantk -v`` / ``hvantk --log-file``.
    """
    ctx.ensure_object(dict)
    logger.info("Starting HGC command")


register_combine_commands(hgc_group)
register_convert_commands(hgc_group)
register_qc_commands(hgc_group)
register_pipeline_command(hgc_group)
