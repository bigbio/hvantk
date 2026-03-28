"""
HGC CLI Commands Module - Main Entry Point

This module provides the main HGC command group and registers all HGC subcommands.

Command Groups:
- Combine: combine GVCF files and VDS datasets
- Convert: convert between VDS, MatrixTable, and VCF formats
- QC: quality control metrics, filtering, and visualization
- Pipeline: end-to-end workflow orchestration
"""

import logging
import click

from .combine_cli import register_combine_commands
from .convert_cli import register_convert_commands
from .qc_cli import register_qc_commands
from .pipeline_cli import register_pipeline_command
logger = logging.getLogger(__name__)


@click.group(
    name="hgc",
    help="HGC (Hail-based Genotype Combiner) commands for joint genotyping workflows",
)
@click.pass_context
def hgc_group(ctx):
    """
    HGC command group for joint genotyping operations.

    Logging is configured centrally via ``hvantk -v`` / ``hvantk --log-file``.
    """
    ctx.ensure_object(dict)
    logger.info("Starting HGC command")


# Register all command groups
register_combine_commands(hgc_group)
register_convert_commands(hgc_group)
register_qc_commands(hgc_group)
register_pipeline_command(hgc_group)
