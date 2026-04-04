"""
CLI command to download PeptideAtlas human phospho build data.

Example:
    hvantk download peptideatlas-phospho --output-dir /data/peptideatlas/
"""

import logging

import click
from hvantk.datasets.peptideatlas_phospho_datasets import PeptideAtlasPhosphoDataset

logger = logging.getLogger(__name__)


@click.command(name="peptideatlas-phospho")
@click.option(
    "-o",
    "--output-dir",
    type=click.Path(),
    required=True,
    help="Directory to save the PeptideAtlas phospho TSV file",
)
@click.option(
    "--build-date",
    type=str,
    default=None,
    help="Build date in YYYYMM format (default: latest build)",
)
@click.option(
    "--build-id",
    type=str,
    default=None,
    help="Build ID number (default: latest build)",
)
@click.option(
    "--overwrite",
    is_flag=True,
    help="Overwrite existing file",
)
@click.pass_context
def peptideatlas_phospho_downloader(ctx, output_dir, build_date, build_id, overwrite):
    """Download human phospho site data from PeptideAtlas.

    Downloads the PeptideAtlas phospho build TSV dump, parses phospho
    site positions with observation counts, and writes an intermediate
    TSV compatible with the PTM pipeline.
    """
    try:
        if build_date and build_id:
            dataset = PeptideAtlasPhosphoDataset.from_build(build_date, build_id)
        else:
            dataset = PeptideAtlasPhosphoDataset.from_latest()

        click.echo(f"Downloading PeptideAtlas phospho data ({dataset})...")
        output_path = dataset.download(output_dir, overwrite=overwrite)
        click.echo(f"Downloaded to: {output_path}")

    except (OSError, ValueError, RuntimeError) as e:
        logger.exception("Download failed: %s", e)
        click.echo(f"Error: {e}", err=True)
        ctx.exit(1)
