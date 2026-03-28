"""
CLI command to download UniProt PTM data for human proteins.

Example:
    hvantk download uniprot-ptm --output-dir /data/ptm/
"""

import logging

import click

logger = logging.getLogger(__name__)


@click.command(name="uniprot-ptm")
@click.option(
    "-o",
    "--output-dir",
    type=click.Path(),
    required=True,
    help="Directory to save the downloaded PTM TSV file",
)
@click.option(
    "--date",
    type=str,
    default=None,
    help="Version date label in YYYY-MM-DD format (default: latest UniProt snapshot)",
)
@click.option(
    "--overwrite",
    is_flag=True,
    help="Overwrite existing file",
)
@click.pass_context
def uniprot_ptm_downloader(ctx, output_dir, date, overwrite):
    """Download human PTM annotations from UniProt REST API.

    Queries UniProt for all reviewed human proteins with MOD_RES
    annotations and writes a TSV with one row per PTM site.
    """
    try:
        from hvantk.datasets.uniprot_ptm_datasets import UniProtPTMDataset

        if date:
            dataset = UniProtPTMDataset.from_date(date)
        else:
            dataset = UniProtPTMDataset.from_latest()

        click.echo(f"Downloading UniProt PTM data ({dataset})...")
        output_path = dataset.download(output_dir, overwrite=overwrite)
        click.echo(f"Downloaded to: {output_path}")

    except Exception as e:
        logger.exception(f"Download failed: {e}")
        click.echo(f"Error: {e}", err=True)
        ctx.exit(1)
