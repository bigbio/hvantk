"""
CLI command for downloading ClinGen Gene-Disease Validity data.

Examples:
    # Download latest version
    hvantk clingen-downloader --output-dir data/clingen

    # Download specific version
    hvantk clingen-downloader --version 2026-01-15 --output-dir data/clingen

    # List available versions
    hvantk clingen-downloader --list-versions
"""

import logging

import click

from hvantk.core.config import CONTEXT_SETTINGS

logger = logging.getLogger(__name__)


def _print_available_versions():
    """Print ClinGen Gene-Disease dataset availability info."""
    from hvantk.datasets.clingen_datasets import get_available_versions

    versions = get_available_versions()
    if not versions:
        click.echo("Unable to reach ClinGen endpoint. Check network connection.")
        return

    click.echo("ClinGen Gene-Disease Validity download is available.")
    click.echo(
        "Note: ClinGen now provides real-time snapshots (no versioned archives)."
    )
    click.echo(f"  Download will be labeled with today's date: {versions[0]}")


@click.command(
    "clingen-downloader", short_help="Download ClinGen Gene-Disease Validity data"
)
@click.option(
    "--version",
    "version_date",
    type=str,
    default="latest",
    show_default=True,
    help="Version date (YYYY-MM-DD) or 'latest' for most recent.",
)
@click.option(
    "--output-dir",
    type=str,
    default="data/clingen",
    show_default=True,
    help="Output directory for downloaded files.",
)
@click.option(
    "--overwrite",
    is_flag=True,
    help="Overwrite existing files if present.",
)
@click.option(
    "--list-versions",
    is_flag=True,
    help="Check download availability and show the current snapshot date.",
)
@click.pass_context
def clingen_downloader(ctx, version_date, output_dir, overwrite, list_versions):
    """
    Download ClinGen Gene-Disease Validity data.

    ClinGen provides curated gene-disease associations with evidence-based
    classifications (Definitive, Strong, Moderate, Limited, etc.).

    The downloaded CSV file contains gene-disease pairs with classification
    levels, mode of inheritance, and other metadata.

    Examples:

        # Download the latest version

        hvantk clingen-downloader --output-dir data/clingen

        # Download a specific version by date

        hvantk clingen-downloader --version 2026-01-15 --output-dir data/clingen

        # List available versions

        hvantk clingen-downloader --list-versions
    """
    from hvantk.datasets.clingen_datasets import ClinGenGeneDiseaseDataset

    if list_versions:
        _print_available_versions()
        ctx.exit(0)

    # Create dataset reference
    try:
        if version_date == "latest":
            click.echo("Fetching latest ClinGen Gene-Disease dataset...")
            dataset = ClinGenGeneDiseaseDataset.from_latest()
        else:
            dataset = ClinGenGeneDiseaseDataset.from_date(version_date)

        click.echo(f"Dataset version: {dataset.version_date}")
        click.echo(f"File: {dataset.file_name}")

    except ValueError as e:
        click.echo(f"Error: {e}", err=True)
        ctx.exit(1)
    except RuntimeError as e:
        click.echo(f"Error: {e}", err=True)
        ctx.exit(1)

    # Download the dataset
    try:
        output_path = dataset.download(output_dir=output_dir, overwrite=overwrite)
        click.echo(f"Downloaded to: {output_path}")
        logger.info(f"ClinGen dataset downloaded to {output_path}")

    except FileExistsError as e:
        click.echo(f"Error: {e}", err=True)
        click.echo("Use --overwrite to replace existing files.")
        ctx.exit(1)
    except RuntimeError as e:
        click.echo(f"Download failed: {e}", err=True)
        ctx.exit(1)


if __name__ == "__main__":
    clingen_downloader()
