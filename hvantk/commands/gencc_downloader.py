"""
CLI command for downloading GenCC (Gene Curation Coalition) submissions data.

Examples:
    # Download today's snapshot
    hvantk gencc-downloader --output-dir data/gencc

    # Check download availability
    hvantk gencc-downloader --list-versions
"""

import logging

import click

from hvantk.core.config import CONTEXT_SETTINGS

logger = logging.getLogger(__name__)


def _print_available_versions():
    """Print GenCC submissions dataset availability info."""
    from hvantk.datasets.gencc_datasets import get_available_versions

    versions = get_available_versions()
    if not versions:
        click.echo("Unable to reach GenCC endpoint. Check network connection.")
        return

    click.echo("GenCC submissions download is available.")
    click.echo("Note: GenCC provides real-time snapshots (no versioned archives).")
    click.echo(f"  Download will be labeled with today's date: {versions[0]}")


@click.command("gencc-downloader", short_help="Download GenCC submissions data")
@click.option(
    "--version",
    "version_date",
    type=str,
    default="latest",
    show_default=True,
    help="'latest' to download today's snapshot (only option; "
    "GenCC does not provide archival versions).",
)
@click.option(
    "--output-dir",
    type=str,
    default="data/gencc",
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
def gencc_downloader(ctx, version_date, output_dir, overwrite, list_versions):
    """
    Download GenCC (Gene Curation Coalition) submissions data.

    GenCC aggregates gene-disease validity assertions from 12+ submitting
    organizations (ClinGen, PanelApp, G2P, Orphanet, etc.) using HGNC
    gene IDs and MONDO disease IDs.

    The downloaded TSV file contains gene-disease-submitter assertions
    with classification levels and modes of inheritance.

    Examples:

        # Download today's snapshot

        hvantk gencc-downloader --output-dir data/gencc

        # Check download availability

        hvantk gencc-downloader --list-versions
    """
    from hvantk.datasets.gencc_datasets import GenCCSubmissionsDataset

    if list_versions:
        _print_available_versions()
        ctx.exit(0)

    try:
        if version_date == "latest":
            click.echo("Fetching latest GenCC submissions dataset...")
            dataset = GenCCSubmissionsDataset.from_latest()
        else:
            dataset = GenCCSubmissionsDataset.from_date(version_date)

        click.echo(f"Dataset version: {dataset.version_date}")
        click.echo(f"File: {dataset.file_name}")

    except (ValueError, RuntimeError) as e:
        click.echo(f"Error: {e}", err=True)
        ctx.exit(1)

    try:
        output_path = dataset.download(output_dir=output_dir, overwrite=overwrite)
        click.echo(f"Downloaded to: {output_path}")
        logger.info(f"GenCC dataset downloaded to {output_path}")

    except FileExistsError as e:
        click.echo(f"Error: {e}", err=True)
        click.echo("Use --overwrite to replace existing files.")
        ctx.exit(1)
    except RuntimeError as e:
        click.echo(f"Download failed: {e}", err=True)
        ctx.exit(1)


if __name__ == "__main__":
    gencc_downloader()
