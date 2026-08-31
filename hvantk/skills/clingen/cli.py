"""
CLI command and lifecycle entry point for downloading ClinGen Gene-Disease
Validity data.

Examples:
    # Download today's snapshot
    hvantk clingen-download --output-dir data/clingen

    # Check download availability
    hvantk clingen-download --list-versions
"""

import logging

import click

from hvantk.skills.clingen.shared.datasets import ClinGenGeneDiseaseDataset

logger = logging.getLogger(__name__)


def _print_available_versions():
    """Print ClinGen Gene-Disease dataset availability info."""
    from hvantk.skills.clingen.shared.datasets import get_available_versions

    versions = get_available_versions()
    if not versions:
        click.echo("Unable to reach ClinGen endpoint. Check network connection.")
        return

    click.echo("ClinGen Gene-Disease Validity download is available.")
    click.echo(
        "Note: ClinGen now provides real-time snapshots (no versioned archives)."
    )
    click.echo(f"  Download will be labeled with today's date: {versions[0]}")


def download_dataset(raw_dir: str, overwrite: bool = False, **kwargs) -> str:
    """Lifecycle entry point for the plugin loader.

    Per the ``DatasetSpec`` contract documented in
    :mod:`hvantk.core.plugin_api`, a lifecycle ``download_fn`` accepts
    ``raw_dir=<path>`` and writes the raw upstream files under that
    directory. For ClinGen there is exactly one artifact, the real-time
    Gene-Disease Validity CSV, labelled with today's date.

    Parameters
    ----------
    raw_dir : str
        Directory under which the downloaded CSV is placed.
    overwrite : bool, optional
        Whether to overwrite an existing file (default: False).
    **kwargs
        Reserved for future lifecycle keyword arguments; ignored today.

    Returns
    -------
    str
        Path to the downloaded raw CSV file.
    """
    dataset = ClinGenGeneDiseaseDataset.from_latest()
    return dataset.download(output_dir=raw_dir, overwrite=overwrite)


@click.command(
    "clingen-download", short_help="Download ClinGen Gene-Disease Validity data"
)
@click.option(
    "--version",
    "version_date",
    type=str,
    default="latest",
    show_default=True,
    help="'latest' to download today's snapshot (only option; "
    "ClinGen does not provide archival versions).",
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
def download_cmd(ctx, version_date, output_dir, overwrite, list_versions):
    """
    Download ClinGen Gene-Disease Validity data.

    ClinGen provides curated gene-disease associations with evidence-based
    classifications (Definitive, Strong, Moderate, Limited, etc.).

    The downloaded CSV file contains gene-disease pairs with classification
    levels, mode of inheritance, and other metadata.

    Examples:

        # Download today's snapshot

        hvantk clingen-download --output-dir data/clingen

        # Check download availability

        hvantk clingen-download --list-versions
    """
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


# Backward-compatible alias for the legacy public name. The umbrella
# ``hvantk/tools/plugins/download_cli.py`` still imports the Click command under
# this name; once the auto-attach for ``cli:`` manifest blocks lands, the
# alias and the umbrella's import can be removed together.
clingen_downloader = download_cmd


if __name__ == "__main__":
    download_cmd()
