"""
CLI command and lifecycle entry point for downloading UniProt PTM
(post-translational modification) data for reviewed human proteins.

Examples:
    # Download today's snapshot
    hvantk uniprot-ptm-download --output-dir data/uniprot_ptm

    # Download with an explicit date label
    hvantk uniprot-ptm-download --output-dir data/uniprot_ptm --version 2026-05-15
"""

import logging

import click

from hvantk.core.config import CONTEXT_SETTINGS  # noqa: F401  (kept for parity with sibling CLIs)
from hvantk.skills.uniprot_ptm.shared.datasets import UniProtPTMDataset

logger = logging.getLogger(__name__)


def download_dataset(raw_dir: str, overwrite: bool = False, **kwargs) -> str:
    """Lifecycle entry point for the plugin loader.

    Per the ``DatasetSpec`` contract documented in
    :mod:`hvantk.core.plugin_api`, a lifecycle ``download_fn`` accepts
    ``raw_dir=<path>`` and writes the raw upstream files under that
    directory. For UniProt PTM there is exactly one artifact, the live
    REST-API-derived MOD_RES TSV, labelled with today's date.

    Parameters
    ----------
    raw_dir : str
        Directory under which the downloaded TSV is placed.
    overwrite : bool, optional
        Whether to overwrite an existing file (default: False).
    **kwargs
        Reserved for future lifecycle keyword arguments; ignored today.

    Returns
    -------
    str
        Path to the downloaded raw TSV file.
    """
    dataset = UniProtPTMDataset.from_latest()
    return dataset.download(output_dir=raw_dir, overwrite=overwrite)


@click.command(
    "uniprot-ptm-download",
    short_help="Download UniProt PTM (MOD_RES) data for human proteins",
)
@click.option(
    "--version",
    "version_date",
    type=str,
    default="latest",
    show_default=True,
    help="Date label in YYYY-MM-DD format, or 'latest' for today's snapshot. "
    "UniProt provides a live search endpoint (no archival versions); the "
    "date is used only for labelling the output file.",
)
@click.option(
    "-o",
    "--output-dir",
    type=click.Path(),
    default="data/uniprot_ptm",
    show_default=True,
    help="Directory to save the downloaded PTM TSV file.",
)
@click.option(
    "--overwrite",
    is_flag=True,
    help="Overwrite existing files if present.",
)
@click.pass_context
def download_cmd(ctx, version_date, output_dir, overwrite):
    """
    Download human PTM annotations from the UniProt REST API.

    Queries UniProt for all reviewed human proteins with MOD_RES
    annotations and writes a TSV with one row per PTM site.

    Examples:

        # Download today's snapshot

        hvantk uniprot-ptm-download --output-dir data/uniprot_ptm

        # Download with an explicit date label

        hvantk uniprot-ptm-download --output-dir data/uniprot_ptm --version 2026-05-15
    """
    try:
        if version_date == "latest":
            click.echo("Fetching latest UniProt PTM dataset...")
            dataset = UniProtPTMDataset.from_latest()
        else:
            dataset = UniProtPTMDataset.from_date(version_date)

        click.echo(f"Dataset version: {dataset.version_date}")
        click.echo(f"File: {dataset.file_name}")

    except ValueError as e:
        click.echo(f"Error: {e}", err=True)
        ctx.exit(1)

    try:
        output_path = dataset.download(output_dir=output_dir, overwrite=overwrite)
        click.echo(f"Downloaded to: {output_path}")
        logger.info(f"UniProt PTM dataset downloaded to {output_path}")

    except FileExistsError as e:
        click.echo(f"Error: {e}", err=True)
        click.echo("Use --overwrite to replace existing files.")
        ctx.exit(1)
    except RuntimeError as e:
        click.echo(f"Download failed: {e}", err=True)
        ctx.exit(1)


# Backward-compatible alias for the legacy public name. The umbrella
# ``hvantk/commands/download_cli.py`` still imports the Click command under
# this name; once the auto-attach for ``cli:`` manifest blocks lands, the
# alias and the umbrella's import can be removed together.
uniprot_ptm_downloader = download_cmd


if __name__ == "__main__":
    download_cmd()
