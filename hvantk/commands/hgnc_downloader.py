"""
CLI command for downloading HGNC gene nomenclature data.

Examples:
    # Download to default location
    hvantk hgnc-downloader --output data/hgnc/hgnc_complete_set.txt

    # Download with overwrite
    hvantk hgnc-downloader --output data/hgnc/hgnc_complete_set.txt --overwrite
"""

import logging
from pathlib import Path

import click

from hvantk.core.config import CONTEXT_SETTINGS
from hvantk.core.constants import HGNC_DOWNLOAD_URL, HGNC_INFO_URL

logger = logging.getLogger(__name__)


def download_hgnc(output_path: str, overwrite: bool = False) -> str:
    """
    Download the HGNC complete gene set.

    Parameters
    ----------
    output_path : str
        Path to save the downloaded file.
    overwrite : bool
        Whether to overwrite existing file.

    Returns
    -------
    str
        Path to the downloaded file.

    Raises
    ------
    FileExistsError
        If file exists and overwrite is False.
    RuntimeError
        If download fails.
    """
    import urllib.request
    import urllib.error

    output = Path(output_path)

    if output.exists() and not overwrite:
        raise FileExistsError(f"File already exists: {output_path}")

    # Create parent directories if needed
    output.parent.mkdir(parents=True, exist_ok=True)

    try:
        logger.info(f"Downloading HGNC data from {HGNC_DOWNLOAD_URL}")
        urllib.request.urlretrieve(HGNC_DOWNLOAD_URL, output)
        logger.info(f"Downloaded to {output_path}")
        return str(output)

    except urllib.error.URLError as e:
        raise RuntimeError(f"Failed to download HGNC data: {e}") from e
    except OSError as e:
        raise RuntimeError(f"Failed to write file: {e}") from e


@click.command("hgnc-downloader", short_help="Download HGNC gene nomenclature data")
@click.option(
    "--output",
    "output_path",
    type=str,
    required=True,
    help="Output path for the downloaded file.",
)
@click.option(
    "--overwrite",
    is_flag=True,
    help="Overwrite existing file if present.",
)
@click.pass_context
def hgnc_downloader(ctx, output_path, overwrite):
    """
    Download the HGNC complete gene nomenclature dataset.

    HGNC (HUGO Gene Nomenclature Committee) provides the authoritative
    source for human gene symbols and cross-references to other databases
    (Ensembl, Entrez, UniProt, OMIM, etc.).

    The downloaded file contains ~43,000 genes with identifiers, symbols,
    aliases, and cross-references.

    For more information: https://www.genenames.org/download/statistics-and-files/

    Examples:

        # Download HGNC data

        hvantk hgnc-downloader --output data/hgnc/hgnc_complete_set.txt

        # Overwrite existing file

        hvantk hgnc-downloader --output data/hgnc/hgnc_complete_set.txt --overwrite
    """
    click.echo("Downloading HGNC gene nomenclature data...")
    click.echo(f"Source: {HGNC_DOWNLOAD_URL}")

    try:
        result_path = download_hgnc(output_path, overwrite=overwrite)
        click.echo(f"Downloaded to: {result_path}")

    except FileExistsError as e:
        click.echo(f"Error: {e}", err=True)
        click.echo("Use --overwrite to replace existing file.")
        ctx.exit(1)
    except RuntimeError as e:
        click.echo(f"Download failed: {e}", err=True)
        ctx.exit(1)


if __name__ == "__main__":
    hgnc_downloader()
