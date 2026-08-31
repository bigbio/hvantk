"""
CLI command and lifecycle entry point for downloading HGNC gene nomenclature data.

Examples:
    # Download via the Click command
    hvantk hgnc-download --output data/hgnc/hgnc_complete_set.txt

    # Download with overwrite
    hvantk hgnc-download --output data/hgnc/hgnc_complete_set.txt --overwrite
"""

import logging
from pathlib import Path

import click

from hvantk.skills.hgnc.shared.constants import HGNC_DOWNLOAD_URL, HGNC_INFO_URL  # noqa: F401

logger = logging.getLogger(__name__)

# Canonical filename for the HGNC complete-set TSV. Used both by the bare CLI
# (as a sensible default file name) and by the lifecycle ``download_dataset``
# wrapper to know where to drop the raw file inside ``raw_dir``.
HGNC_COMPLETE_SET_FILENAME = "hgnc_complete_set.txt"


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


def download_dataset(raw_dir: str, overwrite: bool = False, **kwargs) -> str:
    """Lifecycle entry point for the plugin loader.

    Per the ``DatasetSpec`` contract documented in ``hvantk.core.plugin_api``,
    a lifecycle ``download_fn`` accepts ``raw_dir=<path>`` and writes the raw
    upstream files under that directory. For HGNC there is exactly one
    artifact, the complete-set TSV, written as
    ``<raw_dir>/hgnc_complete_set.txt``.

    Parameters
    ----------
    raw_dir : str
        Directory under which raw upstream files are placed.
    overwrite : bool, optional
        Whether to overwrite an existing file (default: False).
    **kwargs
        Reserved for future lifecycle keyword arguments; ignored today.

    Returns
    -------
    str
        Path to the downloaded raw file.
    """
    raw = Path(raw_dir)
    raw.mkdir(parents=True, exist_ok=True)
    target = raw / HGNC_COMPLETE_SET_FILENAME
    return download_hgnc(str(target), overwrite=overwrite)


@click.command("hgnc-download", short_help="Download HGNC gene nomenclature data")
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
def download_cmd(ctx, output_path, overwrite):
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

        hvantk hgnc-download --output data/hgnc/hgnc_complete_set.txt

        # Overwrite existing file

        hvantk hgnc-download --output data/hgnc/hgnc_complete_set.txt --overwrite
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


# Backward-compatible alias for the legacy public name. The umbrella
# ``hvantk/tools/plugins/download_cli.py`` still imports the Click command under
# this name; once the auto-attach for ``cli:`` manifest blocks lands, the
# alias and the umbrella's import can be removed together.
hgnc_downloader = download_cmd


if __name__ == "__main__":
    download_cmd()
