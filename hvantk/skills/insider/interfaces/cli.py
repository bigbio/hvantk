"""CLI command and lifecycle entry point for downloading the INSIDER per-protein
interface table (``H_sapiens_interfacesALL.txt``).

The sibling ``insider:variants`` dataset stays BYO (its BED is >1 GB), but this
~49 MB pair table is public, small, and served from a stable direct URL (see
``drift_probe.py``), so it qualifies for a real downloader per
``docs_site/guide/data-sources.md#downloader-or-documented-acquisition`` (#386).

Examples:
    # Download via the Click command
    hvantk download insider-interfaces --output data/insider/H_sapiens_interfacesALL.txt

    # Download with overwrite
    hvantk download insider-interfaces --output data/insider/H_sapiens_interfacesALL.txt --overwrite
"""

from __future__ import annotations

import logging
from pathlib import Path

import click
import requests

from hvantk.core.utils.http import stream_to_file
from hvantk.skills.insider.interfaces.builder import DEFAULT_FILENAME
from hvantk.skills.insider.interfaces.drift_probe import INSIDER_INTERFACES_URL

logger = logging.getLogger(__name__)

_TIMEOUT_S = 30


def download_insider_interfaces(output_path: str, overwrite: bool = False) -> str:
    """
    Download the INSIDER per-protein-pair interface table.

    Delegates the transfer to ``hvantk.core.utils.http.stream_to_file``, which
    writes through a sibling temp file and only ``os.replace``s it onto
    *output_path* once the whole body has arrived, so a failed or interrupted
    download never leaves a half-written file where the builder expects one.

    Parameters
    ----------
    output_path : str
        Path to save the downloaded file.
    overwrite : bool
        Whether to overwrite an existing file.

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
    output = Path(output_path)

    if output.exists() and not overwrite:
        raise FileExistsError(f"File already exists: {output_path}")

    try:
        logger.info(
            "Downloading INSIDER interfaces table from %s", INSIDER_INTERFACES_URL
        )
        stream_to_file(INSIDER_INTERFACES_URL, output, timeout=_TIMEOUT_S)
        logger.info("Downloaded to %s", output_path)
        return str(output)

    except requests.exceptions.RequestException as e:
        raise RuntimeError(f"Failed to download INSIDER interfaces table: {e}") from e
    except OSError as e:
        raise RuntimeError(f"Failed to write file: {e}") from e


def download_dataset(raw_dir: str, overwrite: bool = False, **kwargs) -> str:
    """Lifecycle entry point for the plugin loader.

    Per the ``DatasetSpec`` contract documented in ``hvantk.core.plugin.api``, a
    lifecycle ``download_fn`` accepts ``raw_dir=<path>`` and writes the raw
    upstream files under that directory. For ``insider:interfaces`` there is
    exactly one artifact, the per-protein-pair table, written as
    ``<raw_dir>/H_sapiens_interfacesALL.txt`` -- the filename
    ``builder._resolve_input`` looks for when handed a directory.

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
    target = raw / DEFAULT_FILENAME
    return download_insider_interfaces(str(target), overwrite=overwrite)


@click.command(
    "insider-interfaces-download",
    short_help="Download the INSIDER per-protein interface table",
)
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
    Download the INSIDER per-protein-pair interface table
    (H_sapiens_interfacesALL.txt).

    Per-protein-pair interface residue arrays with Source provenance (ECLAIR
    predicted / PDB / I3D homology) from the Yu lab's Interactome Insider
    resource. Backs ``insider:interfaces`` (~49 MB); the sibling
    ``insider:variants`` BED stays a manual download (>1 GB).

    For more information: http://interactomeinsider.yulab.org/downloads.html

    Examples:

        # Download the interfaces table

        hvantk download insider-interfaces --output data/insider/H_sapiens_interfacesALL.txt

        # Overwrite an existing download

        hvantk download insider-interfaces --output data/insider/H_sapiens_interfacesALL.txt --overwrite
    """
    click.echo("Downloading INSIDER per-protein interface table...")
    click.echo(f"Source: {INSIDER_INTERFACES_URL}")

    try:
        result_path = download_insider_interfaces(output_path, overwrite=overwrite)
        click.echo(f"Downloaded to: {result_path}")

    except FileExistsError as e:
        click.echo(f"Error: {e}", err=True)
        click.echo("Use --overwrite to replace existing file.")
        ctx.exit(1)
    except RuntimeError as e:
        click.echo(f"Download failed: {e}", err=True)
        ctx.exit(1)
