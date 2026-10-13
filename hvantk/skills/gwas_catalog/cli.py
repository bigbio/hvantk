"""
CLI command and lifecycle entry point for downloading the EBI GWAS Catalog
v1.0 full-associations release.

Examples:
    # Download the current release
    hvantk download gwas-catalog --output data/gwas_catalog/gwas-catalog-download-associations-v1.0-full.tsv

    # Download with overwrite
    hvantk download gwas-catalog --output data/gwas_catalog/gwas-catalog-download-associations-v1.0-full.tsv --overwrite
"""

from __future__ import annotations

import logging
import os
import shutil
import tempfile
import zipfile
from pathlib import Path

import click
import requests

from hvantk.core.plugin.drift_runner import _current_umask
from hvantk.skills.gwas_catalog.drift_probe import GWAS_CATALOG_FULL_URL

logger = logging.getLogger(__name__)

# Canonical filename for the extracted associations TSV inside the raw dir. Matches the
# upstream member name recorded in hvantk/skills/gwas_catalog/catalog/datasets.json
# (files[0].path), so the raw directory mirrors what SKILL.md documents. Used both by
# the bare CLI (as a sensible default output name in the docstring/examples above) and
# by the lifecycle ``download_dataset`` wrapper, so the gevir-style directory resolver
# in builder.py (``_resolve_gwas_catalog_path``) finds exactly one file there.
GWAS_CATALOG_ASSOCIATIONS_FILENAME = "gwas-catalog-download-associations-v1.0-full.tsv"

_DOWNLOAD_TIMEOUT_S = 120
_CHUNK_SIZE = 1024 * 1024  # 1 MiB


def _stream_to_file(url: str, dest: Path, timeout: int = _DOWNLOAD_TIMEOUT_S) -> None:
    """Stream *url* to *dest* through a sibling temp file, then ``os.replace``.

    A plain ``open(dest, "wb")`` would leave a truncated file at *dest* if the
    network drops mid-download. Writing to a temp file in the same directory
    first and renaming it into place only on success keeps the write atomic --
    *dest* either does not exist or holds the complete download, never a
    partial one.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(dir=str(dest.parent), prefix=dest.name + ".")
    os.close(fd)
    done = False
    try:
        with requests.get(url, stream=True, timeout=timeout) as response:
            response.raise_for_status()
            with open(tmp_name, "wb") as fh:
                for chunk in response.iter_content(chunk_size=_CHUNK_SIZE):
                    if chunk:
                        fh.write(chunk)
        # mkstemp always creates mode 0600; match what a normal file create
        # (or the existing file, on overwrite) would get, mirroring
        # drift_runner.write_fingerprint's atomic-publish idiom.
        if dest.exists():
            shutil.copymode(dest, tmp_name)
        else:
            os.chmod(tmp_name, 0o666 & ~_current_umask())
        os.replace(tmp_name, dest)
        done = True
    finally:
        if not done and os.path.exists(tmp_name):
            os.remove(tmp_name)


def _extract_single_member(zip_path: Path, dest: Path) -> None:
    """Stream the one file inside *zip_path* out to *dest*.

    Uses ``ZipFile.open()`` + ``shutil.copyfileobj`` rather than
    ``ZipFile.extractall()``/``.read()`` so the member (the ~600 MB TSV) is
    copied in chunks and never fully materialized in memory. The release
    archive is expected to contain exactly one file; anything else is a shape
    hvantk does not know how to resolve, so it is a hard error rather than a
    silent guess. The target is written through the same temp-file +
    ``os.replace`` pattern as the zip download, for the same atomicity reason.
    """
    with zipfile.ZipFile(zip_path) as zf:
        members = [info for info in zf.infolist() if not info.is_dir()]
        if len(members) != 1:
            raise ValueError(
                f"expected exactly one file inside {zip_path}, found "
                f"{len(members)}: {[m.filename for m in members]}"
            )
        member = members[0]

        dest.parent.mkdir(parents=True, exist_ok=True)
        fd, tmp_name = tempfile.mkstemp(dir=str(dest.parent), prefix=dest.name + ".")
        os.close(fd)
        done = False
        try:
            with zf.open(member) as src, open(tmp_name, "wb") as out:
                shutil.copyfileobj(src, out, length=_CHUNK_SIZE)
            if dest.exists():
                shutil.copymode(dest, tmp_name)
            else:
                os.chmod(tmp_name, 0o666 & ~_current_umask())
            os.replace(tmp_name, dest)
            done = True
        finally:
            if not done and os.path.exists(tmp_name):
                os.remove(tmp_name)


def download_gwas_catalog_associations(
    output_path: str, overwrite: bool = False
) -> str:
    """
    Download and extract the current EBI GWAS Catalog full-associations TSV.

    Streams the current release zip from the EBI FTP ``releases/latest/``
    alias (``GWAS_CATALOG_FULL_URL``, ~70 MB) to a temp file, then streams the
    single TSV member it contains (~600 MB) out to *output_path*. Both writes
    are atomic (temp file + rename), so a failed or interrupted run never
    leaves a partial file at *output_path*. The downloaded zip is discarded
    once the TSV is extracted.

    Parameters
    ----------
    output_path : str
        Path to save the extracted associations TSV.
    overwrite : bool
        Whether to overwrite an existing file.

    Returns
    -------
    str
        Path to the extracted TSV.

    Raises
    ------
    FileExistsError
        If file exists and overwrite is False.
    RuntimeError
        If the download or extraction fails.
    """
    output = Path(output_path)

    if output.exists() and not overwrite:
        raise FileExistsError(f"File already exists: {output_path}")

    output.parent.mkdir(parents=True, exist_ok=True)

    try:
        with tempfile.TemporaryDirectory(prefix="gwas-catalog-dl-") as staging:
            zip_path = Path(staging) / "gwas-catalog-associations-full.zip"
            logger.info(
                "Downloading GWAS Catalog full-associations release from %s",
                GWAS_CATALOG_FULL_URL,
            )
            _stream_to_file(GWAS_CATALOG_FULL_URL, zip_path)
            logger.info("Extracting associations TSV to %s", output_path)
            _extract_single_member(zip_path, output)
    except requests.exceptions.RequestException as e:
        raise RuntimeError(f"Failed to download GWAS Catalog data: {e}") from e
    except (OSError, zipfile.BadZipFile, ValueError) as e:
        raise RuntimeError(f"Failed to extract GWAS Catalog data: {e}") from e

    logger.info("Downloaded to %s", output_path)
    return str(output)


def download_dataset(raw_dir: str, overwrite: bool = False, **kwargs) -> str:
    """Lifecycle entry point for the plugin loader.

    Per the ``DatasetSpec`` contract documented in
    :mod:`hvantk.core.plugin_api`, a lifecycle ``download_fn`` accepts
    ``raw_dir=<path>`` and writes the raw upstream files under that
    directory. For gwas-catalog there is exactly one artifact after
    extraction, the full-associations TSV, written as
    ``<raw_dir>/gwas-catalog-download-associations-v1.0-full.tsv`` -- the one
    file ``_resolve_gwas_catalog_path`` (``builder.py``) expects to find
    there.

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
        Path to the downloaded raw TSV file.
    """
    raw = Path(raw_dir)
    raw.mkdir(parents=True, exist_ok=True)
    target = raw / GWAS_CATALOG_ASSOCIATIONS_FILENAME
    return download_gwas_catalog_associations(str(target), overwrite=overwrite)


@click.command(
    "gwas-catalog-download",
    short_help="Download the EBI GWAS Catalog full-associations release",
)
@click.option(
    "--output",
    "output_path",
    type=str,
    required=True,
    help="Output path for the extracted associations TSV.",
)
@click.option(
    "--overwrite",
    is_flag=True,
    help="Overwrite existing file if present.",
)
@click.pass_context
def download_cmd(ctx, output_path, overwrite):
    """
    Download and extract the EBI GWAS Catalog v1.0 full-associations TSV.

    Fetches the current release zip from the EBI FTP ``releases/latest/``
    alias (~70 MB) and extracts the single TSV member it contains (~600 MB)
    to OUTPUT. A new release appears every 1-3 weeks; this always fetches
    whatever is current, not the release pinned in the plugin's catalog
    entry (``hvantk/skills/gwas_catalog/catalog/datasets.json``).

    For more information: https://www.ebi.ac.uk/gwas/

    Examples:

        # Download the current release

        hvantk download gwas-catalog --output data/gwas_catalog/gwas-catalog-download-associations-v1.0-full.tsv

        # Overwrite an existing file

        hvantk download gwas-catalog --output data/gwas_catalog/gwas-catalog-download-associations-v1.0-full.tsv --overwrite
    """
    click.echo("Downloading EBI GWAS Catalog full-associations release...")
    click.echo(f"Source: {GWAS_CATALOG_FULL_URL}")

    try:
        result_path = download_gwas_catalog_associations(
            output_path, overwrite=overwrite
        )
        click.echo(f"Downloaded to: {result_path}")

    except FileExistsError as e:
        click.echo(f"Error: {e}", err=True)
        click.echo("Use --overwrite to replace existing file.")
        ctx.exit(1)
    except RuntimeError as e:
        click.echo(f"Download failed: {e}", err=True)
        ctx.exit(1)


if __name__ == "__main__":
    download_cmd()
