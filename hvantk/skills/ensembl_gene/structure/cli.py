"""Download lifecycle for ensembl-gene:structure -- fetches the pinned Ensembl GTF."""

from __future__ import annotations

import logging
import os

import click

from hvantk.core.utils.http import stream_to_file
from hvantk.resources.ensembl_release import (
    ENSEMBL_GTF_FILENAME,
    ENSEMBL_GTF_URL,
)

logger = logging.getLogger(__name__)

# urlretrieve had no timeout; 30s matches this codebase's other short-request
# conventions (clinvar/clingen/gencc's version checks).
_DOWNLOAD_TIMEOUT_S = 30


def download_dataset(raw_dir: str, **params) -> str:
    """Fetch the pinned Ensembl GTF into ``raw_dir`` and return its path."""
    os.makedirs(raw_dir, exist_ok=True)
    target = os.path.join(raw_dir, ENSEMBL_GTF_FILENAME)

    if os.path.exists(target) and not params.get("overwrite", False):
        logger.info("Using cached GTF: %s", target)
        return target

    logger.info("Downloading %s -> %s", ENSEMBL_GTF_URL, target)
    stream_to_file(ENSEMBL_GTF_URL, target, timeout=_DOWNLOAD_TIMEOUT_S)
    return target


@click.command("ensembl-structure-download")
@click.option("--raw-dir", required=True, help="Directory to download the GTF into.")
@click.option("--overwrite", is_flag=True, help="Re-download even if cached.")
def download_cmd(raw_dir, overwrite):
    """Download the pinned Ensembl GTF."""
    path = download_dataset(raw_dir, overwrite=overwrite)
    click.echo(path)
