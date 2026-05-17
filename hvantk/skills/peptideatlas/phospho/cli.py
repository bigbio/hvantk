"""CLI command and lifecycle entry point for downloading PeptideAtlas phospho builds.

Examples:
    # Default: download the latest known build
    hvantk peptideatlas-phospho-download -o /data/peptideatlas

    # Pin a specific build
    hvantk peptideatlas-phospho-download -o /data/peptideatlas --build-date 202512 --build-id 606
"""

import logging

import click

from hvantk.skills.peptideatlas.phospho.shared.datasets import PeptideAtlasPhosphoDataset

logger = logging.getLogger(__name__)


def download_dataset(
    raw_dir: str,
    build_date: str | None = None,
    build_id: str | None = None,
    overwrite: bool = False,
    **kwargs,
) -> str:
    """Lifecycle entry point for the plugin loader.

    Per the ``DatasetSpec`` contract documented in
    :mod:`hvantk.core.plugin_api`, a lifecycle ``download_fn`` accepts
    ``raw_dir=<path>`` and writes the raw upstream files (and, for
    PeptideAtlas, the parsed intermediate TSV) under that directory.

    Parameters
    ----------
    raw_dir : str
        Directory under which the zip + intermediate TSV are placed.
    build_date : str, optional
        Build date in ``YYYYMM`` format. Defaults to the latest known build.
    build_id : str, optional
        Build numeric ID. Defaults to the latest known build.
    overwrite : bool
        If True, re-download and re-parse even if files already exist.
    **kwargs
        Reserved for future lifecycle keyword arguments; ignored today.

    Returns
    -------
    str
        Path to the intermediate TSV file produced by the dataset.
    """
    if build_date and build_id:
        dataset = PeptideAtlasPhosphoDataset.from_build(build_date, build_id)
    else:
        dataset = PeptideAtlasPhosphoDataset.from_latest()
    return dataset.download(raw_dir, overwrite=overwrite)


@click.command(name="peptideatlas-phospho-download")
@click.option(
    "-o",
    "--output-dir",
    type=click.Path(),
    required=True,
    help="Directory to save the PeptideAtlas phospho TSV file",
)
@click.option(
    "--build-date",
    type=str,
    default=None,
    help="Build date in YYYYMM format (default: latest build)",
)
@click.option(
    "--build-id",
    type=str,
    default=None,
    help="Build ID number (default: latest build)",
)
@click.option(
    "--overwrite",
    is_flag=True,
    help="Overwrite existing file",
)
@click.pass_context
def download_cmd(ctx, output_dir, build_date, build_id, overwrite):
    """Download human phospho site data from PeptideAtlas.

    Downloads the PeptideAtlas phospho build TSV dump, parses phospho
    site positions with observation counts, and writes an intermediate
    TSV compatible with the PTM pipeline.
    """
    try:
        if build_date and build_id:
            dataset = PeptideAtlasPhosphoDataset.from_build(build_date, build_id)
        else:
            dataset = PeptideAtlasPhosphoDataset.from_latest()

        click.echo(f"Downloading PeptideAtlas phospho data ({dataset})...")
        output_path = dataset.download(output_dir, overwrite=overwrite)
        click.echo(f"Downloaded to: {output_path}")

    except (OSError, ValueError, RuntimeError) as e:
        logger.exception("Download failed: %s", e)
        click.echo(f"Error: {e}", err=True)
        ctx.exit(1)


# Backwards-compatible alias so existing imports in
# ``hvantk/commands/download_cli.py`` continue to work while still pointing at
# the new plugin module. New code should import ``download_cmd`` directly.
peptideatlas_phospho_downloader = download_cmd
