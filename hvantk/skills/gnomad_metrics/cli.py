"""CLI command and lifecycle entry point for downloading gnomAD constraint
gene metrics.

The gnomAD constraint tables are small (v2.1.1 by-gene ~4.6 MB, v4.0
constraint metrics ~82 MB) and served from the public gnomAD GCS bucket, so
they qualify for a built-in downloader.

Examples:
    # Default: gnomAD v2.1.1 per-gene constraint (the hvantk standard)
    hvantk gnomad-metrics-download --output data/gnomad/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz

    # v4.0 (GRCh38) constraint metrics
    hvantk gnomad-metrics-download --version v4.0 --output data/gnomad/gnomad.v4.0.constraint_metrics.tsv
"""

import logging
from pathlib import Path

import click

from hvantk.core.config import (
    CONTEXT_SETTINGS,
)  # noqa: F401  (parity with sibling CLIs)
from hvantk.skills.gnomad_metrics.shared.constants import (
    DEFAULT_VERSION,
    GNOMAD_CONSTRAINT_TABLES,
    constraint_filename,
    constraint_url,
    resolve_table,
)

logger = logging.getLogger(__name__)


def download_gnomad_metrics(
    output_path: str,
    version: str = DEFAULT_VERSION,
    table: str | None = None,
    overwrite: bool = False,
) -> str:
    """Download a gnomAD constraint table to ``output_path``.

    Parameters
    ----------
    output_path : str
        Path to save the downloaded file.
    version : str
        gnomAD constraint release (``v2.1.1`` or ``v4.0``).
    table : str, optional
        Table within the release. Defaults to the per-release default
        (``by_gene`` for v2.1.1, ``constraint_metrics`` for v4.0).
    overwrite : bool
        Whether to overwrite an existing file.

    Returns
    -------
    str
        Path to the downloaded file.

    Raises
    ------
    ValueError
        If ``version``/``table`` is unknown.
    FileExistsError
        If the file exists and ``overwrite`` is False.
    RuntimeError
        If the download or write fails.
    """
    import urllib.error
    import urllib.request

    url = constraint_url(version, table)  # validates version/table first
    resolved_table = resolve_table(version, table)

    output = Path(output_path)
    if output.exists() and not overwrite:
        raise FileExistsError(f"File already exists: {output_path}")
    output.parent.mkdir(parents=True, exist_ok=True)

    try:
        logger.info(
            "Downloading gnomAD constraint (%s / %s) from %s",
            version,
            resolved_table,
            url,
        )
        urllib.request.urlretrieve(url, output)
        logger.info("Downloaded to %s", output_path)
        return str(output)
    except urllib.error.URLError as e:
        raise RuntimeError(f"Failed to download gnomAD constraint: {e}") from e
    except OSError as e:
        raise RuntimeError(f"Failed to write file: {e}") from e


def download_dataset(
    raw_dir: str,
    version: str = DEFAULT_VERSION,
    table: str | None = None,
    overwrite: bool = False,
    **kwargs,
) -> str:
    """Lifecycle entry point for the plugin loader.

    Per the ``DatasetSpec`` contract, a lifecycle ``download_fn`` accepts
    ``raw_dir=<path>`` and writes the raw upstream file under that directory,
    named with its canonical gnomAD basename.

    Parameters
    ----------
    raw_dir : str
        Directory under which the raw file is placed.
    version : str
        gnomAD constraint release (default ``v2.1.1``).
    table : str, optional
        Table within the release (default: per-release default).
    overwrite : bool
        Whether to overwrite an existing file.
    **kwargs
        Reserved for future lifecycle keyword arguments; ignored today.

    Returns
    -------
    str
        Path to the downloaded raw file.
    """
    raw = Path(raw_dir)
    raw.mkdir(parents=True, exist_ok=True)
    target = raw / constraint_filename(version, table)
    return download_gnomad_metrics(
        str(target), version=version, table=table, overwrite=overwrite
    )


@click.command(
    "gnomad-metrics-download",
    short_help="Download gnomAD constraint gene metrics (pLI, oe_lof/LOEUF, mis_z)",
)
@click.option(
    "--output",
    "output_path",
    type=str,
    required=True,
    help="Output path for the downloaded file.",
)
@click.option(
    "--version",
    type=click.Choice(sorted(GNOMAD_CONSTRAINT_TABLES)),
    default=DEFAULT_VERSION,
    show_default=True,
    help="gnomAD constraint release.",
)
@click.option(
    "--table",
    type=str,
    default=None,
    help=(
        "Table within the release. v2.1.1: by_gene (default), by_transcript. "
        "v4.0: constraint_metrics (default)."
    ),
)
@click.option(
    "--overwrite",
    is_flag=True,
    help="Overwrite existing file if present.",
)
@click.pass_context
def download_cmd(ctx, output_path, version, table, overwrite):
    """Download the gnomAD constraint gene-metrics table.

    gnomAD constraint provides per-gene intolerance scores (pLI, oe_lof and
    its upper bound LOEUF, missense z) widely used to rank genes by their
    tolerance to variation.

    v2.1.1 (GRCh37, ~4.6 MB) is the hvantk default and builds cleanly via
    ``hvantk reprocess gnomad-metrics:metrics`` (keyed by gene_id). v4.0
    (GRCh38, ~82 MB) has a different, transcript-level schema — see the SKILL.

    Upstream: https://gnomad.broadinstitute.org/downloads

    Examples:

        hvantk gnomad-metrics-download --output data/gnomad/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz

        hvantk gnomad-metrics-download --version v4.0 --output data/gnomad/gnomad.v4.0.constraint_metrics.tsv
    """
    try:
        url = constraint_url(version, table)
    except ValueError as e:
        click.echo(f"Error: {e}", err=True)
        ctx.exit(1)

    click.echo(f"Downloading gnomAD constraint ({version}) from {url} ...")
    try:
        result_path = download_gnomad_metrics(
            output_path, version=version, table=table, overwrite=overwrite
        )
        click.echo(f"Downloaded to: {result_path}")
    except FileExistsError as e:
        click.echo(f"Error: {e}", err=True)
        click.echo("Use --overwrite to replace existing file.")
        ctx.exit(1)
    except RuntimeError as e:
        click.echo(f"Download failed: {e}", err=True)
        ctx.exit(1)


# Backward-compatible alias mirroring sibling plugin CLIs.
gnomad_metrics_downloader = download_cmd


if __name__ == "__main__":
    download_cmd()
