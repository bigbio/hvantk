"""CLI command and lifecycle entry point for downloading CPTAC protein-expression data.

Examples:
    hvantk download cptac-expression --cancer-type brca -o data/cptac/
    hvantk download cptac-expression --all -o data/cptac/
    hvantk download cptac-expression --list-cancers
"""

import logging
from typing import Optional

import click

from hvantk.skills.cptac.shared.constants import CPTAC_CANCER_TYPES

logger = logging.getLogger(__name__)


def download_dataset(
    raw_dir: str,
    cancer_type: Optional[str] = None,
    overwrite: bool = False,
    **kwargs,
) -> dict:
    """Lifecycle entry point for the plugin loader.

    Mirrors :func:`hvantk.skills.cptac.phospho.cli.download_dataset` -- see its
    docstring for the ``DatasetSpec`` contract. For CPTAC expression, "raw" means
    the per-cancer-type long-format matrix CSV plus clinical metadata CSV already
    produced by :class:`CPTACExpressionDataset` -- the upstream ``cptac`` Python
    package fetches into memory only.

    Parameters
    ----------
    raw_dir : str
        Directory under which the per-cancer-type files are placed.
    cancer_type : str, optional
        Single cancer type from :data:`CPTAC_CANCER_TYPES`. If omitted, all
        cancer types are downloaded.
    overwrite : bool
        If True, re-download even if files already exist.
    **kwargs
        Reserved for future lifecycle keyword arguments; ignored today.

    Returns
    -------
    dict
        Mapping of cancer-type to per-cancer output-file dict (the dict
        returned by :meth:`CPTACExpressionDataset.download`).
    """
    from hvantk.skills.cptac.shared.datasets import (
        CPTACExpressionDataset,
        download_cancer_types,
    )

    return download_cancer_types(
        CPTACExpressionDataset,
        raw_dir,
        cancer_type=cancer_type,
        overwrite=overwrite,
        label="CPTAC expression",
    )


@click.command(name="cptac-expression-download")
@click.option(
    "-o",
    "--output-dir",
    type=click.Path(),
    default=None,
    help="Directory to save CPTAC expression output files",
)
@click.option(
    "--cancer-type",
    type=click.Choice(CPTAC_CANCER_TYPES, case_sensitive=False),
    default=None,
    help="Cancer type to download",
)
@click.option(
    "--all",
    "download_all",
    is_flag=True,
    help="Download every cancer type (one matrix+metadata CSV pair each)",
)
@click.option(
    "--list-cancers",
    is_flag=True,
    help="List available cancer types and exit",
)
@click.option(
    "--overwrite",
    is_flag=True,
    help="Overwrite existing files",
)
@click.pass_context
def download_cmd(ctx, output_dir, cancer_type, download_all, list_cancers, overwrite):
    """Download CPTAC whole-proteome protein-expression data.

    Downloads protein abundance from the CPTAC Python package and writes a
    long-format matrix CSV (SampleID, GeneID, Gene Name, Expression) plus a
    clinical metadata CSV per cancer type. Unlike ``cptac-phospho``, ``--all``
    does not also produce a pan-cancer merge: each cancer type's matrix keeps
    its own gene set and sample columns, and build_cptac_expression consumes
    one cancer type's pair at a time.

    \b
    Requires: pip install cptac
    """
    if list_cancers:
        click.echo("Available CPTAC cancer types:")
        for ct in CPTAC_CANCER_TYPES:
            click.echo(f"  {ct}")
        return

    if not output_dir:
        click.echo(
            "Error: --output-dir is required (unless using --list-cancers)", err=True
        )
        ctx.exit(1)

    if not cancer_type and not download_all:
        click.echo(
            "Error: specify --cancer-type or --all. "
            "Use --list-cancers to see available types.",
            err=True,
        )
        ctx.exit(1)

    # hvantk class (no cptac import at module scope), so this never fails on a
    # missing cptac package -- that ImportError surfaces lazily inside
    # dataset.download() and is handled per-cancer below.
    from hvantk.skills.cptac.shared.datasets import CPTACExpressionDataset

    cancer_types = CPTAC_CANCER_TYPES if download_all else [cancer_type]

    succeeded, failed = [], {}
    for ct in cancer_types:
        click.echo(f"Processing {ct}...")
        try:
            paths = CPTACExpressionDataset(cancer_type=ct).download(
                output_dir, overwrite=overwrite
            )
        except ImportError as e:
            # The cptac package itself is missing -> nothing can succeed.
            click.echo(f"Error: {e}", err=True)
            ctx.exit(1)
        except Exception as exc:  # noqa: BLE001 - one bad cancer must not lose the rest
            logger.warning("CPTAC expression download failed for %s: %s", ct, exc)
            click.echo(f"  FAILED: {exc}", err=True)
            failed[ct] = str(exc)
            continue
        click.echo(f"  Matrix: {paths['matrix']}")
        click.echo(f"  Metadata: {paths['metadata']}")
        succeeded.append(ct)

    if failed:
        click.echo(
            f"Summary: {len(succeeded)} succeeded, {len(failed)} failed "
            f"({', '.join(sorted(failed))}).",
            err=True,
        )
    if not succeeded:
        ctx.exit(1)
