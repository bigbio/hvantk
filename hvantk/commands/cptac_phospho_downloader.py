"""
CLI command to download CPTAC phosphoproteomics data.

Example:
    hvantk download cptac-phospho --cancer-type brca -o data/cptac/
    hvantk download cptac-phospho --all -o data/cptac/
    hvantk download cptac-phospho --list-cancers
"""

import logging
import os

import click

from hvantk.ptm.constants import CPTAC_CANCER_TYPES

logger = logging.getLogger(__name__)


@click.command(name="cptac-phospho")
@click.option(
    "-o",
    "--output-dir",
    type=click.Path(),
    default=None,
    help="Directory to save CPTAC phospho output files",
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
    help="Download all cancer types and produce pan-cancer output",
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
def cptac_phospho_downloader(ctx, output_dir, cancer_type, download_all, list_cancers, overwrite):
    """Download CPTAC phosphoproteomics data.

    Downloads phospho site intensities from the CPTAC Python package,
    extracts site positions, and writes files for PTM pipeline and
    MatrixTable construction.

    \b
    Requires: pip install cptac
    """
    if list_cancers:
        click.echo("Available CPTAC cancer types:")
        for ct in CPTAC_CANCER_TYPES:
            click.echo(f"  {ct}")
        return

    if not output_dir:
        click.echo("Error: --output-dir is required (unless using --list-cancers)", err=True)
        ctx.exit(1)

    if not cancer_type and not download_all:
        click.echo(
            "Error: specify --cancer-type or --all. "
            "Use --list-cancers to see available types.",
            err=True,
        )
        ctx.exit(1)

    try:
        import csv
        from hvantk.datasets.cptac_phospho_datasets import CPTACPhosphoDataset, _TSV_COLUMNS

        cancer_types = CPTAC_CANCER_TYPES if download_all else [cancer_type]

        all_tsv_paths = []
        for ct in cancer_types:
            click.echo(f"Processing {ct}...")
            dataset = CPTACPhosphoDataset(cancer_type=ct)
            paths = dataset.download(output_dir, overwrite=overwrite)
            click.echo(f"  TSV: {paths['tsv']}")
            click.echo(f"  Matrix: {paths['matrix']}")
            click.echo(f"  Metadata: {paths['metadata']}")
            all_tsv_paths.append(paths["tsv"])

        # Pan-cancer merge if --all
        if download_all and len(all_tsv_paths) > 1:
            pancancer_path = os.path.join(output_dir, "cptac-phospho-pancancer.tsv")
            click.echo(f"Merging {len(all_tsv_paths)} cancer types into {pancancer_path}...")

            with open(pancancer_path, "w", newline="") as fout:
                writer = csv.DictWriter(
                    fout, fieldnames=_TSV_COLUMNS, delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                for tsv_path in all_tsv_paths:
                    with open(tsv_path) as fin:
                        reader = csv.DictReader(fin, delimiter="\t")
                        for row in reader:
                            writer.writerow(row)

            click.echo(f"Pan-cancer TSV: {pancancer_path}")

    except ImportError as e:
        click.echo(f"Error: {e}", err=True)
        ctx.exit(1)
    except Exception as e:
        logger.exception(f"Download failed: {e}")
        click.echo(f"Error: {e}", err=True)
        ctx.exit(1)
