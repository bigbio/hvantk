"""
CLI command for downloading ClinVar VCF data.

Examples:
    # Download latest ClinVar VCF (GRCh38)
    hvantk clinvar-downloader --output-dir data/clinvar

    # Download a specific archived version
    hvantk clinvar-downloader --version 20260101 --output-dir data/clinvar

    # Download GRCh37 build without tabix index
    hvantk clinvar-downloader --genome-build GRCh37 --no-index

    # Download and verify MD5 checksum
    hvantk clinvar-downloader --verify-md5
"""

import logging

import click

from hvantk.core.config import CONTEXT_SETTINGS

logger = logging.getLogger(__name__)


@click.command("clinvar-downloader", short_help="Download ClinVar VCF data from NCBI")
@click.option(
    "--version",
    "version_date",
    type=str,
    default="latest",
    show_default=True,
    help="Version date (YYYYMMDD) or 'latest' for the current release.",
)
@click.option(
    "--output-dir",
    type=str,
    default="data/clinvar",
    show_default=True,
    help="Output directory for downloaded files.",
)
@click.option(
    "--genome-build",
    type=click.Choice(["GRCh38", "GRCh37"], case_sensitive=True),
    default="GRCh38",
    show_default=True,
    help="Reference genome build.",
)
@click.option(
    "--overwrite",
    is_flag=True,
    help="Overwrite existing files if present.",
)
@click.option(
    "--no-index",
    is_flag=True,
    help="Skip downloading the .tbi tabix index.",
)
@click.option(
    "--verify-md5",
    is_flag=True,
    help="Verify MD5 checksum after download.",
)
@click.pass_context
def clinvar_downloader(
    ctx, version_date, output_dir, genome_build, overwrite, no_index, verify_md5
):
    """
    Download ClinVar VCF data from NCBI FTP.

    ClinVar provides clinically relevant variant annotations including
    pathogenicity classifications, review status, and disease associations.

    Downloads the VCF file and optionally the tabix index (.tbi).

    Examples:

        # Download the latest ClinVar VCF

        hvantk clinvar-downloader --output-dir data/clinvar

        # Download a specific archived version

        hvantk clinvar-downloader --version 20260101 --output-dir data/clinvar

        # Download GRCh37 build

        hvantk clinvar-downloader --genome-build GRCh37
    """
    from hvantk.datasets.clinvar_datasets import ClinVarDataset

    # Create dataset reference
    try:
        if version_date == "latest":
            dataset = ClinVarDataset.latest(genome_build=genome_build)
            click.echo(f"Using latest ClinVar release ({genome_build})")
        else:
            dataset = ClinVarDataset.from_date(version_date, genome_build=genome_build)

        click.echo(f"File: {dataset.file_name}")
        click.echo(f"URL: {dataset.download_url}")

    except ValueError as e:
        click.echo(f"Error: {e}", err=True)
        ctx.exit(1)

    # Download
    try:
        output_path = dataset.download(
            output_dir=output_dir,
            overwrite=overwrite,
            download_index=not no_index,
        )
        click.echo(f"Downloaded to: {output_path}")

    except FileExistsError as e:
        click.echo(f"Error: {e}", err=True)
        click.echo("Use --overwrite to replace existing files.")
        ctx.exit(1)
    except RuntimeError as e:
        click.echo(f"Download failed: {e}", err=True)
        ctx.exit(1)

    # Optional MD5 verification
    if verify_md5:
        click.echo("Verifying MD5 checksum...")
        try:
            if dataset.verify_md5(output_path):
                click.echo("MD5 checksum verified successfully.")
            else:
                click.echo("WARNING: MD5 checksum mismatch!", err=True)
                ctx.exit(1)
        except RuntimeError as e:
            click.echo(f"MD5 verification failed: {e}", err=True)
            ctx.exit(1)


if __name__ == "__main__":
    clinvar_downloader()
