"""
Standalone CLI command to convert standard gzip files to block gzip (BGZF).

Examples:
  hvantk convert-bgz input.tsv.gz
  hvantk convert-bgz input.tsv.gz -o output.tsv.bgz --threads 4
"""

import logging

import click

from hvantk.core.config import CONTEXT_SETTINGS

logger = logging.getLogger(__name__)


@click.command("convert-bgz", context_settings=CONTEXT_SETTINGS)
@click.argument("input_file", type=click.Path(exists=True))
@click.option(
    "-o",
    "--output",
    "output_file",
    default=None,
    type=click.Path(),
    help="Output .bgz file path (default: replaces .gz with .bgz)",
)
@click.option(
    "-t",
    "--threads",
    default=4,
    type=int,
    show_default=True,
    help="Number of threads for system bgzip (ignored by other backends)",
)
def convert_bgz_cmd(input_file: str, output_file: str, threads: int):
    """Convert a standard gzip (.gz) file to block gzip (.bgz) for Hail parallel import."""
    if threads < 1:
        raise click.BadParameter("must be a positive integer", param_hint="'--threads'")

    from hvantk.core.utils.file_utils import convert_gz_to_bgz, detect_compression

    compression = detect_compression(input_file)

    if compression == "bgzf":
        click.echo(f"File '{input_file}' is already BGZF compressed. Nothing to do.")
        return

    if compression != "gzip":
        click.echo(
            f"File '{input_file}' is not gzip compressed. "
            "This command only converts standard gzip to BGZF."
        )
        raise SystemExit(1)

    result_path = convert_gz_to_bgz(
        input_path=input_file,
        output_path=output_file,
        threads=threads,
    )
    click.echo(f"Converted to {result_path}")
