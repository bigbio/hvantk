"""
CLI command for building a Hail MatrixTable from local 1000 Genomes VCF files.

Example usage::

    hvantk build-1k-genome \\
        --input-vcfs /data/1kg/vcfs/ \\
        --output-mt /data/1kg/1kg_genomes.mt

    # With sample annotations and chromosome subset
    hvantk build-1k-genome \\
        --input-vcfs /data/1kg/vcfs/ \\
        --output-mt /data/1kg/1kg_genomes.mt \\
        --sample-annotations /data/1kg/igsr_samples.tsv \\
        --chromosomes chr1,chr2,chr22,chrX \\
        --overwrite

    # With space-delimited PED file
    hvantk build-1k-genome \\
        --input-vcfs /data/1kg/vcfs/ \\
        --output-mt /data/1kg/1kg_genomes.mt \\
        --sample-annotations /data/1kg/samples.ped \\
        --sample-annotations-delimiter space
"""

import logging

import click

from hvantk.core.config import CONTEXT_SETTINGS

logger = logging.getLogger(__name__)


def _build_1k_genome_mt(**kwargs):
    """Lazy import wrapper — keeps Hail out of the import-time path."""
    from hvantk.tables.genome_builders import build_1k_genome_mt

    return build_1k_genome_mt(**kwargs)


@click.command("build-1k-genome", context_settings=CONTEXT_SETTINGS)
@click.option(
    "--input-vcfs",
    "input_vcfs",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True, resolve_path=True),
    help=(
        "Directory containing per-chromosome genotype VCF files (.vcf.gz). "
        "Each file must have a corresponding tabix index (.vcf.gz.tbi)."
    ),
)
@click.option(
    "--output-mt",
    "output_mt",
    required=True,
    type=click.Path(),
    help="Output path for the generated Hail MatrixTable (.mt).",
)
@click.option(
    "--sample-annotations",
    "sample_annotations",
    default=None,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, resolve_path=True),
    help=(
        "Optional delimited file with sample metadata/annotations (e.g. population, "
        "sex, family). Must contain a column whose values match the sample IDs in "
        "the VCFs. All columns are joined as a 'sample_annotations' struct on the "
        "MatrixTable columns."
    ),
)
@click.option(
    "--sample-annotations-delimiter",
    "sample_annotations_delimiter",
    default=None,
    type=str,
    help=(
        "Field delimiter for the sample annotations file. "
        "Accepts named aliases: 'space', 'tab', 'comma', 'semicolon', "
        "or any literal character. Defaults to tab."
    ),
)
@click.option(
    "--reference-genome",
    "reference_genome",
    default="GRCh38",
    show_default=True,
    type=click.Choice(["GRCh37", "GRCh38"], case_sensitive=False),
    help="Reference genome to use for VCF import.",
)
@click.option(
    "--chromosomes",
    "chromosomes",
    default=None,
    type=str,
    help=(
        "Comma-separated list of chromosomes to include "
        "(e.g. 'chr1,chr2,chrX'). "
        "When omitted, all discovered chromosomes are imported."
    ),
)
@click.option(
    "--overwrite",
    is_flag=True,
    default=False,
    help="Overwrite the output MatrixTable if it already exists.",
)
@click.option(
    "--auto-convert-bgz",
    "auto_convert_bgz",
    is_flag=True,
    default=False,
    help=(
        "Automatically convert plain-gzip VCF files to BGZF before import. "
        "Required when VCF files are gzip-compressed but not block-gzipped."
    ),
)
@click.option(
    "--force-reconvert-bgz",
    "force_reconvert_bgz",
    is_flag=True,
    default=False,
    help=(
        "Re-convert VCF files to BGZF even if they are already detected as "
        "block-gzipped. Use when files pass header checks but still cause "
        "ZipException errors in Hail."
    ),
)
def build_1k_genome_cmd(
    input_vcfs: str,
    output_mt: str,
    sample_annotations: str | None,
    sample_annotations_delimiter: str | None,
    reference_genome: str,
    chromosomes: str | None,
    overwrite: bool,
    auto_convert_bgz: bool,
    force_reconvert_bgz: bool,
) -> None:
    """Build a Hail MatrixTable from local 1000 Genomes high-coverage VCF files.

    \b
    Input requirements
    ------------------
    * VCF files must be bgzipped (.vcf.gz) with tabix indexes (.vcf.gz.tbi).
    * One file per chromosome is expected (chr1-chr22, chrX, chrY).
    * Files are auto-detected and sorted by chromosome order.

    \b
    Dataset source (reference only — download not included)
    -------------------------------------------------------
    https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/1000G_2504_high_coverage/
    """
    chrom_list = (
        [c.strip() for c in chromosomes.split(",") if c.strip()]
        if chromosomes
        else None
    )

    from hvantk.tables.genome_builders import resolve_delimiter

    logger.info("Starting 1000 Genomes MatrixTable build")
    mt = _build_1k_genome_mt(
        input_vcfs=input_vcfs,
        output_mt=output_mt,
        sample_annotations=sample_annotations,
        sample_annotations_delimiter=resolve_delimiter(sample_annotations_delimiter),
        reference_genome=reference_genome,
        chromosomes=chrom_list,
        overwrite=overwrite,
        auto_convert_bgz=auto_convert_bgz,
        force_reconvert_bgz=force_reconvert_bgz,
    )

    n_variants = mt.count_rows()
    n_samples = mt.count_cols()
    click.echo(f"1000 Genomes MatrixTable created at {output_mt}")
    click.echo(f"  Variants : {n_variants:,}")
    click.echo(f"  Samples  : {n_samples:,}")