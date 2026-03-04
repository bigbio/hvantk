#!/usr/bin/env python3
"""End-to-end CLI wrapper for building a 1000 Genomes MatrixTable via ``hvantk``.

Unlike ``build_1kg_nygc.py`` (which calls the Python API directly), this script
invokes ``hvantk build-1k-genome`` as a subprocess so the full CLI path is
exercised — useful for catching argument parsing bugs, entry-point issues, etc.

Example usage::

    # Basic
    python build_1kg_nygc_cli.py --vcf-dir /data/1kg/vcfs --output-mt /data/1kg/out.mt

    # With chromosomes and sample annotations
    python build_1kg_nygc_cli.py \\
        --vcf-dir /data/1kg/vcfs \\
        --output-mt /data/1kg/out.mt \\
        --chromosomes chr1,chr2,chrX \\
        --sample-annotations samples.ped \\
        --sample-annotations-delimiter " " \\
        --overwrite
"""

from __future__ import annotations

import glob
import logging
import os
import shutil
import subprocess
import tempfile
from datetime import datetime

import click

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
)
logger = logging.getLogger(__name__)


def _stage_vcfs(vcf_dir: str, stage_dir: str) -> int:
    """Symlink recalibrated genotype VCFs into *stage_dir*.

    Excludes ``*.annotated.*`` and ``*_others.*`` files.
    Returns the number of staged VCFs.
    """
    pattern = os.path.join(vcf_dir, "*_chr*.recalibrated_variants.vcf.gz")
    count = 0
    for vcf in sorted(glob.glob(pattern)):
        basename = os.path.basename(vcf)
        if ".annotated." in basename or "_others." in basename:
            continue

        os.symlink(os.path.realpath(vcf), os.path.join(stage_dir, basename))
        tbi = vcf + ".tbi"
        if os.path.isfile(tbi):
            os.symlink(
                os.path.realpath(tbi),
                os.path.join(stage_dir, basename + ".tbi"),
            )
        else:
            logger.warning("Missing tabix index for %s", basename)
        count += 1
    return count


@click.command()
@click.option(
    "--vcf-dir",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True, resolve_path=True),
    help="Directory containing per-chromosome genotype VCF files.",
)
@click.option(
    "--output-mt",
    required=True,
    type=click.Path(),
    help="Output path for the generated Hail MatrixTable.",
)
@click.option(
    "--chromosomes",
    default=None,
    type=str,
    help="Comma-separated list of chromosomes to include (e.g. chr1,chr2,chrX).",
)
@click.option(
    "--sample-annotations",
    default=None,
    type=str,
    help=(
        "Path to a sample annotations file. Can be an absolute path or a "
        "filename relative to --vcf-dir (e.g. samples.tsv, samples.ped)."
    ),
)
@click.option(
    "--sample-annotations-delimiter",
    default=None,
    type=str,
    help='Field delimiter for the annotations file (default: tab). Use " " for PED.',
)
@click.option(
    "--reference-genome",
    default="GRCh38",
    show_default=True,
    type=click.Choice(["GRCh37", "GRCh38"], case_sensitive=False),
    help="Reference genome for VCF import.",
)
@click.option(
    "--overwrite",
    is_flag=True,
    default=False,
    help="Overwrite the output MatrixTable if it already exists.",
)
def main(
    vcf_dir: str,
    output_mt: str,
    chromosomes: str | None,
    sample_annotations: str | None,
    sample_annotations_delimiter: str | None,
    reference_genome: str,
    overwrite: bool,
) -> None:
    """Build a Hail MatrixTable from NYGC 1000 Genomes VCFs via the hvantk CLI.

    \b
    Stages only recalibrated genotype VCFs (*_chr*.recalibrated_variants.vcf.gz)
    into a temporary directory, excluding annotated and "others" contig files,
    then invokes ``hvantk build-1k-genome`` as a subprocess.
    """
    # Resolve sample annotations: absolute path used as-is, otherwise relative to vcf_dir
    sample_annotations_path = None
    if sample_annotations is not None:
        if os.path.isabs(sample_annotations):
            sample_annotations_path = sample_annotations
        else:
            sample_annotations_path = os.path.join(vcf_dir, sample_annotations)
        if not os.path.isfile(sample_annotations_path):
            raise click.BadParameter(
                f"File not found: {sample_annotations_path}",
                param_hint="'--sample-annotations'",
            )

    # -- Banner --
    logger.info("=" * 42)
    logger.info("  1000G NYGC MatrixTable Build (CLI)")
    logger.info("  Started : %s", datetime.now().isoformat(timespec="seconds"))
    logger.info("  VCF dir : %s", vcf_dir)
    logger.info("  Output  : %s", output_mt)
    logger.info("=" * 42)

    # -- Stage VCFs --
    stage_dir = tempfile.mkdtemp(prefix="1kg_stage_")
    try:
        logger.info("Staging recalibrated genotype VCFs into %s ...", stage_dir)
        count = _stage_vcfs(vcf_dir, stage_dir)
        logger.info("Staged %d VCF file(s).", count)

        if count == 0:
            raise click.ClickException(
                f"No recalibrated genotype VCFs found in {vcf_dir}"
            )

        # -- Build CLI command --
        cmd = [
            "hvantk", "build-1k-genome",
            "--input-vcfs", stage_dir,
            "--output-mt", output_mt,
            "--reference-genome", reference_genome,
        ]
        if chromosomes:
            cmd += ["--chromosomes", chromosomes]
        if sample_annotations_path:
            cmd += ["--sample-annotations", sample_annotations_path]
        if sample_annotations_delimiter:
            cmd += ["--sample-annotations-delimiter", sample_annotations_delimiter]
        if overwrite:
            cmd.append("--overwrite")

        logger.info("Running: %s", " ".join(cmd))
        result = subprocess.run(cmd, check=False)

        if result.returncode != 0:
            raise click.ClickException(
                f"hvantk build-1k-genome exited with code {result.returncode}"
            )

        logger.info("=" * 42)
        logger.info("  Completed : %s", datetime.now().isoformat(timespec="seconds"))
        logger.info("  MatrixTable : %s", output_mt)
        logger.info("=" * 42)
    finally:
        shutil.rmtree(stage_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
