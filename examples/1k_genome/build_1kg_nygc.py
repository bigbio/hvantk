#!/usr/bin/env python3
"""Build a Hail MatrixTable from the NYGC/CCDG 1000 Genomes high-coverage VCF files.

The VCF directory should contain the per-chromosome phased genotype VCFs from the
high-coverage NYGC/CCDG release (see
`data-sources.md <https://bigbio.github.io/hvantk/guide/data-sources/#1000-genomes-nygcccdg-high-coverage-callset>`_
for the mirror). Autosomes (and most other contigs) ship as
``*_chr*.filtered.shapeit2-duohmm-phased.vcf.gz``; chrX is phased separately and
ships as ``*_chrX.filtered.eagle2-phased.v2.vcf.gz``. Both name forms are staged;
any ``*.annotated.*`` or ``*_others.*`` companion files present alongside them
are excluded.

Example usage::

    # Basic
    python build_1kg_nygc.py --vcf-dir /data/1kg/vcfs --output-mt /data/1kg/out.mt

    # With chromosomes and sample annotations
    python build_1kg_nygc.py \\
        --vcf-dir /data/1kg/vcfs \\
        --output-mt /data/1kg/out.mt \\
        --chromosomes chr1,chr2,chrX \\
        --sample-annotations samples.ped \\
        --sample-annotations-delimiter space \\
        --overwrite
"""

from __future__ import annotations

import glob
import logging
import os
import re
import shutil
import tempfile
from datetime import datetime

import click

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
)
logger = logging.getLogger(__name__)


# Two name forms ship in the NYGC/CCDG release: autosomes (and most other
# contigs) as shapeit2-duohmm-phased, chrX separately as eagle2-phased.
_VCF_GLOB_PATTERNS = (
    "*_chr*.filtered.shapeit2-duohmm-phased.vcf.gz",
    "*_chrX.filtered.eagle2-phased*.vcf.gz",
)

# Mirrors the builder's own chromosome-token regex (hvantk/skills/onek_genomes/
# builder.py's _CHROM_PATTERN) so "is chrX staged?" agrees with how the plugin
# itself will read these symlinked files.
_CHROM_TOKEN_RE = re.compile(r"(?i)chr[_]?(\d+|X|Y|M|MT)\b")


def _chrom_token(name: str) -> str | None:
    """Extract the chromosome token (e.g. "1", "X") from a filename or a
    --chromosomes entry, so staged files and requests compare equal."""
    match = _CHROM_TOKEN_RE.search(name)
    return match.group(1).upper() if match else None


def _stage_vcfs(vcf_dir: str, stage_dir: str) -> tuple[int, set[str]]:
    """Symlink phased genotype VCFs into *stage_dir*.

    Excludes ``*.annotated.*`` and ``*_others.*`` files, if present.
    Returns the number of staged VCFs and the set of chromosome tokens found.
    """
    vcfs: set[str] = set()
    for pattern in _VCF_GLOB_PATTERNS:
        vcfs.update(glob.glob(os.path.join(vcf_dir, pattern)))

    count = 0
    staged_chroms: set[str] = set()
    for vcf in sorted(vcfs):
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
        token = _chrom_token(basename)
        if token:
            staged_chroms.add(token)
    return count, staged_chroms


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
    help=(
        "Field delimiter for the annotations file. "
        "Accepts named aliases: 'space', 'tab', 'comma', 'semicolon', "
        "or any literal character. Defaults to tab."
    ),
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
@click.option(
    "--auto-convert-bgz",
    is_flag=True,
    default=False,
    help=(
        "Convert all gzip-family VCF files (including ones detected as BGZF) "
        "to a clean BGZF file before import. Use this when VCF files cause "
        "ZipException errors in Hail, even if they appear to be valid BGZF."
    ),
)
def main(
    vcf_dir: str,
    output_mt: str,
    chromosomes: str | None,
    sample_annotations: str | None,
    sample_annotations_delimiter: str | None,
    reference_genome: str,
    overwrite: bool,
    auto_convert_bgz: bool,
) -> None:
    """Build a Hail MatrixTable from NYGC 1000 Genomes VCF files.

    \b
    Stages phased genotype VCFs (*_chr*.filtered.shapeit2-duohmm-phased.vcf.gz,
    plus chrX's separately-named *_chrX.filtered.eagle2-phased.v2.vcf.gz) into a
    temporary directory, excluding any annotated or "others" contig files, then
    invokes the hvantk builder.
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

    chrom_list = None
    if chromosomes:
        chrom_list = [c.strip() for c in chromosomes.split(",") if c.strip()]

    # -- Banner --
    logger.info("=" * 42)
    logger.info("  1000G NYGC MatrixTable Build")
    logger.info("  Started : %s", datetime.now().isoformat(timespec="seconds"))
    logger.info("  VCF dir : %s", vcf_dir)
    logger.info("  Output  : %s", output_mt)
    logger.info("=" * 42)

    # -- Stage VCFs --
    stage_dir = tempfile.mkdtemp(prefix="1kg_stage_")
    try:
        logger.info("Staging phased genotype VCFs into %s ...", stage_dir)
        count, staged_chroms = _stage_vcfs(vcf_dir, stage_dir)
        logger.info("Staged %d VCF file(s).", count)

        if count == 0:
            raise click.ClickException(f"No phased genotype VCFs found in {vcf_dir}")

        if chrom_list:
            missing = sorted(
                {c for c in chrom_list if _chrom_token(c) not in staged_chroms}
            )
            if missing:
                raise click.ClickException(
                    "Requested --chromosomes not found among staged VCFs in "
                    f"{vcf_dir}: {', '.join(missing)}. Check the file naming "
                    "against the module docstring (chrX ships under a "
                    "different name than the autosomes)."
                )

        # -- Build via hvantk reprocess onek-genomes:variants --
        # The legacy `build_1k_genome_mt` Python helper was retired in #116;
        # all builds now route through the plugin system.
        import subprocess

        cmd = [
            "hvantk",
            "reprocess",
            "onek-genomes:variants",
            "--raw-dir",
            stage_dir,
            "--output",
            output_mt,
            "--skip-download",
            "--plugin-arg",
            f"reference_genome={reference_genome}",
            "--plugin-arg",
            f"auto_convert_bgz={'true' if auto_convert_bgz else 'false'}",
        ]
        if chrom_list:
            cmd += ["--plugin-arg", f"chromosomes={','.join(chrom_list)}"]
        if sample_annotations_path is not None:
            logger.warning(
                "--sample-annotations is no longer baked into the variants build "
                "(post-#116). Run `hvantk reprocess onek-genomes:samples` to "
                "build the canonical IGSR samples table, and join post-load. "
                "Provided sample_annotations argument (%s) will be ignored.",
                sample_annotations_path,
            )

        logger.info("Invoking: %s", " ".join(cmd))
        subprocess.run(cmd, check=True)

        logger.info("=" * 42)
        logger.info("  Completed : %s", datetime.now().isoformat(timespec="seconds"))
        logger.info("  MatrixTable : %s", output_mt)
        logger.info("=" * 42)
    finally:
        shutil.rmtree(stage_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
