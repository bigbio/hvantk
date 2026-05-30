"""Builder for `onek-genomes:variants` — imports per-chromosome 1KG VCFs into a VariantMatrix."""
from __future__ import annotations

import glob
import logging
import os
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Optional

logger = logging.getLogger(__name__)

# Chromosome sort order: autosomes first, then sex chromosomes
_CHROM_ORDER = {str(i): i for i in range(1, 23)}
_CHROM_ORDER.update({"X": 23, "Y": 24, "M": 25, "MT": 26})

_CHROM_PATTERN = re.compile(r"(?i)chr[_]?(\d+|X|Y|M|MT)\b")
_STANDARD_CHROMOSOMES = [f"chr{c}" for c in [*range(1, 23), "X", "Y"]]


def _extract_chrom_token(filepath: str) -> Optional[str]:
    basename = os.path.basename(filepath)
    match = _CHROM_PATTERN.search(basename)
    if match:
        return "chr" + match.group(1).upper()
    return None


def _chrom_sort_key(filepath: str) -> int:
    token = _extract_chrom_token(filepath)
    if token is None:
        return 999
    chrom = token.replace("chr", "").upper()
    return _CHROM_ORDER.get(chrom, 998)


def _discover_vcf_files(
    input_dir: str,
    chromosomes: Optional[List[str]] = None,
) -> List[str]:
    """Discover bgzipped VCFs in *input_dir* and return them sorted by chromosome."""
    vcf_pattern = os.path.join(input_dir, "*.vcf.gz")
    vcf_files = sorted(glob.glob(vcf_pattern))

    if not vcf_files:
        raise FileNotFoundError(
            f"No *.vcf.gz files found in '{input_dir}'. "
            "Ensure files are bgzipped (.vcf.gz) and located directly in that directory."
        )

    if chromosomes is not None:
        requested = {c.upper() for c in chromosomes}
        vcf_files = [
            f for f in vcf_files if (_extract_chrom_token(f) or "").upper() in requested
        ]
        if not vcf_files:
            raise ValueError(
                f"No VCF files matched the requested chromosomes: {chromosomes}. "
                "Check that the filenames contain chromosome tokens like 'chr1', 'chrX'."
            )

    missing_indexes = [f for f in vcf_files if not os.path.exists(f + ".tbi")]
    if missing_indexes:
        raise FileNotFoundError(
            "The following VCF files are missing a .tbi tabix index:\n"
            + "\n".join(f"  {f}" for f in missing_indexes)
        )

    if chromosomes is None:
        found_tokens = {(_extract_chrom_token(f) or "").lower() for f in vcf_files}
        missing_chroms = [
            c for c in _STANDARD_CHROMOSOMES if c.lower() not in found_tokens
        ]
        if missing_chroms:
            logger.warning(
                "The following standard chromosomes were not found in '%s': %s",
                input_dir,
                ", ".join(missing_chroms),
            )

    return sorted(vcf_files, key=_chrom_sort_key)


def _resolve_compression_parallel(
    vcf_files: List[str], auto_convert_bgz: bool,
) -> tuple[List[str], bool]:
    """Run hvantk.core.utils.file_utils.resolve_compression in parallel across files."""
    from hvantk.core.utils.file_utils import resolve_compression

    n_files = len(vcf_files)
    max_workers = min(n_files, os.cpu_count() or 4)
    resolved_files: List[Optional[str]] = [None] * n_files
    force_bgz = True

    def _resolve_one(idx: int, vcf: str) -> tuple[int, str, bool]:
        logger.info(
            "Resolving compression for file %d/%d: %s",
            idx + 1, n_files, os.path.basename(vcf),
        )
        path, fbgz = resolve_compression(vcf, force_bgz=True, auto_convert=auto_convert_bgz)
        return idx, path, fbgz

    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = [pool.submit(_resolve_one, i, vcf) for i, vcf in enumerate(vcf_files)]
        for future in as_completed(futures):
            idx, resolved_path, file_force_bgz = future.result()
            resolved_files[idx] = resolved_path
            if not file_force_bgz:
                force_bgz = False

    return [f for f in resolved_files if f is not None], force_bgz


def build_onek_genomes_variants(parsed_input, ctx, **params):
    """Phase B builder — returns a VariantMatrix from a directory of 1KG VCFs.

    Parameters
    ----------
    parsed_input : str | Path
        Directory containing per-chromosome bgzipped 1KG VCFs (*.vcf.gz with .tbi).
    ctx : hvantk.core.models.BuildContext
        Platform-provided context; supplies provenance.
    **params
        chromosomes : list[str] | str — optional chromosome subset. Accepts a
            comma-separated string as a workaround for #119.
        reference_genome : str — defaults to "GRCh38".
        auto_convert_bgz : bool — defaults to False.
    """
    import hail as hl
    from hvantk.core.models import VariantMatrix

    input_vcfs = str(parsed_input)
    chromosomes = params.get("chromosomes")
    if isinstance(chromosomes, str):
        # Workaround for #119 — comma-list arg coercion. Remove this branch when #119 lands.
        chromosomes = [c.strip() for c in chromosomes.split(",") if c.strip()]
    reference_genome = params.get("reference_genome", "GRCh38")
    auto_convert_bgz = params.get("auto_convert_bgz", False)

    logger.info("Discovering VCF files in '%s'", input_vcfs)
    vcf_files = _discover_vcf_files(input_vcfs, chromosomes=chromosomes)
    logger.info("Found %d VCF file(s)", len(vcf_files))

    resolved_files, force_bgz = _resolve_compression_parallel(vcf_files, auto_convert_bgz)
    logger.info(
        "Importing %d VCF file(s) with reference genome '%s'",
        len(resolved_files), reference_genome,
    )
    mt = hl.import_vcf(
        resolved_files,
        force_bgz=force_bgz,
        reference_genome=reference_genome,
        array_elements_required=False,
    )
    return VariantMatrix.from_hail_mt(
        mt, provenance=ctx.provenance(schema_id="onek-genomes-variants-v1"),
    )
