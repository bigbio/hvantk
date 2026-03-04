"""
MatrixTable builder for the 1000 Genomes Project high-coverage genotype VCFs.

Discovers per-chromosome VCF files in a local directory, imports them as a single
Hail MatrixTable, and optionally joins sample annotations/metadata.

Expected input layout:
    /path/to/vcfs/
        20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.filtered.shapeit2-duohmm-phased.vcf.gz
        20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.filtered.shapeit2-duohmm-phased.vcf.gz.tbi
        ...
        20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.shapeit2-duohmm-phased.vcf.gz
        20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.shapeit2-duohmm-phased.vcf.gz.tbi

Any filename containing a recognized chromosome token (chr1-chr22, chrX, chrY) is
accepted, regardless of naming convention.
"""

from __future__ import annotations

import glob
import logging
import os
import re
from typing import List, Optional

logger = logging.getLogger(__name__)

__all__ = [
    "build_1k_genome_mt",
    "discover_vcf_files",
    "resolve_delimiter",
    "DELIMITER_ALIASES",
]

# Named delimiter aliases — used by CLIs to avoid ambiguous quoted strings.
DELIMITER_ALIASES: dict[str, str] = {
    "space": " ",
    "tab": "\t",
    "comma": ",",
    "semicolon": ";",
}


def resolve_delimiter(value: Optional[str]) -> Optional[str]:
    """Resolve a delimiter value that may be a named alias or a literal character.

    Recognised aliases: ``space``, ``tab``, ``comma``, ``semicolon``.
    Any other string is returned unchanged (allows literal characters).
    *None* passes through as *None*.
    """
    if value is None:
        return None
    return DELIMITER_ALIASES.get(value.lower(), value)

# Chromosome sort order: autosomes first, then sex chromosomes
_CHROM_ORDER = {str(i): i for i in range(1, 23)}
_CHROM_ORDER.update({"X": 23, "Y": 24, "M": 25, "MT": 26})

# Recognise tokens like chr1, Chr22, CHR_X, chrX, chrY, chrM, chrMT
_CHROM_PATTERN = re.compile(r"(?i)chr[_]?(\d+|X|Y|M|MT)\b")

_STANDARD_CHROMOSOMES = [f"chr{c}" for c in [*range(1, 23), "X", "Y"]]


def _extract_chrom_token(filepath: str) -> Optional[str]:
    """Return the first 'chrN' token found in the filename, or None."""
    basename = os.path.basename(filepath)
    match = _CHROM_PATTERN.search(basename)
    if match:
        return "chr" + match.group(1).upper()
    return None


def _chrom_sort_key(filepath: str) -> int:
    """Numeric sort key derived from the chromosome token."""
    token = _extract_chrom_token(filepath)
    if token is None:
        return 999
    chrom = token.replace("chr", "").upper()
    return _CHROM_ORDER.get(chrom, 998)


def discover_vcf_files(
    input_dir: str,
    chromosomes: Optional[List[str]] = None,
) -> List[str]:
    """Discover bgzipped VCF files in *input_dir* and return them sorted by chromosome.

    Parameters
    ----------
    input_dir:
        Directory that contains ``*.vcf.gz`` files.
    chromosomes:
        Optional list of chromosome names to restrict to (e.g. ``["chr1", "chrX"]``).
        Names are matched case-insensitively against the token extracted from each
        filename.  If *None*, all detected chromosome files are returned.

    Returns
    -------
    list of str
        Absolute paths to VCF files, sorted in chromosome order.

    Raises
    ------
    FileNotFoundError
        If no ``*.vcf.gz`` files are found, or if any required ``.tbi`` index is
        missing.
    ValueError
        If a chromosome filter is requested but none of the discovered files match.
    """
    vcf_pattern = os.path.join(input_dir, "*.vcf.gz")
    vcf_files = sorted(glob.glob(vcf_pattern))

    if not vcf_files:
        raise FileNotFoundError(
            f"No *.vcf.gz files found in '{input_dir}'. "
            "Ensure files are bgzipped (.vcf.gz) and located directly in that directory."
        )

    # Apply chromosome filter
    if chromosomes is not None:
        requested = {c.upper() for c in chromosomes}
        vcf_files = [
            f for f in vcf_files
            if (_extract_chrom_token(f) or "").upper() in requested
        ]
        if not vcf_files:
            raise ValueError(
                f"No VCF files matched the requested chromosomes: {chromosomes}. "
                "Check that the filenames contain chromosome tokens like 'chr1', 'chrX'."
            )

    # Validate .tbi indexes on the (possibly filtered) file list
    missing_indexes = [f for f in vcf_files if not os.path.exists(f + ".tbi")]
    if missing_indexes:
        raise FileNotFoundError(
            "The following VCF files are missing a .tbi tabix index:\n"
            + "\n".join(f"  {f}" for f in missing_indexes)
        )

    # Warn about missing standard chromosomes only when no explicit filter was applied
    if chromosomes is None:
        found_tokens = {(_extract_chrom_token(f) or "").lower() for f in vcf_files}
        missing_chroms = [c for c in _STANDARD_CHROMOSOMES if c.lower() not in found_tokens]
        if missing_chroms:
            logger.warning(
                "The following standard chromosomes were not found in '%s': %s",
                input_dir,
                ", ".join(missing_chroms),
            )

    vcf_files = sorted(vcf_files, key=_chrom_sort_key)
    return vcf_files


def _join_sample_annotations(
    mt,
    annotations_path: str,
    sample_id_col: Optional[str] = None,
    delimiter: Optional[str] = None,
) -> "hl.MatrixTable":  # noqa: F821
    """Import a sample-annotations file and left-join it to the MatrixTable columns.

    The annotations table is keyed by *sample_id_col*. Each column (sample) in *mt*
    receives a ``sample_annotations`` struct annotation containing all fields.
    Samples absent from the annotations file get a ``null`` struct.

    Parameters
    ----------
    mt:
        Input MatrixTable whose column key is ``s`` (sample ID).
    annotations_path:
        Path to a delimited file with a sample ID column.
    sample_id_col:
        Name of the sample ID column in the annotations file.  When *None*, the
        function tries common names (``sample_id``, ``sample``, ``s``, ``ID``); if
        none match, the first column is used and a warning is emitted.
    delimiter:
        Field delimiter for the annotations file.  When *None*, Hail's default
        tab delimiter is used.

    Returns
    -------
    hl.MatrixTable
        The input MatrixTable with an added ``sample_annotations`` column annotation.
    """
    import hail as hl

    logger.info("Reading sample annotations file from %s", annotations_path)
    import_kwargs = {"impute": True}
    if delimiter is not None:
        import_kwargs["delimiter"] = delimiter
    annot_ht = hl.import_table(annotations_path, **import_kwargs)

    # Detect sample ID column
    available_cols = list(annot_ht.row.dtype)
    if sample_id_col is not None:
        if sample_id_col not in available_cols:
            raise ValueError(
                f"Sample ID column '{sample_id_col}' not found in annotations file. "
                f"Available columns: {available_cols}"
            )
    else:
        candidates = ["sample_id", "sample", "s", "ID", "SampleID", "id"]
        sample_id_col = next((c for c in candidates if c in available_cols), None)
        if sample_id_col is None:
            sample_id_col = available_cols[0]
            logger.warning(
                "No standard sample ID column found in annotations file. "
                "Using the first column: '%s'",
                sample_id_col,
            )
        else:
            logger.info("Using '%s' as sample ID column in annotations file", sample_id_col)

    annot_ht = annot_ht.key_by(sample_id_col)

    # Check overlap between samples
    mt_samples = set(mt.s.collect())
    annot_samples = set(annot_ht[sample_id_col].collect())
    overlap = mt_samples & annot_samples
    if not overlap:
        logger.warning(
            "No samples overlap between the VCF (n=%d) and the annotations file (n=%d). "
            "All sample annotations will be null.",
            len(mt_samples),
            len(annot_samples),
        )
    elif len(overlap) < len(mt_samples):
        missing_n = len(mt_samples) - len(overlap)
        logger.warning(
            "%d sample(s) in the VCF have no matching entry in the annotations file "
            "and will receive null annotations.",
            missing_n,
        )

    mt = mt.annotate_cols(sample_annotations=annot_ht[mt.s])
    return mt


def build_1k_genome_mt(
    input_vcfs: str,
    output_mt: str,
    sample_annotations: Optional[str] = None,
    reference_genome: str = "GRCh38",
    chromosomes: Optional[List[str]] = None,
    overwrite: bool = False,
    sample_id_col: Optional[str] = None,
    auto_convert_bgz: bool = False,
    sample_annotations_delimiter: Optional[str] = None,
) -> "hl.MatrixTable":  # noqa: F821
    """Build a Hail MatrixTable from local 1000 Genomes high-coverage VCF files.

    Discovers per-chromosome bgzipped VCF files in *input_vcfs*, imports them as a
    single MatrixTable, and checkpoints the result to *output_mt*.

    Parameters
    ----------
    input_vcfs:
        Directory containing genotype VCF files (``*.vcf.gz``) with tabix indexes
        (``*.vcf.gz.tbi``).  One file per chromosome is expected.
    output_mt:
        Output path for the generated Hail MatrixTable.
    sample_annotations:
        Optional path to a delimited file with sample metadata / annotations.  The
        file must contain a column whose values match the sample IDs in the VCF.
        All columns are joined as a ``sample_annotations`` struct on the MatrixTable
        columns.
    reference_genome:
        Reference genome for VCF import.  Defaults to ``"GRCh38"``.
    chromosomes:
        Optional list of chromosome names (e.g. ``["chr1", "chr2", "chrX"]``) to
        restrict the import.  When *None*, all discovered chromosomes are included.
    overwrite:
        If *True*, overwrite *output_mt* if it already exists.
    sample_id_col:
        Name of the sample ID column in *sample_annotations* file.  When *None*,
        common column names are tried automatically (see
        :func:`_join_sample_annotations`).
    auto_convert_bgz:
        If *True*, automatically convert plain gzip VCF files to BGZF before
        import.  Default is *False*.
    sample_annotations_delimiter:
        Field delimiter for the sample annotations file.  When *None*, Hail's
        default tab delimiter is used.  Use ``" "`` for space-delimited files
        (e.g. PED format).

    Returns
    -------
    hl.MatrixTable
        The imported and (optionally annotated) MatrixTable, checkpointed
        to *output_mt*.

    Raises
    ------
    FileNotFoundError
        If no VCF files or tabix indexes are found.
    ValueError
        If the chromosome filter matches no files.
    """
    import hail as hl

    from hvantk.data.file_utils import resolve_compression

    # --- Step 1: Discover input files ---
    logger.info("Discovering VCF files in '%s'", input_vcfs)
    vcf_files = discover_vcf_files(input_vcfs, chromosomes=chromosomes)
    logger.info(
        "Found %d VCF file(s): %s",
        len(vcf_files),
        [os.path.basename(f) for f in vcf_files],
    )

    # --- Step 2: Resolve compression for each VCF file ---
    resolved_files = []
    force_bgz = True
    for vcf in vcf_files:
        resolved_path, file_force_bgz = resolve_compression(
            vcf, force_bgz=True, auto_convert=auto_convert_bgz
        )
        resolved_files.append(resolved_path)
        if not file_force_bgz:
            force_bgz = False

    # --- Step 3: Import VCFs ---
    logger.info(
        "Importing %d VCF file(s) with reference genome '%s'",
        len(resolved_files),
        reference_genome,
    )
    mt = hl.import_vcf(
        resolved_files,
        force_bgz=force_bgz,
        reference_genome=reference_genome,
        array_elements_required=False,
    )

    # --- Step 3: Optional sample annotations join ---
    if sample_annotations is not None:
        mt = _join_sample_annotations(
            mt,
            sample_annotations,
            sample_id_col=sample_id_col,
            delimiter=sample_annotations_delimiter,
        )

    # --- Step 4: Checkpoint ---
    logger.info("Writing MatrixTable to '%s'", output_mt)
    mt = mt.checkpoint(output_mt, overwrite=overwrite)

    return mt