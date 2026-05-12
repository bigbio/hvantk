"""
Hail Table builders for converting raw data sources into Hail Tables.

This module provides builder functions that convert raw annotation sources
(TSV, VCF, BED) into checkpointed Hail Tables with standardized patterns.
"""

import hail as hl
import logging
import os
import re
from typing import Optional, List, Callable
from hvantk.utils.table_utils import get_row_fields, build_rename_map, str_to_bool
from hvantk.core.metadata import build_table_metadata

logger = logging.getLogger(__name__)
_FILE_URI_PREFIX = "file://"

from hvantk.core.constants import (
    ENSEMBL_BIOMART_FIELDS,
    CLINGEN_GENE_DISEASE_FIELDS,
    CLINGEN_CLASSIFICATION_LEVELS,
    GENCC_SUBMISSION_FIELDS,
    GENCC_CLASSIFICATION_LEVELS,
    HGNC_GENE_FIELDS,
    HGNC_PIPE_SEPARATED_FIELDS,
    COSMIC_CGC_FIELDS,
    COSMIC_CGC_CLASSIFICATION_LEVELS,
    COSMIC_MUTATION_CONTEXTS,
)
from hvantk.data.file_utils import resolve_compression
from hvantk.utils.genome import contig_recoding  # correct module import


def _create_table_base(
    source_name: str,
    input_path: str,
    output_path: str,
    import_func: Callable[[], hl.Table],
    transform_func: Optional[Callable[[hl.Table], hl.Table]] = None,
    cleanup_func: Optional[Callable[[], None]] = None,
    fields: Optional[List[str]] = None,
    overwrite: bool = False,
    export_tsv: bool = False,
) -> hl.Table:
    """
    Base helper for creating Hail Tables with common import/checkpoint/export pattern.

    This internal helper reduces boilerplate by handling the common workflow:
    1. Log import message
    2. Execute import function
    3. Apply optional transformations
    4. Apply optional field selection
    5. Annotate globals with provenance metadata
    6. Checkpoint with logging
    7. Optional TSV export

    Parameters
    ----------
    source_name : str
        Human-readable name of the data source (for logging).
    input_path : str
        Path to the input file.
    output_path : str
        Path to write the output Hail Table.
    import_func : Callable
        Function that returns the imported Hail Table.
    transform_func : Callable, optional
        Function to transform the table after import (default: None).
    cleanup_func : Callable, optional
        Cleanup hook run after checkpoint/export finish, even on failure
        (default: None).
    fields : list of str, optional
        List of fields to select from the table (default: None, keeps all).
    overwrite : bool, optional
        Whether to overwrite the output file if it exists (default: False).
    export_tsv : bool, optional
        If True, also export a TSV version (default: False).

    Returns
    -------
    hl.Table
        The checkpointed Hail Table.
    """
    try:
        logger.info(f"Creating {source_name} table from {input_path}")
        ht = import_func()

        if transform_func is not None:
            ht = transform_func(ht)

        if fields is not None:
            logger.info(f"Selecting fields: {fields}")
            ht = ht.select(*fields)

        ht = ht.annotate_globals(
            hvantk_metadata=build_table_metadata(source_name, input_path, ht)
        )

        logger.info(f"Checkpointing table to {output_path}")
        ht = ht.checkpoint(output=output_path, overwrite=overwrite)

        if export_tsv:
            logger.info(f"Exporting table to {output_path}.tsv.bgz")
            ht.export(output_path + ".tsv.bgz")

        return ht
    finally:
        if cleanup_func is not None:
            cleanup_func()


__all__ = [
    "create_gnomad_constraint_gene_metrics_tb",
    "create_interactome_tb",
    "create_clinvar_tb",
    "create_gevir_tb",
    "create_ensembl_gene_tb",
    "create_dbnsfp_tb",
    "create_clingen_gene_disease_tb",
    "create_hgnc_gene_tb",
    "create_ptm_sites_tb",
    "create_gwas_catalog_tb",
    "create_msigdb_tb",
]


def _cleanup_temp_file(tmp_path: Optional[str]) -> None:
    """Best-effort cleanup for local or Hadoop/S3/GS temp files."""
    if not tmp_path:
        return
    try:
        import hailtop.fs as hfs

        if hfs.exists(tmp_path):
            if hfs.is_dir(tmp_path):
                hfs.rmtree(tmp_path)
            else:
                hfs.remove(tmp_path)
        return
    except Exception:
        logger.debug(
            "Failed to remove temp path via hailtop.fs: %s", tmp_path, exc_info=True
        )

    try:
        local_path = tmp_path
        if local_path.startswith(_FILE_URI_PREFIX):
            local_path = local_path[len(_FILE_URI_PREFIX):]
        if os.path.exists(local_path):
            os.remove(local_path)
    except Exception:
        logger.debug(
            "Failed to remove temp path via os.remove: %s", tmp_path, exc_info=True
        )


def create_gnomad_constraint_gene_metrics_tb(
    input_path: str,
    output_path: str,
    fields: list = None,
    overwrite: bool = False,
    export_tsv: bool = False,
) -> "hl.Table":
    """
    Create a Hail Table from gnomAD constraint gene metrics TSV file keyed by gene_id.

    Example usage:
        ht = create_gnomad_constraint_gene_metrics_tb(
            input_path="/path/to/gnomad_constraint_metrics.txt",
            output_path="/path/to/output.ht",
            fields=["gene_id", "pLI", "oe_lof"]
        )

    Parameters
    ----------
    input_path : str
        Path to the gnomAD constraint gene metrics TSV input file.
    output_path : str
        Path to write the output Hail Table.
    fields : list, optional
        List of fields to select from the imported table (default: None, selects all fields).
    overwrite : bool, optional
        Whether to overwrite the output file if it exists (default: False).
    export_tsv : bool, optional
        If True, also export a TSV version (default: False).

    Returns
    -------
    hl.Table
        Hail Table keyed by gene_id with selected constraint metrics.

    Notes
    -----
    If schema is stable, consider specifying types=... instead of impute=True.
    """
    return _create_table_base(
        source_name="gnomAD constraint gene metrics",
        input_path=input_path,
        output_path=output_path,
        import_func=lambda: hl.import_table(
            paths=input_path, impute=True, min_partitions=100, key="gene_id"
        ),
        fields=fields,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )


_TRACK_NAME_RE = re.compile(r'name=([^\s]+)')


def _normalize_hadoop_path(path: str) -> str:
    """Normalize local file URIs for filesystem APIs."""
    return path[len(_FILE_URI_PREFIX):] if path.startswith(_FILE_URI_PREFIX) else path


def _parse_insider_bed_to_temp_tsv(input_path: str) -> str:
    """Pre-process an Interactome Insider BED into a TSV with ppi_id column.

    The INSIDER BED is segmented by `track name=<P1>_ppi_<P2> ...` directives,
    each followed by per-residue BED data rows. `hl.import_bed` silently skips
    the track headers, dropping PPI identity. This helper iterates the BED
    line-by-line, tracks the current PPI from each `track name=...` header, and
    writes a 4-column TSV (`contig\\tstart\\tend\\tppi_id`) for downstream
    `hl.import_table`.

    Filters applied:
      - `browser` lines are ignored.
      - Track headers with no parseable `name=...` are skipped (current_ppi_id
        becomes None, so subsequent rows until the next valid track are dropped).
      - Zero-length intervals (`start == end`) are dropped — matches the prior
        builder's `skip_invalid_intervals=True` behavior.

    Returns the path to a Hail-managed temp file (extension `tsv`).
    """
    import hailtop.fs as hfs

    out_path = hl.utils.new_temp_file(extension="tsv")
    current_ppi_id: Optional[str] = None
    n_rows_written = 0
    with hfs.open(_normalize_hadoop_path(input_path), "r") as src:
        with hfs.open(_normalize_hadoop_path(out_path), "w") as dst:
            dst.write("contig\tstart\tend\tppi_id\n")
            for line in src:
                if line.startswith("browser"):
                    continue
                if line.startswith("track"):
                    match = _TRACK_NAME_RE.search(line)
                    current_ppi_id = match.group(1) if match else None
                    continue
                if current_ppi_id is None:
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                try:
                    start = int(fields[1])
                    end = int(fields[2])
                except ValueError:
                    continue
                if start >= end:
                    continue
                dst.write(f"{fields[0]}\t{start}\t{end}\t{current_ppi_id}\n")
                n_rows_written += 1
    logger.info(
        "Parsed INSIDER BED %s into %s (%d data rows after filtering)",
        input_path,
        out_path,
        n_rows_written,
    )
    return out_path


def create_interactome_tb(
    input_path: str,
    output_path: str,
    overwrite: bool = False,
    export_tsv: bool = False,
    reference_genome: str = "GRCh38",
) -> "hl.Table":
    """
    Create a Hail Table from an Interactome Insider per-residue BED file.

    The BED is segmented by `track name=<P1>_ppi_<P2>` directives; this builder
    parses those headers and preserves PPI identity as a `ppi_ids: array<str>`
    field per interval. Intervals appearing in multiple PPI tracks are
    aggregated (collected as a sorted, deduplicated array).

    Example usage:
        ht = create_interactome_tb(
            input_path="/path/to/Whole_Human_Interactome_Interface_hg38.bed",
            output_path="/path/to/output.ht"
        )

    Parameters
    ----------
    input_path : str
        Path to the INSIDER BED input file (must contain `track name=...`
        directives to identify PPIs; plain BEDs without tracks produce
        empty output).
    output_path : str
        Path to write the output Hail Table.
    overwrite : bool, optional
        Whether to overwrite the output file if it exists (default: False).
    export_tsv : bool, optional
        If True, also export a TSV version (default: False).
    reference_genome : str, optional
        Reference genome to use for parsing intervals (default: "GRCh38").

    Returns
    -------
    hl.Table
        Hail Table keyed by `interval<locus<rg>>` with field
        `ppi_ids: array<str>` carrying the PPI identifiers from track headers.
    """
    tsv_path = None

    def _import() -> "hl.Table":
        nonlocal tsv_path
        tsv_path = _parse_insider_bed_to_temp_tsv(input_path)
        ht = hl.import_table(
            tsv_path,
            types={"start": hl.tint32, "end": hl.tint32},
            min_partitions=4,
        )
        # BED is 0-based half-open; Hail loci are 1-based. Match hl.import_bed's
        # conversion by shifting both endpoints by +1 (so a BED row [100, 200)
        # becomes Hail interval [chr:101, chr:201)).
        ht = ht.annotate(
            interval=hl.locus_interval(
                ht.contig,
                ht.start + 1,
                ht.end + 1,
                reference_genome=reference_genome,
            )
        )
        return ht.select("interval", "ppi_id")

    def _transform(ht: "hl.Table") -> "hl.Table":
        grouped = ht.group_by(ht.interval).aggregate(
            ppi_ids=hl.agg.collect_as_set(ht.ppi_id)
        )
        return grouped.annotate(ppi_ids=hl.sorted(hl.array(grouped.ppi_ids)))

    return _create_table_base(
        source_name="interactome",
        input_path=input_path,
        output_path=output_path,
        import_func=_import,
        transform_func=_transform,
        cleanup_func=lambda: _cleanup_temp_file(tsv_path),
        overwrite=overwrite,
        export_tsv=export_tsv,
    )


def create_clinvar_tb(
    input_path: str,
    output_path: str,
    overwrite: bool = False,
    export_tsv: bool = False,
    reference_genome: str = "GRCh38",
) -> "hl.Table":
    """
    Create a Hail Table from a ClinVar VCF file keyed by (locus, alleles).

    Example usage:
        ht = create_clinvar_tb(
            input_path="/path/to/clinvar.vcf.gz",
            output_path="/path/to/output.ht"
        )

    Parameters
    ----------
    input_path : str
        Path to the ClinVar VCF input file.
    output_path : str
        Path to write the output Hail Table.
    overwrite : bool, optional
        Whether to overwrite the output file if it exists (default: False).
    export_tsv : bool, optional
        If True, also export a flattened TSV version (default: False).
    reference_genome : str, optional
        Reference genome to use for parsing variants (default: "GRCh38").

    Returns
    -------
    hl.Table
        Hail Table keyed by (locus, alleles) with ClinVar annotations.

    Notes
    -----
    Uses force=True for import which may impact performance (single-threaded processing).
    """
    recode = contig_recoding()

    clinvar_tb = _create_table_base(
        source_name="ClinVar",
        input_path=input_path,
        output_path=output_path,
        import_func=lambda: hl.import_vcf(
            path=input_path,
            force=True,
            reference_genome=reference_genome,
            contig_recoding=recode,
            skip_invalid_loci=True,
        ).rows(),
        transform_func=lambda ht: ht.repartition(100).key_by("locus", "alleles"),
        overwrite=overwrite,
        export_tsv=False,  # Handle custom export below
    )

    if export_tsv:
        logger.info(f"Exporting flattened table to {output_path}.tsv.bgz")
        clinvar_tb.flatten().export(f"{output_path}.tsv.bgz")

    return clinvar_tb


def create_gevir_tb(
    input_path: str,
    output_path: str,
    fields: list = None,
    overwrite: bool = False,
    export_tsv: bool = False,
) -> "hl.Table":
    """
    Create a Hail Table from GEVIR gene metrics TSV file keyed by gene_id.

    Example usage:
        ht = create_gevir_tb(
            input_path="/path/to/gevir_metrics.txt",
            output_path="/path/to/output.ht",
            fields=["gene_id", "gevir_score", "gevir_rank"]
        )

    Parameters
    ----------
    input_path : str
        Path to the GEVIR gene metrics TSV input file.
    output_path : str
        Path to write the output Hail Table.
    fields : list, optional
        List of fields to select from the imported table (default: None, selects all fields).
    overwrite : bool, optional
        Whether to overwrite the output file if it exists (default: False).
    export_tsv : bool, optional
        If True, also export a TSV version (default: False).

    Returns
    -------
    hl.Table
        Hail Table keyed by gene_id with selected GEVIR metrics.

    Notes
    -----
    If schema is stable, consider specifying types=... instead of impute=True.
    """
    return _create_table_base(
        source_name="GEVIR",
        input_path=input_path,
        output_path=output_path,
        import_func=lambda: hl.import_table(
            paths=input_path, impute=True, min_partitions=100, key="gene_id"
        ),
        fields=fields,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )


def create_ensembl_gene_tb(
    input_path: str,
    output_path: str,
    fields: list = None,
    canonical: bool = True,
    overwrite: bool = False,
    export_tsv: bool = False,
) -> "hl.Table":
    """
    Create a Hail Table from an Ensembl BioMart gene annotation TSV file keyed by gene_id.

    Example usage:
        ht = create_ensembl_gene_tb(
            input_path="/path/to/ensembl_biomart_genes.txt",
            output_path="/path/to/output.ht",
            fields=["gene_id", "gene_name", "gene_type", "chromosome",
                    "gene_start", "gene_end", "transcript_id", "protein_id"]
        )

    Parameters
    ----------
    input_path : str
        Path to the Ensembl BioMart gene annotation TSV input file.
    output_path : str
        Path to write the output Hail Table.
    fields : list, optional
        List of fields to select from the imported table (default: None, selects all fields).
    canonical : bool, optional
        If True, filter to only canonical transcripts (default: True).
    overwrite : bool, optional
        Whether to overwrite the output file if it exists (default: False).
    export_tsv : bool, optional
        If True, also export a TSV version (default: False).

    Returns
    -------
    hl.Table
        Hail Table keyed by gene_id with selected Ensembl gene annotations.
    """

    def transform(ht: hl.Table) -> hl.Table:
        logger.info("Replacing field names")
        ht = ht.rename(ENSEMBL_BIOMART_FIELDS)

        if canonical:
            logger.info("Filtering canonical transcripts")
            ht = ht.filter(ht.canonical == "1", keep=True)

        logger.info("Grouping by gene_id")
        ht = (
            ht.group_by(ht.gene_id)
            .aggregate(
                transcript_id=hl.agg.collect_as_set(ht.transcript_id).filter(
                    lambda x: x != ""
                ),
                protein_id=hl.agg.collect_as_set(ht.protein_id).filter(
                    lambda x: x != ""
                ),
                gene_synonym=hl.agg.collect_as_set(ht.gene_synonym).filter(
                    lambda x: x != ""
                ),
                gene_name=hl.agg.take(hl.or_else(ht.gene_name, ""), 1)[0],
                chromosome=hl.agg.take(hl.or_else(ht.chromosome, ""), 1)[0],
                gene_start=hl.agg.take(
                    hl.or_else(ht.gene_start, hl.null(ht.gene_start.dtype)), 1
                )[0],
                gene_end=hl.agg.take(
                    hl.or_else(ht.gene_end, hl.null(ht.gene_end.dtype)), 1
                )[0],
                gene_type=hl.agg.take(hl.or_else(ht.gene_type, ""), 1)[0],
            )
            .key_by("gene_id")
        )
        return ht

    return _create_table_base(
        source_name="Ensembl gene",
        input_path=input_path,
        output_path=output_path,
        import_func=lambda: hl.import_table(
            paths=input_path, min_partitions=50, impute=True
        ),
        transform_func=transform,
        fields=fields,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )


def create_dbnsfp_tb(
    input_path: str,
    output_path: str,
    reference_genome: str = "GRCh38",
    overwrite: bool = False,
    export_tsv: bool = False,
    min_partitions: int = 200,
    force_bgz: bool = True,
    parse_transcript_scores: bool = True,
    group_prefixes: Optional[List[str]] = None,
    auto_convert_bgz: bool = False,
) -> "hl.Table":
    """
    Create a Hail Table from a dbNSFP variant TSV/BGZ file keyed by (locus, alleles).

    Example usage:
        ht = create_dbnsfp_tb(
            input_path="/path/to/dbNSFP.tsv.bgz",
            output_path="/path/to/output.ht"
        )

    Parameters
    ----------
    input_path : str
        Path to the dbNSFP TSV/BGZ input file.
    output_path : str
        Path to write the output Hail Table.
    reference_genome : str, optional
        Reference genome to use for parsing variants (default: "GRCh38").
    overwrite : bool, optional
        Whether to overwrite the output file if it exists (default: False).
    export_tsv : bool, optional
        If True, also export a flattened TSV version (default: False).
    min_partitions : int, optional
        Minimum number of partitions for import (default: 200).
    force_bgz : bool, optional
        If True, force bgzip compression for input (default: True).
    parse_transcript_scores : bool, optional
        If True, parse transcript-specific scores into dicts (default: True).
    group_prefixes : list of str, optional
        List of field prefixes to group into structs (default: common population/annotation prefixes).
    auto_convert_bgz : bool, optional
        If True, automatically convert plain gzip files to BGZF before import (default: False).

    Returns
    -------
    hl.Table
        Hail Table keyed by (locus, alleles) with parsed and grouped annotations.

    Notes
    -----
    Steps performed:
    - Import table with missing '.' and no type imputation
    - Build a variant key from '#chr', 'pos(1-based)', 'ref', 'alt' and parse to (locus, alleles)
    - Key the table by (locus, alleles)
    - Optionally map transcript-specific scores ending with '_score' or 'CADD_phred' to dict
    - Optionally group common prefixes (e.g., gnomAD, ExAC) into structs
    """

    # Resolve compression: detect gz vs bgzf, optionally convert
    input_path, force_bgz = resolve_compression(
        input_path,
        force_bgz=force_bgz,
        auto_convert=auto_convert_bgz,
    )

    def transform(ht: hl.Table) -> hl.Table:
        # Normalize chromosome field and construct variant key
        row_fields = get_row_fields(ht)
        if "#chr" in row_fields:
            ht = ht.rename({"#chr": "chr"})
        else:
            if "chr" not in row_fields:
                raise ValueError("dbNSFP input missing '#chr' or 'chr' column")

        _chr_str = hl.str(ht["chr"])
        ht = ht.annotate(
            chr=hl.if_else(
                _chr_str.lower().startswith("chr"), _chr_str, hl.str("chr") + _chr_str
            )
        )

        # Build variant_key: chr:pos:ref:alt
        row_fields = get_row_fields(ht)
        if (
            "pos(1-based)" not in row_fields
            or "ref" not in row_fields
            or "alt" not in row_fields
        ):
            raise ValueError(
                "dbNSFP input missing required columns: 'pos(1-based)', 'ref', or 'alt'"
            )

        variant_key_expr = hl.array(
            [ht.chr, hl.str(ht["pos(1-based)"]), ht.ref, ht.alt]
        )
        ht = ht.annotate(variant_key=hl.delimit(variant_key_expr, ":"))

        # Parse to locus/alleles
        ht = ht.annotate(
            **hl.parse_variant(ht.variant_key, reference_genome=reference_genome)
        )

        # Key the table and cleanup staging columns
        ht = ht.key_by("locus", "alleles")
        ht = ht.drop("variant_key", "chr", "pos(1-based)", "ref", "alt")

        # Transcript-specific score parsing
        row_fields = get_row_fields(ht)
        if parse_transcript_scores and "Ensembl_transcriptid" in row_fields:
            logger.info(
                "Parsing transcript-specific scores into dicts keyed by Ensembl_transcriptid"
            )
            ht = ht.annotate(Ensembl_transcriptid=hl.str(ht.Ensembl_transcriptid))
            ht = ht.annotate(Ensembl_transcriptid=ht.Ensembl_transcriptid.split(";"))

            row_fields_list = list(get_row_fields(ht))
            score_fields = [
                f for f in row_fields_list if f.endswith("_score") or f == "CADD_phred"
            ]

            def _to_float_array(s):
                s_def = hl.or_else(s, "")
                arr = s_def.split(";")
                return hl.map(hl.parse_float, arr)

            def _single_to_dict(val):
                return hl.dict(
                    hl.zip(
                        ht.Ensembl_transcriptid,
                        hl.map(lambda _x: hl.parse_float(val), ht.Ensembl_transcriptid),
                    )
                )

            ann = {}
            for f in score_fields:
                is_multi = hl.is_defined(ht[f]) & ht[f].contains(";")
                ann[f] = hl.if_else(
                    is_multi,
                    hl.dict(hl.zip(ht.Ensembl_transcriptid, _to_float_array(ht[f]))),
                    _single_to_dict(ht[f]),
                )
            if ann:
                ht = ht.annotate(**ann)

        # Group common prefixes into structs
        prefixes = group_prefixes or [
            "gnomAD",
            "ExAC",
            "1000Gp3",
            "ESP6500",
            "clinvar",
        ]
        for prefix in prefixes:
            row_fields_list = list(get_row_fields(ht))
            pref_fields = [
                f for f in row_fields_list if f != prefix and f.startswith(prefix)
            ]
            if pref_fields:
                logger.info(
                    f"Grouping {prefix}* fields into struct '{prefix}' ({len(pref_fields)} fields)"
                )
                ht = ht.annotate(
                    **{prefix: hl.struct(**{f: ht[f] for f in pref_fields})}
                )
                ht = ht.drop(*pref_fields)

        return ht

    dbnsfp_tb = _create_table_base(
        source_name="dbNSFP",
        input_path=input_path,
        output_path=output_path,
        import_func=lambda: hl.import_table(
            paths=input_path,
            min_partitions=min_partitions,
            impute=False,
            missing=".",
            force_bgz=force_bgz,
        ),
        transform_func=transform,
        overwrite=overwrite,
        export_tsv=False,  # Handle custom export below
    )

    if export_tsv:
        logger.info(f"Exporting flattened dbNSFP table to {output_path}.tsv.bgz")
        dbnsfp_tb.flatten().export(f"{output_path}.tsv.bgz")

    return dbnsfp_tb


def create_clingen_gene_disease_tb(
    input_path: str,
    output_path: str,
    key_by: str = "gene_disease",
    min_classification: Optional[str] = None,
    fields: Optional[List[str]] = None,
    overwrite: bool = False,
    export_tsv: bool = False,
) -> "hl.Table":
    """
    Create a Hail Table from a ClinGen Gene-Disease Validity CSV file.

    ClinGen provides curated gene-disease associations with evidence-based
    classifications (Definitive, Strong, Moderate, Limited, etc.).

    Example usage:
        # Default: keyed by (hgnc_id, mondo_id)
        ht = create_clingen_gene_disease_tb(
            input_path="/path/to/clingen.csv",
            output_path="/path/to/output.ht"
        )

        # Gene-level aggregation (for joining with other gene tables)
        ht = create_clingen_gene_disease_tb(
            input_path="/path/to/clingen.csv",
            output_path="/path/to/output.ht",
            key_by="gene",
            min_classification="Moderate"
        )

    Parameters
    ----------
    input_path : str
        Path to the ClinGen Gene-Disease Validity CSV file.
    output_path : str
        Path to write the output Hail Table.
    key_by : str, optional
        Keying strategy:
        - "gene_disease" (default): Key by (hgnc_id, mondo_id) - preserves full granularity
        - "gene": Aggregate diseases per gene, key by hgnc_id
    min_classification : str, optional
        Filter to classifications at or above this level. Valid values:
        "Definitive", "Strong", "Moderate", "Limited", "Disputed", "Refuted".
        If None, includes all classifications (default: None).
    fields : list of str, optional
        List of fields to select from the table (default: None, keeps all).
    overwrite : bool, optional
        Whether to overwrite the output file if it exists (default: False).
    export_tsv : bool, optional
        If True, also export a TSV version (default: False).

    Returns
    -------
    hl.Table
        Hail Table with ClinGen gene-disease validity annotations.

    Notes
    -----
    The ClinGen CSV file has a 6-line metadata header that is skipped during import.

    Classification hierarchy (from strongest to weakest):
    1. Definitive
    2. Strong
    3. Moderate
    4. Limited
    5. Disputed
    6. Refuted
    7. No Known Disease Relationship

    When key_by="gene", diseases are aggregated per gene with fields:
    - disease_labels: set of disease labels
    - disease_mondo_pairs: set of (disease_label, mondo_id) pairs
    - mondo_ids: set of MONDO IDs
    - classifications: set of classification levels
    - max_classification_level: numeric level of highest classification
    - max_classification_label: label of highest classification
    - n_diseases: count of associated diseases
    """
    if key_by not in ("gene_disease", "gene"):
        raise ValueError(f"key_by must be 'gene_disease' or 'gene', got: {key_by}")

    if (
        min_classification is not None
        and min_classification not in CLINGEN_CLASSIFICATION_LEVELS
    ):
        raise ValueError(
            f"min_classification must be one of {CLINGEN_CLASSIFICATION_LEVELS}, "
            f"got: {min_classification}"
        )

    # Preprocess CSV to extract header and data rows
    # ClinGen files have varying metadata lines before the actual header
    # The header row starts with "GENE SYMBOL" and separator rows contain "++++++"
    # Use Hadoop API to support cloud URIs (gs://, s3://) and distributed Spark clusters
    logger.info("Preprocessing ClinGen CSV: finding header and filtering metadata")

    # Use Hail's temp file utility to create a Hadoop-accessible temp file
    tmp_path = hl.utils.new_temp_file(prefix="clingen_", extension=".csv")
    try:
        # Stream the file line-by-line to find header and skip metadata/separators
        with hl.hadoop_open(input_path, "r") as f:
            with hl.hadoop_open(tmp_path, "w") as out:
                found_header = False
                for line in f:
                    # Skip separator lines (contain "++++++")
                    if "++++++" in line:
                        continue
                    # Look for header row (starts with "GENE SYMBOL" in quotes or unquoted)
                    if not found_header:
                        if '"GENE SYMBOL"' in line or line.startswith("GENE SYMBOL"):
                            found_header = True
                            out.write(line)
                        # Skip metadata lines before header
                        continue
                    # Write all data lines after header
                    out.write(line)

                if not found_header:
                    raise RuntimeError(
                        f'ClinGen header "GENE SYMBOL" not found in {input_path}'
                    )

        logger.info(f"Preprocessed file written to {tmp_path}")

    except Exception as e:
        # Clean up temp file if preprocessing fails
        _cleanup_temp_file(tmp_path)
        raise RuntimeError(f"Failed to preprocess ClinGen CSV: {e}") from e

    try:

        def import_func():
            return hl.import_table(
                paths=tmp_path,
                delimiter=",",
                quote='"',
                impute=False,
                min_partitions=10,
            )

        def transform(ht: hl.Table) -> hl.Table:
            # Rename fields to standardized names
            logger.info("Renaming fields to standardized names")
            rename_map = {
                k: v
                for k, v in CLINGEN_GENE_DISEASE_FIELDS.items()
                if k in get_row_fields(ht)
            }
            ht = ht.rename(rename_map)

            # Clean HGNC ID (strip "HGNC:" prefix)
            row_fields = get_row_fields(ht)
            if "hgnc_id" in row_fields:
                ht = ht.annotate(
                    hgnc_id=hl.if_else(
                        ht.hgnc_id.startswith("HGNC:"),
                        ht.hgnc_id.replace("HGNC:", ""),
                        ht.hgnc_id,
                    )
                )

            # Clean MONDO ID (strip "MONDO:" prefix if present)
            if "mondo_id" in row_fields:
                ht = ht.annotate(
                    mondo_id=hl.if_else(
                        ht.mondo_id.startswith("MONDO:"),
                        ht.mondo_id.replace("MONDO:", ""),
                        ht.mondo_id,
                    )
                )

            # Add classification level as numeric for filtering/sorting
            classification_order = {
                level: i for i, level in enumerate(CLINGEN_CLASSIFICATION_LEVELS)
            }
            ht = ht.annotate(
                classification_level=hl.literal(classification_order).get(
                    ht.classification, hl.len(CLINGEN_CLASSIFICATION_LEVELS)
                )
            )

            # Apply min_classification filter if specified
            if min_classification is not None:
                min_level = classification_order[min_classification]
                logger.info(
                    f"Filtering to classifications >= {min_classification} (level {min_level})"
                )
                ht = ht.filter(ht.classification_level <= min_level)

            # Apply keying strategy
            if key_by == "gene_disease":
                logger.info("Keying by (hgnc_id, mondo_id)")
                ht = ht.key_by("hgnc_id", "mondo_id")
            else:  # key_by == "gene"
                logger.info("Aggregating by gene (hgnc_id)")
                ht = (
                    ht.group_by("hgnc_id", "gene_symbol")
                    .aggregate(
                        disease_labels=hl.agg.collect_as_set(ht.disease_label),
                        disease_mondo_pairs=hl.agg.collect_as_set(
                            hl.struct(
                                disease_label=ht.disease_label, mondo_id=ht.mondo_id
                            )
                        ),
                        mondo_ids=hl.agg.collect_as_set(ht.mondo_id),
                        classifications=hl.agg.collect_as_set(ht.classification),
                        modes_of_inheritance=hl.agg.collect_as_set(
                            ht.mode_of_inheritance
                        ),
                        max_classification_level=hl.agg.min(ht.classification_level),
                        n_diseases=hl.agg.count(),
                    )
                    .key_by("hgnc_id")
                )
                # Add max classification label with safe index clamping
                # max_classification_level can be len(CLINGEN_CLASSIFICATION_LEVELS) for
                # unknown classifications, so we clamp to valid index range
                classification_labels = hl.literal(CLINGEN_CLASSIFICATION_LEVELS)
                safe_index = hl.min(
                    ht.max_classification_level, hl.len(classification_labels) - 1
                )
                ht = ht.annotate(
                    max_classification_label=classification_labels[safe_index]
                )

            return ht

        clingen_tb = _create_table_base(
            source_name="ClinGen Gene-Disease Validity",
            input_path=input_path,
            output_path=output_path,
            import_func=import_func,
            transform_func=transform,
            fields=fields,
            overwrite=overwrite,
            export_tsv=export_tsv,
        )

        return clingen_tb

    finally:
        # Clean up temp file
        _cleanup_temp_file(tmp_path)


def create_gencc_submissions_tb(
    input_path: str,
    output_path: str,
    key_by: str = "gene_disease_submitter",
    min_classification: Optional[str] = None,
    fields: Optional[List[str]] = None,
    overwrite: bool = False,
    export_tsv: bool = False,
) -> "hl.Table":
    """
    Create a Hail Table from a GenCC submissions TSV file.

    GenCC aggregates gene-disease validity assertions from 12+ submitting
    organizations (ClinGen, PanelApp, G2P, Orphanet, etc.).

    Parameters
    ----------
    input_path : str
        Path to the GenCC submissions TSV file.
    output_path : str
        Path to write the output Hail Table.
    key_by : str, optional
        Keying strategy:
        - "gene_disease_submitter" (default): Key by (hgnc_id, mondo_id, submitter)
        - "gene_disease": Aggregate across submitters, key by (hgnc_id, mondo_id)
        - "gene": Aggregate all diseases per gene, key by hgnc_id
    min_classification : str, optional
        Filter to classifications at or above this level.
    fields : list of str, optional
        List of fields to select from the table (default: None, keeps all).
    overwrite : bool, optional
        Whether to overwrite the output file if it exists (default: False).
    export_tsv : bool, optional
        If True, also export a TSV version (default: False).

    Returns
    -------
    hl.Table
        Hail Table with GenCC gene-disease-submitter validity annotations.
    """
    valid_keys = ("gene_disease_submitter", "gene_disease", "gene")
    if key_by not in valid_keys:
        raise ValueError(f"key_by must be one of {valid_keys}, got: {key_by}")

    if (
        min_classification is not None
        and min_classification not in GENCC_CLASSIFICATION_LEVELS
    ):
        raise ValueError(
            f"min_classification must be one of {GENCC_CLASSIFICATION_LEVELS}, "
            f"got: {min_classification}"
        )

    def import_func():
        return hl.import_table(
            paths=input_path,
            delimiter="\t",
            impute=False,
            min_partitions=10,
        )

    def transform(ht: hl.Table) -> hl.Table:
        # Rename fields to standardized names
        logger.info("Renaming GenCC fields to standardized names")
        rename_map = {
            k: v for k, v in GENCC_SUBMISSION_FIELDS.items() if k in get_row_fields(ht)
        }
        ht = ht.rename(rename_map)

        # Clean HGNC ID (strip "HGNC:" prefix)
        row_fields = get_row_fields(ht)
        if "hgnc_id" in row_fields:
            ht = ht.annotate(
                hgnc_id=hl.if_else(
                    ht.hgnc_id.startswith("HGNC:"),
                    ht.hgnc_id.replace("HGNC:", ""),
                    ht.hgnc_id,
                )
            )

        # Clean MONDO ID (strip "MONDO:" prefix if present)
        if "mondo_id" in row_fields:
            ht = ht.annotate(
                mondo_id=hl.if_else(
                    ht.mondo_id.startswith("MONDO:"),
                    ht.mondo_id.replace("MONDO:", ""),
                    ht.mondo_id,
                )
            )

        # Add classification level as numeric for filtering/sorting
        classification_order = {
            level: i for i, level in enumerate(GENCC_CLASSIFICATION_LEVELS)
        }
        ht = ht.annotate(
            classification_level=hl.literal(classification_order).get(
                ht.classification, hl.len(GENCC_CLASSIFICATION_LEVELS)
            )
        )

        # Apply min_classification filter if specified
        if min_classification is not None:
            min_level = classification_order[min_classification]
            logger.info(
                f"Filtering to classifications >= {min_classification} "
                f"(level {min_level})"
            )
            ht = ht.filter(ht.classification_level <= min_level)

        # Apply keying strategy
        if key_by == "gene_disease_submitter":
            logger.info("Keying by (hgnc_id, mondo_id, submitter)")
            ht = ht.key_by("hgnc_id", "mondo_id", "submitter")
        elif key_by == "gene_disease":
            logger.info("Aggregating by gene-disease (hgnc_id, mondo_id)")
            ht = (
                ht.group_by("hgnc_id", "mondo_id", "gene_symbol", "disease_label")
                .aggregate(
                    submitters=hl.agg.collect_as_set(ht.submitter),
                    classifications=hl.agg.collect_as_set(ht.classification),
                    modes_of_inheritance=hl.agg.collect_as_set(ht.mode_of_inheritance),
                    max_classification_level=hl.agg.min(ht.classification_level),
                )
                .key_by("hgnc_id", "mondo_id")
            )
            ht = ht.annotate(n_submitters=hl.len(ht.submitters))
            classification_labels = hl.literal(GENCC_CLASSIFICATION_LEVELS)
            safe_index = hl.min(
                ht.max_classification_level, hl.len(classification_labels) - 1
            )
            ht = ht.annotate(
                max_classification_label=classification_labels[safe_index],
                classification=classification_labels[safe_index],
                classification_level=ht.max_classification_level,
            )
        else:  # key_by == "gene"
            logger.info("Aggregating by gene (hgnc_id)")
            ht = (
                ht.group_by("hgnc_id", "gene_symbol")
                .aggregate(
                    disease_labels=hl.agg.collect_as_set(ht.disease_label),
                    disease_mondo_pairs=hl.agg.collect_as_set(
                        hl.struct(disease_label=ht.disease_label, mondo_id=ht.mondo_id)
                    ),
                    mondo_ids=hl.agg.collect_as_set(ht.mondo_id),
                    classifications=hl.agg.collect_as_set(ht.classification),
                    modes_of_inheritance=hl.agg.collect_as_set(ht.mode_of_inheritance),
                    submitters=hl.agg.collect_as_set(ht.submitter),
                    max_classification_level=hl.agg.min(ht.classification_level),
                )
                .key_by("hgnc_id")
            )
            ht = ht.annotate(
                n_diseases=hl.len(ht.mondo_ids),
                n_submitters=hl.len(ht.submitters),
            )
            classification_labels = hl.literal(GENCC_CLASSIFICATION_LEVELS)
            safe_index = hl.min(
                ht.max_classification_level, hl.len(classification_labels) - 1
            )
            ht = ht.annotate(max_classification_label=classification_labels[safe_index])

        return ht

    gencc_tb = _create_table_base(
        source_name="GenCC Submissions",
        input_path=input_path,
        output_path=output_path,
        import_func=import_func,
        transform_func=transform,
        fields=fields,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )

    return gencc_tb


def create_cosmic_cgc_tb(
    input_path: str,
    output_path: str,
    hgnc_path: Optional[str] = None,
    min_classification: Optional[str] = None,
    mutation_context: str = "both",
    fields: Optional[List[str]] = None,
    overwrite: bool = False,
    export_tsv: bool = False,
) -> "hl.Table":
    """Create a Hail Table from COSMIC Cancer Gene Census (CGC) data.

    Parameters
    ----------
    input_path : str
        Path to the CGC TSV file (e.g. ``Cosmic_Genes_v98_GRCh38.tsv.gz``).
    output_path : str
        Path to write the output Hail Table.
    hgnc_path : str, optional
        Path to an HGNC Hail Table for gene symbol -> hgnc_id resolution.
        If not provided, the table is keyed by gene_symbol.
    min_classification : str, optional
        Filter to Tier at or above this level (``"Tier 1"`` or ``"Tier 2"``).
    mutation_context : str
        Filter genes by mutation context: ``"somatic"``, ``"germline"``, or
        ``"both"`` (default).
    fields : list of str, optional
        Fields to retain in the output table.
    overwrite : bool
        Overwrite existing output.
    export_tsv : bool
        Also export a TSV version.

    Returns
    -------
    hl.Table
        Gene-keyed Hail Table with COSMIC CGC annotations.
    """
    if mutation_context not in COSMIC_MUTATION_CONTEXTS:
        raise ValueError(
            f"mutation_context must be one of {COSMIC_MUTATION_CONTEXTS}, "
            f"got: {mutation_context}"
        )

    if (
        min_classification is not None
        and min_classification not in COSMIC_CGC_CLASSIFICATION_LEVELS
    ):
        raise ValueError(
            f"min_classification must be one of "
            f"{COSMIC_CGC_CLASSIFICATION_LEVELS}, got: {min_classification}"
        )

    resolved_path, force_bgz = resolve_compression(input_path)

    def import_func():
        kwargs = dict(
            paths=resolved_path,
            delimiter="\t",
            impute=False,
            min_partitions=10,
        )
        if force_bgz:
            kwargs["force_bgz"] = True
        return hl.import_table(**kwargs)

    def transform(ht: hl.Table) -> hl.Table:
        # Rename fields to standardized names using flexible matching
        logger.info("Renaming COSMIC CGC fields to standardized names")
        rename_map = build_rename_map(COSMIC_CGC_FIELDS, get_row_fields(ht))
        ht = ht.rename(rename_map)

        # Normalize Tier: raw "1"/"2" -> "Tier 1"/"Tier 2"
        logger.info("Normalizing tier classification values")
        ht = ht.annotate(
            classification=hl.if_else(
                ht.classification.matches(r"^\d+$"),
                hl.literal("Tier ") + ht.classification,
                ht.classification,
            )
        )

        # Add classification_level numeric field
        classification_order = {
            level: i for i, level in enumerate(COSMIC_CGC_CLASSIFICATION_LEVELS)
        }
        ht = ht.annotate(
            classification_level=hl.literal(classification_order).get(
                ht.classification,
                hl.len(COSMIC_CGC_CLASSIFICATION_LEVELS),
            )
        )

        # Normalize boolean fields using general-purpose str_to_bool
        for bool_field in ("somatic", "germline", "hallmark"):
            if bool_field in get_row_fields(ht):
                ht = ht.annotate(**{bool_field: str_to_bool(ht[bool_field])})

        # Parse comma-separated multi-value fields into arrays
        multi_value_fields = [
            "tumour_types_somatic",
            "tumour_types_germline",
            "role_in_cancer",
            "mutation_types",
        ]
        for mv_field in multi_value_fields:
            if mv_field in get_row_fields(ht):
                ht = ht.annotate(
                    **{
                        mv_field: hl.if_else(
                            hl.is_defined(ht[mv_field]) & (ht[mv_field] != ""),
                            ht[mv_field]
                            .split(",")
                            .map(lambda x: x.strip())
                            .filter(lambda x: x != ""),
                            hl.empty_array(hl.tstr),
                        )
                    }
                )

        # Apply mutation_context filter
        if mutation_context == "somatic":
            logger.info("Filtering to somatic genes")
            ht = ht.filter(ht.somatic)
        elif mutation_context == "germline":
            logger.info("Filtering to germline genes")
            ht = ht.filter(ht.germline)

        # Apply min_classification filter
        if min_classification is not None:
            min_level = classification_order[min_classification]
            logger.info(
                f"Filtering to classifications >= {min_classification} "
                f"(level {min_level})"
            )
            ht = ht.filter(ht.classification_level <= min_level)

        # Resolve gene_symbol -> hgnc_id if HGNC table is available
        if hgnc_path is not None:
            logger.info(f"Resolving gene symbols to HGNC IDs using {hgnc_path}")
            from hvantk.data.gene_mapper import GeneMapper

            hgnc_ht = hl.read_table(hgnc_path)
            mapper = GeneMapper(hgnc_ht)
            symbols = set(ht.aggregate(hl.agg.collect_as_set(ht.gene_symbol)))
            mapping = mapper.map_to_hgnc(list(symbols), source_type="gene_symbol")
            mapping_literal = hl.literal(mapping)
            ht = ht.annotate(hgnc_id=mapping_literal.get(ht.gene_symbol))
            # Strip "HGNC:" prefix to match convention used by other builders
            ht = ht.annotate(
                hgnc_id=hl.if_else(
                    hl.is_defined(ht.hgnc_id) & ht.hgnc_id.startswith("HGNC:"),
                    ht.hgnc_id.replace("HGNC:", ""),
                    ht.hgnc_id,
                )
            )
            n_mapped = len([v for v in mapping.values() if v])
            n_unmapped = len(symbols) - n_mapped
            logger.info(f"Mapped {n_mapped}/{len(symbols)} gene symbols to HGNC IDs")
            if n_unmapped > 0:
                logger.warning(
                    f"{n_unmapped} genes could not be mapped to HGNC IDs "
                    "and will be excluded from the keyed output"
                )
            ht = ht.filter(hl.is_defined(ht.hgnc_id) & (ht.hgnc_id != ""))
            ht = ht.key_by("hgnc_id")
        else:
            logger.warning(
                "No HGNC path provided; keying by gene_symbol. "
                "Provide --hgnc-path for HGNC ID resolution."
            )
            ht = ht.key_by("gene_symbol")

        return ht

    return _create_table_base(
        source_name="COSMIC CGC",
        input_path=input_path,
        output_path=output_path,
        import_func=import_func,
        transform_func=transform,
        fields=fields,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )


def create_hgnc_gene_tb(
    input_path: str,
    output_path: str,
    include_withdrawn: bool = False,
    fields: Optional[List[str]] = None,
    overwrite: bool = False,
    export_tsv: bool = False,
) -> "hl.Table":
    """
    Create a Hail Table from HGNC gene nomenclature data keyed by hgnc_id.

    HGNC (HUGO Gene Nomenclature Committee) provides the authoritative source
    for human gene symbols and cross-references to other databases.

    Example usage:
        ht = create_hgnc_gene_tb(
            input_path="/data/hgnc_complete_set.txt",
            output_path="/tables/hgnc.ht"
        )

    Parameters
    ----------
    input_path : str
        Path to the HGNC complete set TSV file (hgnc_complete_set.txt).
    output_path : str
        Path to write the output Hail Table.
    include_withdrawn : bool, optional
        If True, include withdrawn/non-approved genes (default: False).
    fields : list of str, optional
        List of fields to select from the table (default: None, keeps all).
    overwrite : bool, optional
        Whether to overwrite the output file if it exists (default: False).
    export_tsv : bool, optional
        If True, also export a TSV version (default: False).

    Returns
    -------
    hl.Table
        Hail Table keyed by hgnc_id with gene nomenclature and cross-references.

    Notes
    -----
    The table includes:
    - Core identifiers: hgnc_id, gene_symbol, gene_name, status
    - Symbol history: alias_symbols, prev_symbols (as arrays)
    - Cross-references: ensembl_gene_id, entrez_id, uniprot_ids, etc.
    - Gene classification: locus_group, locus_type, gene_group
    - Disease/clinical links: omim_id, orphanet_id, gencc, mane_select

    Pipe-separated fields (alias_symbol, prev_symbol, uniprot_ids, etc.) are
    automatically parsed into arrays.
    """

    def transform(ht: hl.Table) -> hl.Table:
        # Rename fields to standardized names
        logger.info("Renaming HGNC fields to standardized names")
        row_fields = get_row_fields(ht)
        rename_map = {k: v for k, v in HGNC_GENE_FIELDS.items() if k in row_fields}
        ht = ht.rename(rename_map)

        # Filter to approved genes unless include_withdrawn is True
        if not include_withdrawn:
            logger.info("Filtering to approved genes only")
            ht = ht.filter(ht.status == "Approved")

        # Parse pipe-separated fields into arrays
        logger.info("Parsing pipe-separated fields into arrays")
        row_fields = get_row_fields(ht)
        for field in HGNC_PIPE_SEPARATED_FIELDS:
            if field in row_fields:
                # Split on pipe, filter empty strings
                ht = ht.annotate(
                    **{
                        field: hl.if_else(
                            hl.is_defined(ht[field]) & (ht[field] != ""),
                            ht[field].split("\\|").filter(lambda x: x != ""),
                            hl.empty_array(hl.tstr),
                        )
                    }
                )

        # Key by hgnc_id
        logger.info("Keying table by hgnc_id")
        ht = ht.key_by("hgnc_id")

        return ht

    return _create_table_base(
        source_name="HGNC gene nomenclature",
        input_path=input_path,
        output_path=output_path,
        import_func=lambda: hl.import_table(
            paths=input_path,
            impute=False,
            min_partitions=10,
            missing="",
        ),
        transform_func=transform,
        fields=fields,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )


def create_ptm_sites_tb(
    input_path: str,
    output_path: str,
    reference_genome: str = "GRCh38",
    flanking_codons: int = 5,
    overwrite: bool = False,
    export_tsv: bool = False,
    fields: Optional[List[str]] = None,
) -> "hl.Table":
    """
    Create a Hail Table of PTM sites in genomic coordinates.

    Input is a TSV produced by the PTM coordinate mapper
    (see ``hvantk.ptm.pipeline.map_ptm_sites``) with columns: chrom,
    codon_start, codon_end, strand, uniprot_id, gene_symbol, residue_pos,
    amino_acid, ptm_type, ptm_category, source_db, evidence_type,
    n_observations, tissue_type.

    ``tissue_type`` carries sample provenance for sources that distinguish
    it (e.g. CPTAC ``"normal"``/``"tumor"``); curated or bulk-MS sources
    emit an empty string.

    The table is keyed by locus (codon start position), with a
    ``flanking_interval`` field for proximity-based annotation joins.

    Parameters
    ----------
    input_path : str
        Path to the mapped PTM sites TSV file.
    output_path : str
        Path to write the output Hail Table.
    reference_genome : str, optional
        Reference genome (default: "GRCh38").
    flanking_codons : int, optional
        Number of flanking codons for proximal window (default: 5).
    overwrite : bool, optional
        Whether to overwrite existing file (default: False).
    export_tsv : bool, optional
        If True, also export TSV version (default: False).
    fields : list of str, optional
        List of fields to select (default: None, keeps all).

    Returns
    -------
    hl.Table
        The checkpointed Hail Table keyed by locus.
    """
    if flanking_codons < 0:
        raise ValueError(f"flanking_codons must be >= 0, got {flanking_codons}")

    def transform(ht):
        # Remap contig names to match GRCh38 (e.g., MT -> M for chrM)
        contig_remap = hl.dict({"MT": "M"})
        ht = ht.annotate(
            _contig=hl.str("chr") + contig_remap.get(ht.chrom, ht.chrom),
        )

        # Filter to valid contigs in the reference genome
        valid_contigs = hl.set(hl.literal(hl.get_reference(reference_genome).contigs))
        ht = ht.filter(valid_contigs.contains(ht._contig))

        # Parse locus from contig + codon_start
        ht = ht.annotate(
            locus=hl.locus(
                ht._contig,
                hl.int32(ht.codon_start),
                reference_genome=reference_genome,
            ),
        )

        # Cast numeric fields
        ht = ht.annotate(
            codon_start=hl.int32(ht.codon_start),
            codon_end=hl.int32(ht.codon_end),
            residue_pos=hl.int32(ht.residue_pos),
            n_observations=hl.int32(ht.n_observations),
        )

        # Add flanking interval (codon ± flanking_codons * 3 bp)
        flank_bp = flanking_codons * 3
        ref = hl.get_reference(reference_genome)
        chrom_lengths = hl.dict(hl.literal({c: ref.lengths[c] for c in ref.contigs}))
        ht = ht.annotate(
            flanking_interval=hl.locus_interval(
                ht._contig,
                hl.max(1, ht.codon_start - flank_bp),
                hl.min(chrom_lengths.get(ht._contig), ht.codon_end + flank_bp),
                reference_genome=reference_genome,
                includes_end=True,
            ),
        )

        ht = ht.drop("_contig")
        ht = ht.key_by("locus")
        return ht

    return _create_table_base(
        source_name="PTM sites",
        input_path=input_path,
        output_path=output_path,
        import_func=lambda: hl.import_table(
            paths=input_path,
            impute=False,
            min_partitions=16,
            types={
                "codon_start": hl.tstr,
                "codon_end": hl.tstr,
                "residue_pos": hl.tstr,
                "n_observations": hl.tstr,
            },
        ),
        transform_func=transform,
        fields=fields,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )


# ---------------------------------------------------------------------------
# QTL Table Builder Helpers
# ---------------------------------------------------------------------------


def _parse_gtex_variant_id(
    ht, variant_id_field="variant_id", reference_genome="GRCh38"
):
    """Parse GTEx variant IDs into ``locus`` and ``alleles``.

    Format: ``chr1_1000050_C_T_b38`` — the build suffix is discarded.
    Used by both eQTL and pQTL builders (GTEx/Fang share the same ID format).

    Contig names are normalised to match the reference genome:
    GRCh38 contigs use ``chr`` prefix, GRCh37 contigs omit it.
    """
    parts = ht[variant_id_field].split("_")
    raw_contig = parts[0]
    # Normalise contig for the target reference genome
    bare = hl.if_else(raw_contig.startswith("chr"), raw_contig[3:], raw_contig)
    contig = hl.if_else(
        hl.literal(reference_genome).startswith("GRCh38"),
        "chr" + bare,
        bare,
    )
    return ht.annotate(
        locus=hl.locus(contig, hl.int32(parts[1]), reference_genome=reference_genome),
        alleles=hl.array([parts[2], parts[3]]),
    )


def _strip_ensembl_version(gene_id_expr):
    """Strip Ensembl version suffix (``ENSG00000000003.15`` → ``ENSG00000000003``)."""
    return gene_id_expr.split("\\.")[0]


def _scan_tissue_files(input_path, extensions):
    """Return ``[(file_path, tissue_name), ...]`` from *input_path*.

    Tissue name is inferred from the filename prefix before the first dot
    (e.g. ``Brain_Cortex.v8.signif_variant_gene_pairs.txt.gz`` → ``Brain_Cortex``).
    Accepts a single file or a directory.
    """
    from pathlib import Path

    p = Path(input_path)
    if p.is_file():
        return [(str(p), p.stem.split(".")[0])]
    if not p.is_dir():
        raise FileNotFoundError(f"Not a file or directory: {input_path}")
    matches = []
    for ext in extensions:
        matches.extend(sorted(p.glob(f"*{ext}")))
    if not matches:
        raise FileNotFoundError(f"No files matching {extensions} in {input_path}")
    return [(str(f), f.stem.split(".")[0]) for f in matches]


# ---------------------------------------------------------------------------
# eQTL Builder
# ---------------------------------------------------------------------------


def _import_eqtl_gtex_parquet(input_path, tissue, reference_genome):
    """Import GTEx v11 eQTL parquet files via ``spark.read.parquet``."""
    from pathlib import Path
    from pyspark.sql import SparkSession

    files = _scan_tissue_files(input_path, [".parquet"])
    spark = SparkSession.builder.getOrCreate()

    tables = []
    for fp, tname in files:
        if tissue and tname != tissue:
            continue
        logger.info("Importing eQTL parquet: %s (tissue: %s)", fp, tname)
        # Spark requires absolute paths for parquet (prototype lesson #4).
        sdf = spark.read.parquet(str(Path(fp).resolve()))
        ht_part = hl.Table.from_spark(sdf)
        row_fields = list(ht_part.row)
        # v11 parquet has 'af' (ALT allele frequency, in-sample), not 'maf'.
        # 'maf' is derived as min(af, 1-af). Both are exposed so downstream
        # consumers can choose: af preserves direction (slope is per-ALT);
        # maf is the symmetric population frequency used in coloc / filtering.
        af_expr = (
            hl.float64(ht_part.af)
            if "af" in row_fields
            else (
                # Forward-compat: if a future release renames af → maf, fall
                # back to maf as the AF proxy (acknowledging the direction
                # loss — same caveat as v8).
                hl.float64(ht_part.maf)
                if "maf" in row_fields
                else hl.missing(hl.tfloat64)
            )
        )
        ht_part = ht_part.select(
            gene_id_raw=ht_part.phenotype_id,
            variant_id=ht_part.variant_id,
            beta=hl.float64(ht_part.slope),
            se=hl.float64(ht_part.slope_se),
            p_value=hl.float64(ht_part.pval_nominal),
            af=af_expr,
            maf=hl.min(af_expr, 1.0 - af_expr),
            tissue=tname,
            gene_symbol=hl.missing(hl.tstr),
        )
        tables.append(ht_part)

    if not tables:
        raise FileNotFoundError(f"No eQTL parquet files matched (tissue={tissue})")
    return tables[0].union(*tables[1:]) if len(tables) > 1 else tables[0]


def _import_eqtl_gtex_tsv(input_path, tissue):
    """Import GTEx v8 eQTL TSV files via ``hl.import_table``."""
    files = _scan_tissue_files(input_path, [".txt.gz", ".tsv.gz"])

    tables = []
    for fp, tname in files:
        if tissue and tname != tissue:
            continue
        logger.info("Importing eQTL TSV: %s (tissue: %s)", fp, tname)
        ht_part = hl.import_table(
            fp,
            force=True,
            types={
                "slope": hl.tfloat64,
                "slope_se": hl.tfloat64,
                "pval_nominal": hl.tfloat64,
                "maf": hl.tfloat64,
            },
        )
        row_fields = list(ht_part.row)
        # v8 TSV has 'maf' (minor allele frequency) directly, NOT 'af'.
        # We expose both fields for schema consistency with v11, but af is
        # set equal to maf as an approximation — v8 doesn't carry directional
        # (REF/ALT) allele-frequency info. Downstream consumers that need
        # true ALT-direction info should use v11 inputs. The skill §4
        # documents this caveat.
        maf_expr = (
            ht_part.maf if "maf" in row_fields else hl.missing(hl.tfloat64)
        )
        ht_part = ht_part.select(
            gene_id_raw=ht_part.gene_id,
            variant_id=ht_part.variant_id,
            beta=ht_part.slope,
            se=ht_part.slope_se,
            p_value=ht_part.pval_nominal,
            af=maf_expr,
            maf=maf_expr,
            tissue=tname,
            gene_symbol=hl.missing(hl.tstr),
        )
        tables.append(ht_part)

    if not tables:
        raise FileNotFoundError(f"No eQTL TSV files matched (tissue={tissue})")
    return tables[0].union(*tables[1:]) if len(tables) > 1 else tables[0]


def _import_eqtl_eqtlgen(input_path, reference_genome):
    """Import eQTLGen cis-eQTL summary statistics."""
    logger.info("Importing eQTLGen: %s", input_path)
    ht = hl.import_table(
        input_path,
        force=True,
        types={"Pvalue": hl.tfloat64, "Zscore": hl.tfloat64},
    )
    # Construct a GTEx-format variant_id so the shared parser can handle it.
    if reference_genome.startswith("GRCh38"):
        contig = hl.format("chr%s", ht.SNPChr)
    else:
        contig = ht.SNPChr
    ht = ht.select(
        gene_id_raw=ht.Gene,
        variant_id=hl.delimit(
            [contig, ht.SNPPos, ht.OtherAllele, ht.AssessedAllele, "b37"],
            "_",
        ),
        beta=hl.missing(hl.tfloat64),  # eQTLGen provides Z-score, not beta
        se=hl.missing(hl.tfloat64),
        p_value=ht.Pvalue,
        af=hl.missing(hl.tfloat64),  # eQTLGen doesn't distribute af / maf
        maf=hl.missing(hl.tfloat64),
        tissue="Blood",
        gene_symbol=ht.GeneSymbol,
    )
    return ht


def create_eqtl_tb(
    input_path: str,
    output_path: str,
    reference_genome: str = "GRCh38",
    source: str = "gtex_v11",
    tissue: Optional[str] = None,
    p_threshold: float = 5e-8,
    overwrite: bool = False,
    export_tsv: bool = False,
    fields: Optional[List[str]] = None,
) -> "hl.Table":
    """Build an eQTL Hail Table keyed by ``(locus, alleles, gene_id)``.

    One variant can be an eQTL for multiple genes; the ``gene_id`` key
    prevents information loss and enables correct cascade joins.

    Input can be a single file or a directory of per-tissue files.  Tissue
    name is inferred from the filename prefix before the first dot (e.g.
    ``Brain_Cortex.v8.signif_variant_gene_pairs.txt.gz`` → ``Brain_Cortex``).

    Supported sources:

    * ``gtex_v11`` — Parquet signif_pairs (``spark.read.parquet`` →
      ``hl.Table.from_spark``).
    * ``gtex_v8``  — TSV ``signif_variant_gene_pairs.txt.gz``
      (``hl.import_table``).
    * ``eqtlgen`` — TSV cis-eQTLs (single file, different column names).

    Gene-ID version suffixes are stripped for cross-table compatibility
    (``ENSG00000000003.15`` → ``ENSG00000000003``).

    Parameters
    ----------
    input_path : str
        Single file or directory of per-tissue eQTL files.
    output_path : str
        Output Hail Table path.
    reference_genome : str
        ``GRCh38`` or ``GRCh37``.
    source : str
        Data-source identifier.
    tissue : str, optional
        Restrict import to files matching this tissue name.
    p_threshold : float
        P-value cutoff.  Set to ``0`` to keep all pairs (for coloc).
    overwrite, export_tsv, fields
        Standard builder parameters.
    """
    from hvantk.qtlcascade.constants import EQTL_SOURCES

    if source not in EQTL_SOURCES:
        raise ValueError(f"Unknown eQTL source: {source!r}. Supported: {EQTL_SOURCES}")

    def import_func():
        if source == "gtex_v11":
            return _import_eqtl_gtex_parquet(input_path, tissue, reference_genome)
        if source == "gtex_v8":
            return _import_eqtl_gtex_tsv(input_path, tissue)
        return _import_eqtl_eqtlgen(input_path, reference_genome)

    def transform(ht):
        ht = _parse_gtex_variant_id(ht, "variant_id", reference_genome)
        ht = ht.annotate(gene_id=_strip_ensembl_version(ht.gene_id_raw))
        ht = ht.drop("gene_id_raw", "variant_id")
        if p_threshold > 0:
            ht = ht.filter(ht.p_value <= p_threshold)
        ht = ht.annotate(source=source, is_cis=True)
        ht = ht.key_by("locus", "alleles", "gene_id")
        return ht

    return _create_table_base(
        source_name=f"eQTL ({source})",
        input_path=input_path,
        output_path=output_path,
        import_func=import_func,
        transform_func=transform,
        fields=fields,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )


# ---------------------------------------------------------------------------
# pQTL Builder
# ---------------------------------------------------------------------------


def _import_pqtl_gtex_fang(input_path, tissue):
    """Import Fang et al. (2025) pQTL allpairs (space-delimited gzip).

    Columns: ``gene_name  SNP  CHR  BP  A1  NMISS  BETA  STAT  P``.
    Rows where ``STAT = 0`` are removed (cannot derive SE).
    """
    files = _scan_tissue_files(input_path, [".txt.gz", ".tsv.gz"])

    tables = []
    for fp, tname in files:
        if tissue and tname != tissue:
            continue
        logger.info("Importing pQTL allpairs: %s (tissue: %s)", fp, tname)
        ht_part = hl.import_table(
            fp,
            delimiter=" ",
            force=True,
            types={
                "BETA": hl.tfloat64,
                "STAT": hl.tfloat64,
                "P": hl.tfloat64,
            },
        )
        # STAT = 0 → SE undefined
        ht_part = ht_part.filter(ht_part.STAT != 0.0)
        ht_part = ht_part.select(
            gene_symbol=ht_part.gene_name,
            variant_id=ht_part.SNP,
            beta=ht_part.BETA,
            stat=ht_part.STAT,  # kept for SE derivation in transform
            p_value=ht_part.P,
            tissue=tname,
        )
        tables.append(ht_part)

    if not tables:
        raise FileNotFoundError(f"No pQTL allpairs files matched (tissue={tissue})")
    return tables[0].union(*tables[1:]) if len(tables) > 1 else tables[0]


def create_pqtl_tb(
    input_path: str,
    output_path: str,
    reference_genome: str = "GRCh38",
    source: str = "gtex_fang",
    tissue: Optional[str] = None,
    hgnc_ht: Optional[str] = None,
    no_gene_map: bool = False,
    p_threshold: Optional[float] = None,
    overwrite: bool = False,
    export_tsv: bool = False,
    fields: Optional[List[str]] = None,
) -> "hl.Table":
    """Build a pQTL Hail Table keyed by ``(locus, alleles, gene_id)``.

    For ``source='gtex_fang'``: Fang et al. (2025) allpairs files
    (space-delimited gzip, TMT mass spectrometry, 5 tissues).
    SE is derived as ``|BETA / STAT|`` (Fang files lack an SE column).

    Gene symbols are mapped to Ensembl gene IDs via the HGNC table
    and :class:`~hvantk.data.gene_mapper.GeneMapper`.  This is
    **required** because the cascade join uses ``(locus, alleles,
    gene_id)`` with Ensembl IDs on the eQTL side; raw gene symbols
    would produce zero matches.

    Pass ``no_gene_map=True`` to opt out of mapping for non-cascade
    use cases (the table will be keyed by raw gene symbol).

    Parameters
    ----------
    input_path : str
        Single file or directory of per-tissue pQTL files.
    output_path : str
        Output Hail Table path.
    reference_genome : str
        ``GRCh38`` or ``GRCh37``.
    source : str
        Data-source identifier.  Only ``gtex_fang`` is currently supported.
    tissue : str, optional
        Restrict to this tissue.
    hgnc_ht : str, optional
        Path to HGNC Hail Table (built by ``create_hgnc_gene_tb``).
        Required unless ``no_gene_map=True``.
    no_gene_map : bool
        Skip Ensembl mapping and key by raw gene symbol.  The resulting
        table will **not** join with eQTL tables in cascade analysis.
    p_threshold : float, optional
        P-value cutoff.  ``None`` keeps all pairs (for coloc).
    overwrite, export_tsv, fields
        Standard builder parameters.

    Raises
    ------
    ValueError
        If ``hgnc_ht`` is not provided and ``no_gene_map`` is ``False``.
    """
    from hvantk.qtlcascade.constants import PQTL_SOURCES

    if source not in PQTL_SOURCES:
        raise ValueError(f"Unknown pQTL source: {source!r}. Supported: {PQTL_SOURCES}")
    if source != "gtex_fang":
        raise NotImplementedError(
            f"pQTL source {source!r} is not yet implemented. "
            "Only 'gtex_fang' (Fang et al. 2025) is currently supported."
        )

    if not hgnc_ht and not no_gene_map:
        raise ValueError(
            "Ensembl gene mapping is required for cascade-compatible pQTL "
            "tables. Provide --hgnc-ht <path> (HGNC Hail Table built by "
            "'hvantk mktable hgnc-gene'). If you intentionally want a "
            "symbol-keyed table for non-cascade use, pass --no-gene-map."
        )

    def import_func():
        return _import_pqtl_gtex_fang(input_path, tissue)

    def transform(ht):
        ht = _parse_gtex_variant_id(ht, "variant_id", reference_genome)

        # SE = |BETA / STAT| (prototype lesson #5)
        ht = ht.annotate(se=hl.abs(ht.beta / ht.stat))
        ht = ht.drop("stat", "variant_id")

        # Gene-symbol → Ensembl-ID mapping via GeneMapper
        # (prototype lesson #2: use Hail Table join, NOT hl.literal,
        # to avoid IR poisoning).
        if hgnc_ht:
            from hvantk.data.gene_mapper import GeneMapper

            logger.info(
                "Mapping gene symbols → Ensembl IDs via GeneMapper (%s)",
                hgnc_ht,
            )
            hgnc_table = hl.read_table(hgnc_ht)
            mapper = GeneMapper(hgnc_table)
            ht = mapper.annotate_table(
                ht,
                source_field="gene_symbol",
                source_type="gene_symbol",
                fields_to_add=["ensembl_gene_id"],
            )
            ht = ht.annotate(
                gene_id=hl.or_else(ht.hgnc_ensembl_gene_id, ht.gene_symbol),
            )
            ht = ht.drop("hgnc_ensembl_gene_id")
        else:
            # no_gene_map=True path
            logger.warning(
                "--no-gene-map: using gene symbols as gene_id. "
                "This table will NOT join with eQTL tables in cascade "
                "analysis."
            )
            ht = ht.annotate(gene_id=ht.gene_symbol)

        if p_threshold is not None and p_threshold > 0:
            ht = ht.filter(ht.p_value <= p_threshold)

        ht = ht.annotate(source=source, is_cis=True)
        ht = ht.key_by("locus", "alleles", "gene_id")
        return ht

    return _create_table_base(
        source_name=f"pQTL ({source})",
        input_path=input_path,
        output_path=output_path,
        import_func=import_func,
        transform_func=transform,
        fields=fields,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )


# ---------------------------------------------------------------------------
# AlphaGenome Builder
# ---------------------------------------------------------------------------


def create_alphagenome_tb(
    input_path: str,
    output_path: str,
    config_path: str,
    no_resume: bool = False,
    overwrite: bool = False,
) -> "hl.Table":
    """Run AlphaGenome variant effect predictions and return a checkpointed table.

    Runs the AlphaGenomeStreamer to call the API for each variant, then
    writes predictions.json and checkpoint files to output_path (a
    directory), then builds a minimal checkpointed Hail Table keyed by
    ``(locus, alleles)`` from the input variants for builder-protocol
    compatibility.

    Parameters
    ----------
    input_path : str
        Path to Hail Table (.ht) or TSV with chrom/pos/ref/alt columns.
    output_path : str
        Output directory for prediction JSON outputs.
    config_path : str
        Path to AlphaGenome YAML config file.
    no_resume : bool
        If True, discard existing checkpoints and restart.
    overwrite : bool
        If True, overwrite existing output directory contents.
    Returns
    -------
    hl.Table
        Checkpointed Hail Table keyed by ``(locus, alleles)``.
    """
    from hvantk.data.alphagenome_streamer import AlphaGenomeStreamer

    if overwrite and os.path.isdir(output_path):
        import shutil

        shutil.rmtree(output_path)

    streamer = AlphaGenomeStreamer(
        input_path=input_path,
        output_dir=output_path,
        config_path=config_path,
        no_resume=no_resume,
    )
    streamer.setup()
    try:
        for _batch in streamer.stream():
            pass  # checkpointing handled internally
    finally:
        streamer.teardown()

    logger.info("Creating AlphaGenome variants table from %s", input_path)

    if input_path.endswith(".ht"):
        ht = hl.read_table(input_path)
        if "locus" not in ht.row or "alleles" not in ht.row:
            raise ValueError(
                "Input Hail Table must contain 'locus' and 'alleles' fields for "
                "AlphaGenome table builder output."
            )
    else:
        resolved_path, force_bgz = resolve_compression(input_path)
        import_kwargs = {"impute": True}
        if force_bgz:
            import_kwargs["force_bgz"] = True
        ht = hl.import_table(
            resolved_path,
            **import_kwargs,
        )
        required = {"chrom", "pos", "ref", "alt"}
        missing = required.difference(set(ht.row))
        if missing:
            raise ValueError(
                f"TSV input missing required columns for AlphaGenome: "
                f"{', '.join(sorted(missing))}"
            )
        ht = ht.annotate(
            locus=hl.locus(ht.chrom, hl.int(ht.pos), reference_genome="GRCh38"),
            alleles=[ht.ref, ht.alt],
        ).key_by("locus", "alleles")

    ht = ht.annotate_globals(
        hvantk_metadata=build_table_metadata("AlphaGenome", input_path, ht)
    )
    if os.path.isdir(output_path):
        table_output_path = os.path.join(output_path, "alphagenome_variants.ht")
    else:
        table_output_path = output_path
    logger.info("Checkpointing AlphaGenome variants table to %s", table_output_path)
    return ht.checkpoint(output=table_output_path, overwrite=overwrite)


# Map from raw GWAS Catalog v1.0 column names (whitespace + slash + brackets)
# to snake_case fields. See hvantk/skills/gwas-catalog/SKILL.md §4 + §7.
_GWAS_CATALOG_RENAME = {
    "DATE ADDED TO CATALOG": "date_added_to_catalog",
    "PUBMEDID": "pubmedid",
    "FIRST AUTHOR": "first_author",
    "DATE": "date",
    "JOURNAL": "journal",
    "LINK": "link",
    "STUDY": "study",
    "DISEASE/TRAIT": "disease_trait",
    "INITIAL SAMPLE SIZE": "initial_sample_size",
    "REPLICATION SAMPLE SIZE": "replication_sample_size",
    "REGION": "region",
    "CHR_ID": "chr_id",
    "CHR_POS": "chr_pos",
    "REPORTED GENE(S)": "reported_genes",
    "MAPPED_GENE": "mapped_gene",
    "UPSTREAM_GENE_ID": "upstream_gene_id",
    "DOWNSTREAM_GENE_ID": "downstream_gene_id",
    "SNP_GENE_IDS": "snp_gene_ids",
    "UPSTREAM_GENE_DISTANCE": "upstream_gene_distance",
    "DOWNSTREAM_GENE_DISTANCE": "downstream_gene_distance",
    "STRONGEST SNP-RISK ALLELE": "strongest_snp_risk_allele",
    "SNPS": "snps",
    "MERGED": "merged",
    "SNP_ID_CURRENT": "snp_id_current",
    "CONTEXT": "context",
    "INTERGENIC": "intergenic",
    "RISK ALLELE FREQUENCY": "risk_allele_frequency",
    "P-VALUE": "p_value",
    "PVALUE_MLOG": "pvalue_mlog",
    "P-VALUE (TEXT)": "p_value_text",
    "OR or BETA": "or_or_beta",
    "95% CI (TEXT)": "ci_95_text",
    "PLATFORM [SNPS PASSING QC]": "platform",
    "CNV": "cnv",
}


def create_gwas_catalog_tb(
    input_path: str,
    output_path: str,
    overwrite: bool = False,
    export_tsv: bool = False,
    reference_genome: str = "GRCh38",
) -> "hl.Table":
    """Create a Hail Table from the EBI GWAS Catalog v1.0 full-associations TSV.

    Implements the contract in ``hvantk/skills/gwas-catalog/SKILL.md`` (v1.0,
    34 columns, no ``MAPPED_TRAIT_URI``). Rows are keyed by ``(locus, alleles)``
    with ``alleles = [<risk_allele>, "N"]`` (sentinel ALT, judgment call #1).

    Two filters drop rows the schema cannot express cleanly (judgment calls
    #2 + #3): ``STRONGEST SNP-RISK ALLELE`` ending in ``-?`` (no risk allele
    recorded) and ``CHR_ID`` containing ``;`` (multi-chromosome / haplotype
    associations). Both losses are documented in the skill.

    Parameters
    ----------
    input_path : str
        Path to the unzipped GWAS Catalog full-associations TSV.
    output_path : str
        Path to write the output Hail Table (``.ht`` directory).
    overwrite : bool, optional
        Overwrite the output if present (default: False).
    export_tsv : bool, optional
        If True, also export a flattened TSV alongside the HT (default: False).
    reference_genome : str, optional
        Reference genome for ``hl.parse_locus`` (default: "GRCh38").

    Returns
    -------
    hl.Table
        Hail Table keyed by ``(locus, alleles)`` with 34 snake_case fields
        (plus the synthesized ``locus``, ``alleles``, and ``risk_allele``).
    """

    def transform(ht: hl.Table) -> hl.Table:
        logger.info("Renaming GWAS Catalog columns to snake_case")
        ht = ht.rename(_GWAS_CATALOG_RENAME)

        # Judgment call #2 (skill §4): drop rows with no risk allele.
        ht = ht.filter(~ht.strongest_snp_risk_allele.endswith("-?"))
        # Judgment call #3 (skill §4): drop non-canonical contigs — covers
        # ';'-separated haplotype rows ("6;7"), interaction pairs ("1 x 10"),
        # and any other malformed shapes. Replaces the narrower contains(';')
        # check from the initial tier-3 implementation; an empty CHR_ID also
        # fails the regex so no separate empty-string guard is needed.
        ht = ht.filter(ht.chr_id.matches("^(chr)?(\\d+|X|Y|MT?)$"))

        # Type coercions (everything arrives as string from impute=False).
        ht = ht.annotate(
            chr_pos=hl.int32(ht.chr_pos),
            p_value=hl.float64(ht.p_value),
            pvalue_mlog=hl.if_else(
                ht.pvalue_mlog == "",
                hl.missing(hl.tfloat64),
                hl.float64(ht.pvalue_mlog),
            ),
            or_or_beta=hl.if_else(
                ht.or_or_beta == "", hl.missing(hl.tfloat64), hl.float64(ht.or_or_beta)
            ),
            risk_allele_frequency=hl.if_else(
                ht.risk_allele_frequency == "",
                hl.missing(hl.tfloat64),
                # Coerce non-numeric tokens (e.g. "NR") to NA.
                hl.parse_float64(ht.risk_allele_frequency),
            ),
            upstream_gene_distance=hl.if_else(
                ht.upstream_gene_distance == "",
                hl.missing(hl.tfloat64),
                hl.float64(ht.upstream_gene_distance),
            ),
            downstream_gene_distance=hl.if_else(
                ht.downstream_gene_distance == "",
                hl.missing(hl.tfloat64),
                hl.float64(ht.downstream_gene_distance),
            ),
            pubmedid=hl.if_else(
                ht.pubmedid == "", hl.missing(hl.tint32), hl.int32(ht.pubmedid)
            ),
            snp_id_current=hl.if_else(
                ht.snp_id_current == "",
                hl.missing(hl.tint32),
                # snp_id_current sometimes holds non-int tokens; coerce gracefully.
                hl.parse_int32(ht.snp_id_current),
            ),
            merged=hl.if_else(
                ht.merged == "", hl.missing(hl.tint32), hl.int32(ht.merged)
            ),
            intergenic=str_to_bool(ht.intergenic),
            cnv=str_to_bool(ht.cnv),
        )

        # Synthesize key: locus + (risk_allele, "N") sentinel ALT.
        risk_allele = ht.strongest_snp_risk_allele.split("-")[-1]
        # GRCh38 contigs are 'chrN'; the catalog stores bare 'N' — prepend
        # 'chr' for autosomes/sex chroms unless the file already uses it.
        # Mitochondria: GRCh38 uses 'chrM' (not 'chrMT'); the catalog stores
        # 'MT', so normalize before prefixing.
        chr_id_norm = hl.case() \
            .when((ht.chr_id == "MT") | (ht.chr_id == "chrMT"), "M") \
            .default(ht.chr_id)
        contig = hl.if_else(chr_id_norm.startswith("chr"), chr_id_norm, "chr" + chr_id_norm)
        ht = ht.annotate(
            risk_allele=risk_allele,
            locus=hl.parse_locus(
                contig + ":" + hl.str(ht.chr_pos), reference_genome=reference_genome
            ),
        )
        ht = ht.annotate(alleles=[ht.risk_allele, "N"])
        ht = ht.key_by("locus", "alleles")
        return ht

    return _create_table_base(
        source_name="GWAS Catalog",
        input_path=input_path,
        output_path=output_path,
        import_func=lambda: hl.import_table(
            paths=input_path,
            delimiter="\t",
            quote=None,
            missing="",
            impute=False,
            min_partitions=4,
        ),
        transform_func=transform,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )


def create_msigdb_tb(
    input_path: str,
    output_path: str,
    overwrite: bool = False,
    export_tsv: bool = False,
) -> "hl.Table":
    """Create a Hail Table from an MSigDB GMT gene-set file keyed by set_name.

    Implements the contract in ``hvantk/skills/msigdb/SKILL.md``. The GMT
    format is tab-separated with variable-width rows: each row is one gene
    set, column 1 is the set name, column 2 is a description (a gsea-msigdb
    URL in MSigDB-issued files), and columns 3..N are gene members.

    Imported via ``hl.import_lines`` (one row per line, single ``text``
    field) because ``hl.import_table`` rejects variable column counts.
    The transform splits the line on ``\\t`` and slices the gene members
    into an ``array<str>``.

    Parameters
    ----------
    input_path : str
        Path to the GMT file (e.g., ``c2.cp.v2026.1.Hs.symbols.gmt``).
    output_path : str
        Path to write the output Hail Table (``.ht`` directory).
    overwrite : bool, optional
        Overwrite the output if present (default: False).
    export_tsv : bool, optional
        If True, also export a TSV alongside the HT (default: False).

    Returns
    -------
    hl.Table
        Hail Table keyed by ``set_name`` with fields:
          - ``set_name: str``
          - ``source_url: str`` (the GMT description column, verbatim)
          - ``genes: array<str>`` (gene members; symbols for ``.Hs.symbols.gmt``)
    """

    def transform(ht: hl.Table) -> hl.Table:
        # hl.import_lines yields rows with `file: str` and `text: str`.
        # Split on tab; slice [2:] for the variable-width gene-member tail.
        parts = ht.text.split("\t")
        ht = ht.select(
            set_name=parts[0],
            source_url=parts[1],
            genes=parts[2:],
        )
        # Defensive: drop blank lines (parts would be a 1-element array).
        ht = ht.filter(ht.set_name != "")
        ht = ht.key_by("set_name")
        return ht

    return _create_table_base(
        source_name="MSigDB gene sets",
        input_path=input_path,
        output_path=output_path,
        import_func=lambda: hl.import_lines(paths=input_path, min_partitions=4),
        transform_func=transform,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )
