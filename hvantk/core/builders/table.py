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
from hvantk.core.utils.table_utils import get_row_fields
from hvantk.core.models.metadata import build_table_metadata

logger = logging.getLogger(__name__)
_FILE_URI_PREFIX = "file://"

from hvantk.core.constants import ENSEMBL_BIOMART_FIELDS
from hvantk.core.utils.file_utils import resolve_compression
from hvantk.core.utils.genome import contig_recoding  # correct module import


def create_table_base(
    source_name: str,
    input_path: str,
    output_path: str,
    import_func: Callable[[], hl.Table],
    transform_func: Optional[Callable[[hl.Table], hl.Table]] = None,
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


__all__ = [
    "create_table_base",
    "cleanup_temp_file",
    "create_gnomad_constraint_gene_metrics_tb",
    "create_gevir_tb",
    "create_ensembl_gene_tb",
    "create_dbnsfp_tb",
]


def cleanup_temp_file(tmp_path: Optional[str]) -> None:
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
    return create_table_base(
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
    return create_table_base(
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

    return create_table_base(
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

    dbnsfp_tb = create_table_base(
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
    from hvantk.core.qtl_constants import PQTL_SOURCES

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
            "'hvantk reprocess hgnc:lookup'). If you intentionally want a "
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
            from hvantk.core.utils.gene_mapper import GeneMapper

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

    return create_table_base(
        source_name=f"pQTL ({source})",
        input_path=input_path,
        output_path=output_path,
        import_func=import_func,
        transform_func=transform,
        fields=fields,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )


# AlphaGenome builder moved to hvantk/skills/alphagenome/builder.py — the
# external API streamer drives the build end-to-end there, so core/builders
# no longer holds the alphagenome-specific path.
