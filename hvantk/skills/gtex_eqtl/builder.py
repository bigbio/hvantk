"""Hail Table builder for GTEx (and compatible) cis-eQTL summary statistics.

Owns the Phase B ``build_eqtl_associations`` builder. Turns a directory of
per-tissue GTEx v11 (parquet), GTEx v8 (TSV), or eQTLGen (TSV) cis-eQTL
summary statistics into an ``AnnotationTable`` triple-keyed by
``(locus, alleles, gene_id)``.

Shared GTEx variant-ID parsing helpers live in
``hvantk.core.utils.qtl_helpers`` (also used by the pQTL builder). The
source-specific import helpers below stay local to this plugin.
"""

from __future__ import annotations

import logging

import hail as hl

from hvantk.core.utils.qtl_helpers import (
    parse_gtex_variant_id,
    scan_tissue_files,
    strip_ensembl_version,
)

logger = logging.getLogger(__name__)


def _import_gtex_parquet(input_path, tissue, reference_genome):
    """Import GTEx v11 eQTL parquet files via ``spark.read.parquet``."""
    from pathlib import Path
    from pyspark.sql import SparkSession

    files = scan_tissue_files(input_path, [".parquet"])
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


def _import_gtex_tsv(input_path, tissue):
    """Import GTEx v8 eQTL TSV files via ``hl.import_table``."""
    files = scan_tissue_files(input_path, [".txt.gz", ".tsv.gz"])

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
        # true ALT-direction info should use v11 inputs.
        maf_expr = ht_part.maf if "maf" in row_fields else hl.missing(hl.tfloat64)
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


def _import_eqtlgen(input_path, reference_genome):
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


def build_eqtl_associations(
    parsed_input,
    ctx,
    *,
    reference_genome: str = "GRCh38",
    source: str = "gtex_v11",
    tissue=None,
    p_threshold: float = 5e-8,
    fields=None,
):
    """Phase B builder — returns an AnnotationTable keyed by
    ``(locus, alleles, gene_id)``.

    Supported sources:

    * ``gtex_v11`` — Parquet signif_pairs.
    * ``gtex_v8`` — TSV ``signif_variant_gene_pairs.txt.gz``.
    * ``eqtlgen`` — TSV cis-eQTLs.

    Gene-ID version suffixes are stripped for cross-table compatibility
    (``ENSG00000000003.15`` → ``ENSG00000000003``).
    """
    from hvantk.core.models import AnnotationTable
    from hvantk.core.qtl_constants import EQTL_SOURCES

    if source not in EQTL_SOURCES:
        raise ValueError(f"Unknown eQTL source: {source!r}. Supported: {EQTL_SOURCES}")

    if source == "gtex_v11":
        ht = _import_gtex_parquet(str(parsed_input), tissue, reference_genome)
    elif source == "gtex_v8":
        ht = _import_gtex_tsv(str(parsed_input), tissue)
    else:
        ht = _import_eqtlgen(str(parsed_input), reference_genome)

    ht = parse_gtex_variant_id(ht, "variant_id", reference_genome)
    ht = ht.annotate(gene_id=strip_ensembl_version(ht.gene_id_raw))
    ht = ht.drop("gene_id_raw", "variant_id")
    if p_threshold > 0:
        ht = ht.filter(ht.p_value <= p_threshold)
    ht = ht.annotate(source=source, is_cis=True)
    ht = ht.key_by("locus", "alleles", "gene_id")

    if fields is not None:
        ht = ht.select(*fields)

    return AnnotationTable.from_hail(
        ht, provenance=ctx.provenance(schema_id="gtex-eqtl-eqtls-v1")
    )
