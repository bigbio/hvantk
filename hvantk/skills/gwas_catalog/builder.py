"""Hail Table builder for the EBI GWAS Catalog v1.0 full-associations TSV.

This module owns ``create_gwas_catalog_tb``, the canonical builder that turns
the GWAS Catalog v1.0 full-associations TSV into a Hail Table keyed by
``(locus, alleles)`` with ``alleles = [<risk_allele>, "N"]`` (sentinel ALT,
judgment call #1 in the skill). It was migrated out of
:mod:`hvantk.core.builders.table` so that everything gwas-catalog-specific
(builder, drift probe, tests, fixtures, SKILL) lives under the plugin folder
at :mod:`hvantk.skills.gwas_catalog`.

The shared helper ``_create_table_base`` intentionally stays in
``hvantk.core.builders.table`` because it is reused by builders across many
data sources.
"""

from __future__ import annotations

import logging

import hail as hl

from hvantk.core.builders.table import _create_table_base
from hvantk.core.utils.table_utils import str_to_bool

logger = logging.getLogger(__name__)


# Map from raw GWAS Catalog v1.0 column names (whitespace + slash + brackets)
# to snake_case fields. See hvantk/skills/gwas_catalog/SKILL.md sections 4 + 7.
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

    Implements the contract in ``hvantk/skills/gwas_catalog/SKILL.md`` (v1.0,
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

    def transform(ht: "hl.Table") -> "hl.Table":
        logger.info("Renaming GWAS Catalog columns to snake_case")
        ht = ht.rename(_GWAS_CATALOG_RENAME)

        # Judgment call #2 (skill section 4): drop rows with no risk allele.
        ht = ht.filter(~ht.strongest_snp_risk_allele.endswith("-?"))
        # Judgment call #3 (skill section 4): drop non-canonical contigs --
        # covers ';'-separated haplotype rows ("6;7"), interaction pairs
        # ("1 x 10"), and any other malformed shapes. Replaces the narrower
        # contains(';') check from the initial tier-3 implementation; an empty
        # CHR_ID also fails the regex so no separate empty-string guard is
        # needed.
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
        # GRCh38 contigs are 'chrN'; the catalog stores bare 'N' -- prepend
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
