"""Hail Table builder for the EBI GWAS Catalog v1.0 full-associations TSV.

Owns the Phase B ``build_gwas_catalog_associations`` builder. Turns the
GWAS Catalog v1.0 full-associations TSV into an ``AnnotationTable`` keyed by
``(locus, alleles)`` with ``alleles = [<risk_allele>, "N"]`` (sentinel ALT,
judgment call #1 in the skill).
"""

from __future__ import annotations

import logging

import hail as hl

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


def build_gwas_catalog_associations(
    parsed_input,
    ctx,
    *,
    reference_genome: str = "GRCh38",
):
    """Phase B builder — returns an AnnotationTable keyed by
    ``(locus, alleles)`` with ``alleles = [<risk_allele>, "N"]`` (sentinel ALT).

    Two filters drop rows the schema cannot express cleanly (skill judgment
    calls #2 + #3): ``STRONGEST SNP-RISK ALLELE`` ending in ``-?`` (no risk
    allele recorded) and non-canonical ``CHR_ID`` (multi-chromosome /
    haplotype associations).
    """
    from hvantk.core.models import AnnotationTable

    ht = hl.import_table(
        paths=str(parsed_input),
        delimiter="\t",
        quote=None,
        missing="",
        impute=False,
        min_partitions=4,
    )

    ht = ht.rename(_GWAS_CATALOG_RENAME)
    ht = ht.filter(~ht.strongest_snp_risk_allele.endswith("-?"))
    ht = ht.filter(ht.chr_id.matches("^(chr)?(\\d+|X|Y|MT?)$"))

    ht = ht.annotate(
        chr_pos=hl.int32(ht.chr_pos),
        p_value=hl.float64(ht.p_value),
        pvalue_mlog=hl.if_else(
            ht.pvalue_mlog == "",
            hl.missing(hl.tfloat64),
            hl.float64(ht.pvalue_mlog),
        ),
        or_or_beta=hl.if_else(
            ht.or_or_beta == "",
            hl.missing(hl.tfloat64),
            hl.float64(ht.or_or_beta),
        ),
        risk_allele_frequency=hl.if_else(
            ht.risk_allele_frequency == "",
            hl.missing(hl.tfloat64),
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
            hl.parse_int32(ht.snp_id_current),
        ),
        merged=hl.if_else(
            ht.merged == "", hl.missing(hl.tint32), hl.int32(ht.merged)
        ),
        intergenic=str_to_bool(ht.intergenic),
        cnv=str_to_bool(ht.cnv),
    )

    risk_allele = ht.strongest_snp_risk_allele.split("-")[-1]
    chr_id_norm = hl.case() \
        .when((ht.chr_id == "MT") | (ht.chr_id == "chrMT"), "M") \
        .default(ht.chr_id)
    contig = hl.if_else(
        chr_id_norm.startswith("chr"), chr_id_norm, "chr" + chr_id_norm
    )
    ht = ht.annotate(
        risk_allele=risk_allele,
        locus=hl.parse_locus(
            contig + ":" + hl.str(ht.chr_pos), reference_genome=reference_genome
        ),
    )
    ht = ht.annotate(alleles=[ht.risk_allele, "N"])
    ht = ht.key_by("locus", "alleles")

    return AnnotationTable.from_hail(
        ht, provenance=ctx.provenance(schema_id="gwas-catalog-associations-v1")
    )
