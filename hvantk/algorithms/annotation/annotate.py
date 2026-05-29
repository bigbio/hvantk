# eam
# 08.04.22

import hail as hl
import logging

logger = logging.getLogger(__name__)

from hvantk.core.io.legacy_artifacts import (
    load_legacy_gene_expression_table,
    load_legacy_table,
)


def annotate_clinvar_clnsig(t: hl.Table) -> hl.Table:
    """
    Annotates variants with simplified ClinVar clinical significance labels.

    Variants are annotated with a clinical significance label based on ClinVar data: "P" for pathogenic, "B" for benign, or missing if neither applies. The annotation is determined by matching ClinVar CLNSIG values against predefined sets of pathogenic and benign labels.
    """
    # TEMP duplication — tracked by https://github.com/bigbio/hvantk/issues/133
    #
    # These label sets are also defined in
    # hvantk/skills/clinvar/shared/constants.py. Importing them from
    # there would violate the algorithms-must-not-import-from-skills
    # dependency guard.
    #
    # The proper fix — parameterizing the schema and vocabulary so this
    # function accepts any conformant pathogenicity-labeled table, not
    # just ClinVar — is tracked by issue #133. Remove this duplication
    # when that parameterization lands.
    CLINVAR_PATHOGENIC_LABELS = [
        "Pathogenic/Likely_pathogenic",
        "Likely_pathogenic",
        "Pathogenic",
    ]
    CLINVAR_BENIGN_LABELS = [
        "Benign/Likely_benign",
        "Likely_benign",
        "Benign",
    ]

    logger.info("Annotating ClinVar CLNSIG")
    clinvar_ht = load_legacy_table("clinvar")

    # First annotate clinvar_clnsig from the ClinVar table
    t = t.annotate(clinvar_clnsig=clinvar_ht[t.key].info.CLNSIG)

    # Now compute conditions using the annotated field
    is_pathogenic = t.clinvar_clnsig.any(
        lambda x: hl.set(CLINVAR_PATHOGENIC_LABELS).contains(x)
    )
    is_benign = t.clinvar_clnsig.any(
        lambda x: hl.set(CLINVAR_BENIGN_LABELS).contains(x)
    )

    t = t.annotate(
        clinvar_clnsig=hl.case()
        .when(is_pathogenic, "P")
        .when(is_benign, "B")
        .or_missing()
    )

    return t


def annotate_ccr(t: hl.Table) -> hl.Table:
    """
    Annotates variants with constrained coding region (CCR) percentile scores.

    Adds a `ccr_pct` field to the input Hail Table by joining on the variant locus.
    """
    logger.info("Annotating CCR")
    ccr_ht = load_legacy_table("ccr")
    t = t.annotate(ccr_pct=ccr_ht[t.locus].ccr_pct)
    return t


def annotate_gevir(
    t: hl.Table,
    gene_id_col: str,
) -> hl.Table:
    """
    Annotates a Hail Table with GEVIR gene-level constraint metrics.

    Adds the `gevir_pct` and `virlof_pct` fields to each row by joining on the specified gene ID column.
    """
    logger.info("Annotating GEVIR")
    gevir_ht = load_legacy_table("gevir").select("gevir_pct", "virlof_pct")
    t = t.annotate(**gevir_ht[t[gene_id_col]])
    return t


def annotate_rnaseq_expression(
    t: hl.Table, gene_id_col: str, organ: str = "Heart"
) -> hl.Table:
    """
    Annotates a Hail Table with RNA-seq gene expression data for a specified organ.

    Args:
        gene_id_col: Name of the column containing gene IDs to match against the expression dataset.
        organ: Name of the organ for which gene expression data should be used (default is "Heart").

    Returns:
        The input table annotated with gene expression values from the selected organ.
    """
    logger.info("Annotating RNAseq expression")
    gene_expression_ht = load_legacy_gene_expression_table(organ=organ)
    t = t.annotate(**gene_expression_ht[t[gene_id_col]])
    return t


def annotate_ppi(t: hl.Table) -> hl.Table:
    """
    Annotates variants with a binary flag indicating presence at a protein-protein interaction (PPI) site.

    Adds a `ppi_site` field set to 1 if the variant locus is present in the PPI dataset, or 0 otherwise.

    Returns:
        A Hail Table with the `ppi_site` annotation.
    """
    logger.info("Annotating PPI")
    ppi_ht = load_legacy_table("ppi")
    t = t.annotate(ppi_site=hl.int(hl.is_defined(ppi_ht[t.locus])))
    return t


def annotate_ensembl_gene(t: hl.Table, gene_symbol_col: str) -> hl.Table:
    """
    Annotates a table with Ensembl gene and canonical transcript IDs based on gene symbols.

    Searches for gene symbols and their synonyms to map each entry to its corresponding Ensembl gene and transcript IDs.
    """
    logger.info("Annotating Ensembl gene")

    # Import and prepare gene table for annotation
    gene_ht = load_legacy_table("gene_ann")
    gene_ht = (
        gene_ht.transmute(gene_aliases=gene_ht.Gene_Synonym.add(gene_ht.Gene))
        .explode("gene_aliases", name="Gene")
        .key_by("Gene")
        .select("GeneID", "TranscriptID")
    )

    # Annotate table
    t = t.annotate(**gene_ht[t[gene_symbol_col]])
    return t


def annotate_dbnsfp_scores(t: hl.Table, transcript_id_col: str) -> hl.Table:
    """
    Annotates variants with transcript-specific deleteriousness scores from the dbNSFP database.

    The function joins transcript-level scores (e.g., CADD, REVEL) from dbNSFP to the input Hail Table keyed by variant, then extracts the score corresponding to the specified transcript ID column.

    Args:
        t: Hail Table keyed by 'locus' and 'alleles'.
        transcript_id_col: Name of the column containing Ensembl transcript IDs.

    Returns:
        Hail Table annotated with transcript-specific deleteriousness scores.
    """
    logger.info("Annotating dbNSFP scores")

    # Import and parse dbNSFP dataset with annotation scores
    ht_scores = load_legacy_table("dbnsfp_scores")
    scores_fields = [
        f for f in ht_scores.row if f.endswith("_score") or f == "CADD_phred"
    ]
    ht_scores = ht_scores.select(*scores_fields)

    # Annotate scores taking into account the affected transcript.
    t = t.annotate(**ht_scores[t.key])
    t = t.annotate(**{f: t[f].get(t[transcript_id_col]) for f in scores_fields})

    return t


def annotate_gnomad_constraint_metrics(t: hl.Table, transcript_id_col: str) -> hl.Table:
    """
    Annotates transcript-level loss-of-function and missense constraint metrics from gnomAD.

    Args:
        t: Input Hail Table to annotate.
        transcript_id_col: Name of the column containing Ensembl transcript IDs.

    Returns:
        Hail Table annotated with gnomAD constraint metrics for each transcript.
    """
    logger.info("Annotating gnomAD constraint metrics")
    gnomad_metrics = load_legacy_table("gnomad_metrics")
    t = t.annotate(**gnomad_metrics[t[transcript_id_col]])
    return t


def annotate_degs(
    t: hl.Table, gene_symbol_col: str, clusters: list = ["C0", "C5", "C7", "C10", "C14"]
) -> hl.Table:
    """
    Annotates genes with binary flags indicating differential expression in specified cardiac cell clusters.

    Args:
        t: Input Hail Table containing gene information.
        gene_symbol_col: Name of the column containing gene symbols.
        clusters: List of cardiac cell cluster IDs to annotate (default: ["C0", "C5", "C7", "C10", "C14"]).

    Returns:
        Hail Table with additional fields for each cluster, set to 1 if the gene is differentially expressed in that cluster, 0 otherwise.
    """
    logger.info("Annotating DEGs")
    degs = load_legacy_table("deg")

    t = t.annotate(sc_cluster_id=degs[t[gene_symbol_col]].cluster_id)

    t = t.transmute(
        **{
            f"sc_cluster_{c}": hl.if_else(
                hl.is_defined(t.sc_cluster_id) & t.sc_cluster_id.contains(c), 1, 0
            )
            for c in clusters
        }
    )

    return t


def annotate_hca(
    t: hl.Table,
    gene_id_col: str,
    cell_categories: tuple = (
        "atrial_cardiomyocyte",
        "endothelial",
        "fibroblast",
        "neuronal",
        "smooth_muscle_cell",
        "ventricular_cardiomyocyte",
    ),
) -> hl.Table:
    """
    Annotates gene expression levels across specified cell categories from the Human Cell Atlas.

    Args:
        t: Input Hail Table.
        gene_id_col: Name of the column containing gene IDs.
        cell_categories: Tuple of cell category names to annotate (default includes major cardiac and related cell types).

    Returns:
        Hail Table annotated with a struct of mean UMI per cell for each specified cell category under the 'hca' field.
    """
    logger.info("Annotating HCA")
    hca_tb = load_legacy_table("hca").select(*cell_categories)

    t = t.annotate(hca=hl.struct(**hca_tb[t[gene_id_col]]))

    return t


def annotate_gnomad_af(t: hl.Table) -> hl.Table:
    """
    Annotates variants with allele frequencies from gnomAD v3.0 whole-genome data.

    Missing allele frequencies are annotated as zero.
    """
    logger.info("Annotating gnomAD AF")

    # import gnomad table with allele frequency annotation
    gnomad_af = load_legacy_table("gnomad_af")

    # define allele frequency annotation expression
    ann_expr = gnomad_af[t.key].AF

    t = t.annotate(
        gnomad_af_genomes=hl.if_else(hl.is_defined(ann_expr), ann_expr, hl.float(0))
    )
    return t


def annotate_variant_id(t: hl.Table, field_name: str = "vid") -> hl.Table:
    """
    Annotates each variant with a string identifier in the format 'chr:position:ref:alt'.

    Args:
        t: Hail Table containing bi-allelic variants with 'locus' and 'alleles' fields.
        field_name: Name of the field to store the variant ID (default is "vid").

    Returns:
        Hail Table with an added field containing the variant ID string.
    """
    logger.info("Annotating variant ID")

    variant_id_ann_exp = {
        field_name: hl.delimit(
            [
                hl.str(t.locus.contig),
                hl.str(t.locus.position),
                hl.str(t.alleles[0]),
                hl.str(t.alleles[1]),
            ],
            delimiter=":",
        )
    }

    return t.annotate(**variant_id_ann_exp)
