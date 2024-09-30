# eam
# 08.04.22

import hail as hl

from pyvatk.utils.dataset import (get_ccr_ht,
                                  get_gevir_ht,
                                  get_gene_expression_ht,
                                  get_ppi_ht,
                                  get_gene_ann_ht,
                                  get_dbnsfp_scores_ht,
                                  get_gnomad_metrics_ht,
                                  get_gnomad_af_ht,
                                  get_deg_ht,
                                  get_clinvar_ht,
                                  get_hca_ht)


def annotate_clinvar_clnsig(t: hl.Table) -> hl.Table:
    clinvar_ht = get_clinvar_ht()
    # Benign labels from Clinvar
    benign_label_clinvar = ['Benign/Likely_benign',
                            'Likely_benign',
                            'Benign']
    # Pathogenic labels from Clinvar
    pathogenic_label_clinvar = ['Pathogenic/Likely_pathogenic',
                                'Likely_pathogenic',
                                'Pathogenic']
    
    # Pre-calculate conditions for better readability and performance
    is_pathogenic = t.clinvar_clnsig.any(lambda x: hl.set(pathogenic_label_clinvar).contains(x))
    is_benign = t.clinvar_clnsig.any(lambda x: hl.set(benign_label_clinvar).contains(x))
    
    t = t.annotate(clinvar_clnsig=clinvar_ht[t.key].info.CLNSIG)
    t = t.annotate(clinvar_clnsig=hl.case()
                   .when(is_pathogenic, 'P')
                   .when(is_benign, 'B')
                   .or_missing()
                   )

    return t


def annotate_ccr(t: hl.Table) -> hl.Table:
    ccr_ht = get_ccr_ht()
    t = t.annotate(ccr_pct=ccr_ht[t.locus].ccr_pct)
    return t


def annotate_gevir(t: hl.Table,
                   gene_id_col: str, ) -> hl.Table:
    gevir_ht = (get_gevir_ht()
                .select('gevir_pct', 'virlof_pct')
                )
    t = t.annotate(**gevir_ht[t[gene_id_col]])
    return t


def annotate_rnaseq_expression(t: hl.Table,
                               gene_id_col: str,
                               organ: str = 'Heart') -> hl.Table:
    gene_expression_ht = get_gene_expression_ht(organ=organ)
    t = t.annotate(**gene_expression_ht[t[gene_id_col]])
    return t


def annotate_ppi(t: hl.Table) -> hl.Table:
    ppi_ht = get_ppi_ht()
    t = t.annotate(ppi_site=hl.int(hl.is_defined(ppi_ht[t.locus])))
    return t


def annotate_ensembl_gene(t: hl.Table,
                          gene_symbol_col: str) -> hl.Table:
    """
    Annotate gene and transcript ensembl canonical IDs given a gene symbol name.
    Includes searching for synonymous gene name.

    :param t: Hail Table
    :param gene_symbol_col: Column name with gene symbols
    :return: Hail Table
    """

    # Import and prepare gene table for annotation
    gene_ht = get_gene_ann_ht()
    gene_ht = (gene_ht
               .transmute(gene_aliases=gene_ht.Gene_Synonym.add(gene_ht.Gene))
               .explode('gene_aliases', name='Gene')
               .key_by('Gene')
               .select('GeneID', 'TranscriptID')
               )

    # Annotate table
    t = t.annotate(**gene_ht[t[gene_symbol_col]])
    return t


def annotate_dbnsfp_scores(t: hl.Table,
                           transcript_id_col: str) -> hl.Table:
    """
    Annotate transcript-specific deleterious scores from dbNSFP database.

    :param t: Hail Table keyed by `locus` and `alleles`
    :param transcript_id_col: Ensembl transcript ID column
    :return: Hail Table
    """

    # Import and parse dbNSFP dataset with annotation scores
    ht_scores = get_dbnsfp_scores_ht()
    scores_fields = [f for f in ht_scores.row if f.endswith('_score') or f == 'CADD_phred']
    ht_scores = (ht_scores
                 .select(*scores_fields)
                 )

    # Annotate scores taking into account the affected transcript.
    t = t.annotate(**ht_scores[t.key])
    t = t.annotate(**{f: t[f].get(t[transcript_id_col]) for f in scores_fields})

    return t


def annotate_gnomad_constraint_metrics(t: hl.Table,
                                       transcript_id_col: str) -> hl.Table:
    """
    Annotate transcript-specific loss-of-function and missense constraint metrics from gnomad.

    :param t: Hail Table
    :param transcript_id_col: Ensembl transcript ID column
    :return: Hail Table
    """
    gnomad_metrics = get_gnomad_metrics_ht()
    t = t.annotate(**gnomad_metrics[t[transcript_id_col]])
    return t


def annotate_degs(t: hl.Table,
                  gene_symbol_col: str,
                  clusters: list = ['C0', 'C5', 'C7', 'C10', 'C14']) -> hl.Table:
    """
    Annotate (1-True, 0-False) whether the gene is differentially expressed in
    cardiac-specific cell clusters.

    :param t: Hail Table
    :param gene_symbol_col: Column name with gene symbols
    :param clusters: Cell cluster ids to query (e.g. C0, C1...)

    :return: Hail Table
    """
    degs = get_deg_ht()

    t = t.annotate(sc_cluster_id=degs[t[gene_symbol_col]].cluster_id)

    t = (t
         .transmute(**{f'sc_cluster_{c}': hl.if_else(hl.is_defined(t.sc_cluster_id) & t.sc_cluster_id.contains(c),
                                                     1, 0)
                       for c in clusters})
         )

    return t


def annotate_hca(t: hl.Table,
                 gene_id_col: str,
                 cell_categories: tuple = ('atrial_cardiomyocyte',
                                           'endothelial',
                                           'fibroblast',
                                           'neuronal',
                                           'smooth_muscle_cell',
                                           'ventricular_cardiomyocyte')) -> hl.Table:
    """
    Annotate gene expression levels (mean umi/cell) per cell categories from HCA dataset (UCSC)

    :param t: Hail Table
    :param gene_id_col: Column name with gene symbols
    :param cell_categories: Cell categories to annotate (e.g. ...) TODO: make cell categories a constant

    :return: Hail Table
    """
    hca_tb = (get_hca_ht()
              .select(*cell_categories)
              )

    t = (t
         .annotate(hca=hl.struct(**hca_tb[t[gene_id_col]]))
         )

    return t


def annotate_gnomad_af(t: hl.Table) -> hl.Table:
    """
    Annotate allele frequencies from gnomad v3.0 (whole-genome).
    Annotate missing (absent) AF values as zero.

    :param t: Hail Table keyed by `locus` and `alleles`
    :return: Hail Table
    """

    # import gnomad table with allele frequency annotation
    gnomad_af = get_gnomad_af_ht()

    # define allele frequency annotation expression
    ann_expr = gnomad_af[t.key].AF

    t = t.annotate(gnomad_af_genomes=hl.if_else(hl.is_defined(ann_expr),
                                                ann_expr,
                                                hl.float(0)))
    return t


def annotate_variant_id(t: hl.Table,
                        field_name: str = 'vid') -> hl.Table:
    """
    Expected input dataset with bi-allelic variant, and fields `locus` and `alleles`.
    Annotate variant ids as follows 'chr:position:ref:alt'.

    :param field_name: variant id field name
    :param t: Hail table
    :return: Hail Table
    """

    variant_id_ann_exp = {
        field_name: hl.delimit([hl.str(t.locus.contig),
                                hl.str(t.locus.position),
                                hl.str(t.alleles[0]),
                                hl.str(t.alleles[1])],
                               delimiter=":")
    }

    return t.annotate(**variant_id_ann_exp)
