# eam
# 06.04.22

import hail as hl
import logging

logger = logging.getLogger(__name__)

source_dir = None


def get_chd_denovo_ht() -> hl.Table:
    """
    Returns a Hail Table of de novo mutations identified in CHD trios.

    The dataset is curated from Jin 2017 and Sifrim-Hitz 2016 studies and contains de novo mutation calls relevant to congenital heart disease research.
    """
    logger.info("Getting CHD de novo HT")
    return hl.read_table(f"{source_dir}/data/ht/DNM_Jin2017_Sifrim2016_GRCh38_lift.ht")


def get_clinvar_ht() -> hl.Table:
    """
    Returns the Clinvar dataset as a Hail Table (hg38 reference).

    The table contains variant annotations from the Clinvar database mapped to the GRCh38 genome build.
    """
    logger.info("Getting Clinvar HT")
    return hl.read_table(f"{source_dir}/data/ht/clinvar.GRCh38.ht")


def get_gene_expression_ht(
    organ: str = "Heart", tp_col: str = "mean_expr_time_point"
) -> hl.Table:
    """
    Retrieves a Hail Table of gene expression values for a specified organ across time points.

    Args:
        organ: Name of the organ to extract expression data for (default is "Heart").
        tp_col: Column containing averaged expression values per time point (default is "mean_expr_time_point").

    Returns:
        A Hail Table keyed by gene, with columns for each time point containing expression values for the specified organ.
    """
    logger.info(f"Getting gene expression HT for organ: {organ}")

    # Import Hail Table with annotated expression values
    t = hl.read_table(f"{source_dir}/data/ht/rnaseq.human.ht")

    # getting available time point for the specified Organ
    tps = t[tp_col].key_set().filter(lambda x: x.organ == organ).time_point.collect()[0]

    # annotate expression values per time point
    t = t.annotate(
        **{
            f"{organ}.{tp}": t[tp_col].get(hl.struct(organ=organ, time_point=tp))
            for tp in tps
        }
    )

    t = t.drop(t["mean_expr_time_point"], t["mean_expr_dev_stage"]).key_by("Gene")

    return t


def get_chd_gene_set() -> hl.expr.SetExpression:
    """
    Retrieves the set of congenital heart disease (CHD) genes.

    Reads a gene set file containing CHD-associated genes and returns them as a Hail set expression.
    """
    logger.info("Getting CHD gene set")

    path = f"{source_dir}/resources/geneset/CHD_genes_all.tsv"
    t = hl.import_table(path, no_header=True)
    chd_gene_set = t.aggregate(hl.agg.collect_as_set(t.f0))

    return chd_gene_set


def get_gene_ann_ht() -> hl.Table:
    """
    Returns a Hail Table containing gene annotation data.

    The table includes gene-level annotations sourced from Ensembl.
    """
    logger.info("Getting gene annotation HT")
    return hl.read_table(f"{source_dir}/data/ht/gene.ann.ensembl.ht")


def get_ccr_ht() -> hl.Table:
    """
    Returns a Hail Table containing constrained coding region (CCR) data.

    The table provides CCR annotations for genomic regions, useful for variant interpretation.
    """
    logger.info("Getting CCR HT")
    return hl.read_table(f"{source_dir}/data/ht/ccr.GRCh38.ht")


def get_gevir_ht() -> hl.Table:
    """
    Returns the GEVIR metrics dataset as a Hail Table.

    The table contains gene-level GEVIR scores used for variant interpretation.
    """
    logger.info("Getting GEVIR HT")
    return hl.read_table(f"{source_dir}/data/ht/gevir.metrics.ht")


def get_ppi_ht() -> hl.Table:
    """
    Returns a Hail Table containing protein-protein interaction (PPI) data.

    The table provides curated interactome information mapped to the GRCh38 reference genome.
    """
    logger.info("Getting PPI HT")
    return hl.read_table(f"{source_dir}/data/ht/interactome.GRCh38.ht")


def get_dbnsfp_scores_ht() -> hl.Table:
    """
    Returns a Hail Table containing dbNSFP variant scores.

    The table includes functional prediction scores and annotations for genetic variants from the dbNSFP database.
    """
    logger.info("Getting dbNSFP scores HT")
    return hl.read_table(f"{source_dir}/data/ht/dbNSFP4.1a_variant.ht")


def get_gnomad_metrics_ht() -> hl.Table:
    """
    Returns a Hail Table containing gnomAD variant metrics.

    The table includes various metrics from the Genome Aggregation Database (gnomAD) for use in downstream analyses.
    """
    logger.info("Getting gnomAD metrics HT")
    return hl.read_table(f"{source_dir}/data/ht/gnomad.metrics.ht")


def get_gnomad_af_ht() -> hl.Table:
    """
    Returns a Hail Table containing gnomAD v3.0 allele frequency data.
    """
    logger.info("Getting gnomAD AF HT")
    return hl.read_table(f"{source_dir}/data/ht/gnomad_3.0_sites_AF.ht")


def get_deg_ht() -> hl.Table:
    """
    Returns a Hail Table containing differentially expressed genes (DEGs) from heart single-cell data.
    """
    logger.info("Getting DEG HT")
    return hl.read_table(f"{source_dir}/data/ht/scell.heart.degs.ht")


def get_hca_ht() -> hl.Table:
    """
    Returns a Hail Table containing Human Cell Atlas (HCA) heart data.

    The table includes single-cell transcriptomic data from the HCA heart project.
    """
    logger.info("Getting HCA HT")
    return hl.read_table(f"{source_dir}/data/ht/hca.heart.ht")


# Define a class to handle data exceptions and errors
class DataException(Exception):
    """
    Exception class for handling data-related errors.
    """

    def __init__(self, message):
        super().__init__(message)

    def __str__(self):
        base = super().__str__()
        if base:
            return f"DataException: {base}"
        return "DataException: An error occurred while processing the data."

