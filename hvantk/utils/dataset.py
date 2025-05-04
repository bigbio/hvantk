# eam
# 06.04.22

import hail as hl
import logging

logger = logging.getLogger(__name__)

source_dir = None


def get_chd_denovo_ht() -> hl.Table:
    """
    Return a list of de novo mutations called from CHD trios.
    Curated from two studies, Jin 2017 and Sifrim-Hitz 2016.

    :return: Hail Table
    """
    logger.info("Getting CHD de novo HT")
    return hl.read_table(f"{source_dir}/data/ht/DNM_Jin2017_Sifrim2016_GRCh38_lift.ht")


def get_clinvar_ht() -> hl.Table:
    """
    Return Clinvar dataset (hg38)

    :return: Hail Table
    """
    logger.info("Getting Clinvar HT")
    return hl.read_table(f"{source_dir}/data/ht/clinvar.GRCh38.ht")


def get_gene_expression_ht(
    organ: str = "Heart", tp_col: str = "mean_expr_time_point"
) -> hl.Table:
    """
    Extract organ-specific gene expression levels per time points.

    :param organ: Organ of interest (e.g. Heart, Brain, Liver, ...)
    :param tp_col: Field with averaged expression value per time points

    :return: Hail Table
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
    logger.info("Getting CHD gene set")

    path = f"{source_dir}/resources/geneset/CHD_genes_all.tsv"
    t = hl.import_table(path, no_header=True)
    chd_gene_set = t.aggregate(hl.agg.collect_as_set(t.f0))

    return chd_gene_set


def get_gene_ann_ht() -> hl.Table:
    logger.info("Getting gene annotation HT")
    return hl.read_table(f"{source_dir}/data/ht/gene.ann.ensembl.ht")


def get_ccr_ht() -> hl.Table:
    logger.info("Getting CCR HT")
    return hl.read_table(f"{source_dir}/data/ht/ccr.GRCh38.ht")


def get_gevir_ht() -> hl.Table:
    logger.info("Getting GEVIR HT")
    return hl.read_table(f"{source_dir}/data/ht/gevir.metrics.ht")


def get_ppi_ht() -> hl.Table:
    logger.info("Getting PPI HT")
    return hl.read_table(f"{source_dir}/data/ht/interactome.GRCh38.ht")


def get_dbnsfp_scores_ht() -> hl.Table:
    logger.info("Getting dbNSFP scores HT")
    return hl.read_table(f"{source_dir}/data/ht/dbNSFP4.1a_variant.ht")


def get_gnomad_metrics_ht() -> hl.Table:
    logger.info("Getting gnomAD metrics HT")
    return hl.read_table(f"{source_dir}/data/ht/gnomad.metrics.ht")


def get_gnomad_af_ht() -> hl.Table:
    logger.info("Getting gnomAD AF HT")
    return hl.read_table(f"{source_dir}/data/ht/gnomad_3.0_sites_AF.ht")


def get_deg_ht() -> hl.Table:
    logger.info("Getting DEG HT")
    return hl.read_table(f"{source_dir}/data/ht/scell.heart.degs.ht")


def get_hca_ht() -> hl.Table:
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
        return "DataException: An error occurred while processing the data."
