import os
import logging

logger = logging.getLogger(__name__)

# context settings for click
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])
logger.debug(f"Context settings: {CONTEXT_SETTINGS}")

# The global variables RAW_DATA_PATH and ANNOTATION_DATA_PATH are used to store
# the paths to the raw data and annotation data, respectively.
RAW_DATA_PATH = None
ANNOTATION_DATA_PATH = None


def set_raw_data_path(raw_data_path: str):
    """
    Sets the global RAW_DATA_PATH variable to the specified directory path.

    Args:
        raw_data_path: Directory path to use for raw data.

    Returns:
        The updated RAW_DATA_PATH.

    Raises:
        ValueError: If the provided path is not a valid directory.
    """
    global RAW_DATA_PATH
    if os.path.isdir(raw_data_path):
        RAW_DATA_PATH = raw_data_path
        logger.info(f"Setting RAW_DATA_PATH to {raw_data_path}")
        return RAW_DATA_PATH
    else:
        logger.error(f"Invalid raw_data_path: {raw_data_path}")
        raise ValueError("Invalid raw_data_path: {}".format(raw_data_path))


def set_annotation_data_path(annotation_data_path: str):
    """
    Sets the global annotation data path if the provided directory exists.

    Args:
        annotation_data_path: Path to the annotation data directory.

    Returns:
        The updated annotation data path.

    Raises:
        ValueError: If the provided path is not a valid directory.
    """
    global ANNOTATION_DATA_PATH
    if os.path.isdir(annotation_data_path):
        ANNOTATION_DATA_PATH = annotation_data_path
        logger.info(f"Setting ANNOTATION_DATA_PATH to {annotation_data_path}")
        return ANNOTATION_DATA_PATH
    else:
        logger.error(f"Invalid annotation_data_path: {annotation_data_path}")
        raise ValueError(
            "Invalid annotation_data_path: {}".format(annotation_data_path)
        )


# A dictionary of raw data paths
RAW_DATA_PATHS = {
    "interactome_path": f"{RAW_DATA_PATH}/interactome/Interactome_INSIDER_hg38_stripped.bed",
    "clinvar_path": f"{RAW_DATA_PATH}/clinvar/clinvar_20220403.vcf.gz",
    "rnaseq_path": f"{RAW_DATA_PATH}/rnaseq-expression/E-MTAB-6814.Human.CPM.txt",
    "gene_ann_path": f"{RAW_DATA_PATH}/ensembl/gene.ensembl.canonical.042022.tsv",
    "gnomad_metrics_path": f"{RAW_DATA_PATH}/gnomad/gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz",
    "gevir_path": f"{RAW_DATA_PATH}/gevir/gevir_metrics_pmid31873297.tsv.txt",
    "scell_heart_path": f"{RAW_DATA_PATH}/rnaseq-expression/deg_scell_heart_pmid31835037.tsv",
    "scell_hca_path": f"{RAW_DATA_PATH}/rnaseq-expression/hca_cells_ucsc_042022.tsv",
}
logger.debug(f"Raw data paths: {RAW_DATA_PATHS}")


# A dictionary of annotation data paths
ANNOTATION_DATA_PATHS = {}
logger.debug(f"Annotation data paths: {ANNOTATION_DATA_PATHS}")
