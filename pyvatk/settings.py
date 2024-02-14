import os


# context settings for click
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

# The global variables RAW_DATA_PATH and ANNOTATION_DATA_PATH are used to store
# the paths to the raw data and annotation data, respectively.
RAW_DATA_PATH = None
ANNOTATION_DATA_PATH = None


def set_raw_data_path(raw_data_path: str):
    """
    Set the global variable RAW_DATA_PATH to the specified raw data path.

    Args:
        raw_data_path (str): The path to the raw data.

    Raises:
        ValueError: If the raw_data_path is not a valid directory.

    Returns:
        str: The updated RAW_DATA_PATH.
    """
    global RAW_DATA_PATH
    if os.path.isdir(raw_data_path):
        RAW_DATA_PATH = raw_data_path
        return RAW_DATA_PATH
    else:
        raise ValueError("Invalid raw_data_path: {}".format(raw_data_path))


def set_annotation_data_path(annotation_data_path: str):
    """
    Set the global variable ANNOTATION_DATA_PATH to the specified annotation data path.

    Args:
        annotation_data_path (str): The path to the annotation data.

    Raises:
        ValueError: If the annotation_data_path is not a valid directory.

    Returns:
        str: The updated ANNOTATION_DATA_PATH.
    """
    global ANNOTATION_DATA_PATH
    if os.path.isdir(annotation_data_path):
        ANNOTATION_DATA_PATH = annotation_data_path
        return ANNOTATION_DATA_PATH
    else:
        raise ValueError("Invalid annotation_data_path: {}".format(annotation_data_path))


# A dictionary of raw data paths
RAW_DATA_PATHS = {
   'interactome_path':      f'{RAW_DATA_PATH}/interactome/Interactome_INSIDER_hg38_stripped.bed',
   'clinvar_path':          f'{RAW_DATA_PATH}/clinvar/clinvar_20220403.vcf.gz',
   'rnaseq_path':           f'{RAW_DATA_PATH}/rnaseq-expression/E-MTAB-6814.Human.CPM.txt',
   'gene_ann_path':         f'{RAW_DATA_PATH}/ensembl/gene.ensembl.canonical.042022.tsv',
   'gnomad_metrics_path':   f'{RAW_DATA_PATH}/gnomad/gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz',
   'gevir_path':            f'{RAW_DATA_PATH}/gevir/gevir_metrics_pmid31873297.tsv.txt',
   'scell_heart_path':      f'{RAW_DATA_PATH}/rnaseq-expression/deg_scell_heart_pmid31835037.tsv',
   'scell_hca_path':        f'{RAW_DATA_PATH}/rnaseq-expression/hca_cells_ucsc_042022.tsv'
}


# A dictionary of annotation data paths
ANNOTATION_DATA_PATHS = {}



