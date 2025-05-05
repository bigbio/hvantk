# Constants for the project

# Ensembl biomart fields mapping
# ['Gene stable ID', 'Transcript stable ID', 'Protein stable ID', 'Chromosome/scaffold name',
# 'Gene start (bp)', 'Gene end (bp)', 'Ensembl Canonical', 'Gene name', 'Gene type', 'Gene Synonym']

from pathlib import Path
import logging

logger = logging.getLogger(__name__)

# Base directory of the project
BASE_DIR = Path(__file__).resolve().parent
logger.debug(f"Base directory: {BASE_DIR}")


ENSEMBL_BIOMART_FIELDS = {
    "Gene stable ID": "gene_id",
    "Transcript stable ID": "transcript_id",
    "Protein stable ID": "protein_id",
    "Chromosome/scaffold name": "chromosome",
    "Gene start (bp)": "gene_start",
    "Gene end (bp)": "gene_end",
    "Ensembl Canonical": "canonical",
    "Gene name": "gene_name",
    "Gene type": "gene_type",
    "Gene Synonym": "gene_synonym",
}
logger.debug(f"Ensembl biomart fields: {ENSEMBL_BIOMART_FIELDS}")

# UCSC Cell Browser base URL
UCSC_CELL_BROWSER_BASE_URL = "https://cells.ucsc.edu"
logger.debug(f"UCSC Cell Browser base URL: {UCSC_CELL_BROWSER_BASE_URL}")
EXPRESSION_MATRIX_FILE_NAME = "exprMatrix.tsv.gz"
logger.debug(f"Expression matrix file name: {EXPRESSION_MATRIX_FILE_NAME}")
METADATA_FILE_NAME = "meta.tsv"
logger.debug(f"Metadata file name: {METADATA_FILE_NAME}")

# Path to the JSON file containing the UCSC cell datasets
UCSC_JSON_FILE_PATH = BASE_DIR.parent / "resources" / "cells_ucsc_datasets.json"
logger.debug(f"UCSC JSON file path: {UCSC_JSON_FILE_PATH}")

# UCSC gene and cell ID columns
UCSC_CELL_ID_COLUMN = "cell_id"
UCSC_GENE_COLUMN = "gene"
