# Core constants for the project

from pathlib import Path
import logging

logger = logging.getLogger(__name__)

# Base directory of the package (hvantk)
BASE_DIR = Path(__file__).resolve().parent
logger.debug(f"Base directory: {BASE_DIR}")

## Ensembl/biomart constants
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

## UCSC Cell Browser base URL
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

## Expression Atlas base URL
EXPRESSION_ATLAS_BASE_URL = "https://www.ebi.ac.uk/gxa/experiments-content"
logger.debug(f"Expression Atlas base URL: {EXPRESSION_ATLAS_BASE_URL}")

# Path to the new unified registry system
REGISTRY_ROOT_PATH = BASE_DIR.parent / "resources" / "registry"

# Backward compatibility - path to transcriptomics datasets (replaces expression_atlas.json)
TRANSCRIPTOMICS_DATASETS_PATH = REGISTRY_ROOT_PATH / "transcriptomics" / "datasets.json"

# Legacy path for compatibility (deprecated)
EXPRESSION_ATLAS_JSON_FILE_PATH = (
    BASE_DIR.parent / "resources" / "expression_atlas.json"
)

logger.debug(f"Registry root path: {REGISTRY_ROOT_PATH}")
logger.debug(f"Transcriptomics datasets path: {TRANSCRIPTOMICS_DATASETS_PATH}")
logger.debug(f"Legacy Expression Atlas path: {EXPRESSION_ATLAS_JSON_FILE_PATH}")

# Explicit public API for this module
__all__ = [
    "BASE_DIR",
    "ENSEMBL_BIOMART_FIELDS",
    "UCSC_CELL_BROWSER_BASE_URL",
    "EXPRESSION_MATRIX_FILE_NAME",
    "METADATA_FILE_NAME",
    "UCSC_JSON_FILE_PATH",
    "UCSC_CELL_ID_COLUMN",
    "UCSC_GENE_COLUMN",
    "EXPRESSION_ATLAS_BASE_URL",
    "REGISTRY_ROOT_PATH",
    "TRANSCRIPTOMICS_DATASETS_PATH",
    "EXPRESSION_ATLAS_JSON_FILE_PATH",  # Backward compatibility
]
