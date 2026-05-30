"""Plugin-local constants for the ucsc_cellbrowser skill.

Moved from hvantk/core/constants.py per issue #120 (clean-core principle).
"""
from pathlib import Path

# Base directory of the hvantk package (used to compute UCSC_JSON_FILE_PATH).
_PACKAGE_DIR = Path(__file__).resolve().parents[3]

# UCSC Cell Browser base URL
UCSC_CELL_BROWSER_BASE_URL = "https://cells.ucsc.edu"

# File names used by UCSC Cell Browser downloads
EXPRESSION_MATRIX_FILE_NAME = "exprMatrix.tsv.gz"
METADATA_FILE_NAME = "meta.tsv"

# Path to the JSON file containing the UCSC cell datasets
UCSC_JSON_FILE_PATH = _PACKAGE_DIR / "resources" / "cells_ucsc_datasets.json"

# UCSC gene and cell ID columns
UCSC_CELL_ID_COLUMN = "cell_id"
UCSC_GENE_COLUMN = "gene"
