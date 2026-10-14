"""Plugin-local constants for the ucsc_cellbrowser skill.

Moved from hvantk/core/constants.py per issue #120 (clean-core principle).
"""

from pathlib import Path

# Plugin root directory. This skill owns its operational dataset manifest
# (the pinned UCSC Cell Browser collection list) under its own data/ folder.
_PLUGIN_DIR = Path(__file__).resolve().parents[1]

# UCSC Cell Browser base URL
UCSC_CELL_BROWSER_BASE_URL = "https://cells.ucsc.edu"

# File names used by UCSC Cell Browser downloads
EXPRESSION_MATRIX_FILE_NAME = "exprMatrix.tsv.gz"
METADATA_FILE_NAME = "meta.tsv"

# Path to the plugin-owned JSON file containing the pinned UCSC cell datasets
UCSC_JSON_FILE_PATH = _PLUGIN_DIR / "data" / "cells_ucsc_datasets.json"

# UCSC gene and cell ID columns
UCSC_CELL_ID_COLUMN = "cell_id"
UCSC_GENE_COLUMN = "gene"

# Groups with fewer cells are dropped by `summarize_ucsc_streaming`; this is its default,
# and the default of `hvantk expression summarize-ucsc --min-cells`. The adult-ctx and
# dev-ctx drift probes apply the same cut-off to their parent collection's class labels.
SUMMARIZE_MIN_CELLS_PER_GROUP = 10
