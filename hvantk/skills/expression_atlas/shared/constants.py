"""Plugin-local constants for the expression_atlas skill.

Moved from hvantk/core/constants.py per issue #120 (clean-core principle).
"""
from pathlib import Path

# Base directory of the hvantk package (used to compute legacy JSON path).
_PACKAGE_DIR = Path(__file__).resolve().parents[3]

# Expression Atlas base URL
EXPRESSION_ATLAS_BASE_URL = "https://www.ebi.ac.uk/gxa/experiments-content"

# Legacy path for compatibility (deprecated): resources/expression_atlas.json
EXPRESSION_ATLAS_JSON_FILE_PATH = _PACKAGE_DIR / "resources" / "expression_atlas.json"
