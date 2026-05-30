"""Plugin-local constants for the cptac skill.

Moved from hvantk/core/ptm_constants.py per issue #122 (clean-core principle).
"""

# CPTAC Phospho (via cptac Python package)
CPTAC_CANCER_TYPES = [
    "brca", "ccrcc", "coad", "gbm",
    "hnscc", "lscc", "luad", "ov", "pdac", "ucec",
]
CPTAC_CANCER_CLASS_MAP = {
    "brca": "Brca",
    "ccrcc": "Ccrcc",
    "coad": "Coad",
    "gbm": "Gbm",
    "hnscc": "Hnscc",
    "lscc": "Lscc",
    "luad": "Luad",
    "ov": "Ov",
    "pdac": "Pdac",
    "ucec": "Ucec",
}
