"""Plugin-local constants for the pqtl skill.

Moved from hvantk/core/qtl_constants.py per issue #122 (clean-core principle).
"""

# ---------------------------------------------------------------------------
# Supported data source identifiers
# ---------------------------------------------------------------------------

# Only gtex_fang is currently implemented; others are planned.
PQTL_SOURCES = ("gtex_fang",)

# ---------------------------------------------------------------------------
# Fang et al. (2025) tissue-matched pQTL tissues
# ---------------------------------------------------------------------------

FANG_TISSUES = ("Colon", "Heart", "Liver", "Lung", "Thyroid")

# Mapping from Fang pQTL tissue names to GTEx eQTL sub-tissue names.
# Some pQTL tissues span multiple eQTL sub-tissues.
FANG_TISSUE_EQTL_MAPPING = {
    "Colon": ("Colon_Transverse", "Colon_Sigmoid"),
    "Heart": ("Heart_Left_Ventricle", "Heart_Atrial_Appendage"),
    "Liver": ("Liver",),
    "Lung": ("Lung",),
    "Thyroid": ("Thyroid",),
}
