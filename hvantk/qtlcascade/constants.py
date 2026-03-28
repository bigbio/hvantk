"""
Constants and default parameters for QTL cascade analysis.

References
----------
- Giambartolomei et al. (2014) PLoS Genet 10(5):e1004383 — coloc priors
- Wakefield (2009) Am J Hum Genet 84(1):60-71 — ABF W parameter
- Fang et al. (2025) — tissue-matched eQTL/pQTL data (eGTEx)
"""

# ---------------------------------------------------------------------------
# Cascade classification
# ---------------------------------------------------------------------------

# Ordered for consistent plotting; tuple to prevent accidental mutation
CASCADE_CLASSES = ("eqtl_mediated", "discordant", "eqtl_only", "pqtl_only")

CASCADE_CLASS_LABELS = {
    "eqtl_mediated": "eQTL-mediated pQTL",
    "discordant": "Discordant (eQTL \u2260 pQTL direction)",
    "eqtl_only": "eQTL only (no pQTL)",
    "pqtl_only": "pQTL only (no eQTL)",
}

CASCADE_CLASS_COLORS = {
    "eqtl_mediated": "#2196F3",
    "discordant": "#F44336",
    "eqtl_only": "#9E9E9E",
    "pqtl_only": "#FF9800",
}

# ---------------------------------------------------------------------------
# Significance thresholds
# ---------------------------------------------------------------------------

DEFAULT_EQTL_P_THRESHOLD = 5e-8
DEFAULT_PQTL_P_THRESHOLD = 5e-8

# ---------------------------------------------------------------------------
# Colocalization priors — Giambartolomei et al. (2014) Table 1
# ---------------------------------------------------------------------------

DEFAULT_COLOC_P1 = 1e-4  # P(variant causal for trait 1 only)
DEFAULT_COLOC_P2 = 1e-4  # P(variant causal for trait 2 only)
DEFAULT_COLOC_P12 = 1e-5  # P(variant causal for both traits)

# Prior variance on true effect size — Wakefield (2009)
# 0.04 is appropriate for quantitative-trait QTLs
DEFAULT_COLOC_W = 0.04

# Posterior threshold for declaring colocalization
DEFAULT_COLOC_H4_THRESHOLD = 0.8

# Regional window for coloc (±kb from lead variant)
DEFAULT_COLOC_WINDOW_KB = 500

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

# ---------------------------------------------------------------------------
# Supported data source identifiers
# ---------------------------------------------------------------------------

EQTL_SOURCES = ("gtex_v11", "gtex_v8", "eqtlgen")
# Only gtex_fang is currently implemented; others are planned.
PQTL_SOURCES = ("gtex_fang",)
