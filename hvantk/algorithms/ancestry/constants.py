"""Default parameters and constants for ancestry inference.

This module defines default values for all configurable parameters in the
ancestry inference pipeline, as well as standard column names and visualization
settings.
"""

from typing import Dict, Tuple

# =============================================================================
# Variant Filtering Defaults
# =============================================================================

#: Minimum allele frequency threshold for variant inclusion
DEFAULT_MIN_AF: float = 0.01

#: Maximum allele frequency threshold for variant inclusion
DEFAULT_MAX_AF: float = 0.99

#: Minimum call rate threshold for variant inclusion
DEFAULT_MIN_CALL_RATE: float = 0.98

#: Hardy-Weinberg equilibrium p-value threshold (variants below are excluded)
DEFAULT_HWE_P: float = 1e-6

# =============================================================================
# LD Pruning Defaults
# =============================================================================

#: LD r-squared threshold for pruning correlated variants
DEFAULT_LD_R2: float = 0.2

#: Window size in base pairs for LD pruning
DEFAULT_LD_WINDOW: int = 500_000  # 500 Kb

# =============================================================================
# PCA Defaults
# =============================================================================

#: Number of principal components to compute
DEFAULT_N_PCS: int = 20

#: Number of principal components to use for classification
DEFAULT_N_PCS_CLASSIFY: int = 10

# =============================================================================
# Classification Defaults
# =============================================================================

#: Number of trees in Random Forest classifier
DEFAULT_N_ESTIMATORS: int = 100

#: Minimum probability threshold for ancestry assignment
#: Samples with max probability below this are labeled "unassigned"
DEFAULT_MIN_PROB: float = 0.75

#: Random seed for reproducibility
DEFAULT_RANDOM_SEED: int = 42

# =============================================================================
# Validation Defaults
# =============================================================================

#: Minimum number of samples per population for training
MIN_SAMPLES_PER_POP: int = 10

#: Default number of cross-validation folds
DEFAULT_N_CV_FOLDS: int = 5

# =============================================================================
# Column Names
# =============================================================================

#: Default column name for ancestry labels in reference MatrixTable
ANCESTRY_COL: str = "ancestry"

#: Column name for predicted ancestry in output
PREDICTED_ANCESTRY_COL: str = "predicted_ancestry"

#: Column name for ancestry prediction probability
ANCESTRY_PROB_COL: str = "ancestry_probability"

#: Internal column name for tracking sample source (query vs reference)
SOURCE_COL: str = "_ancestry_source"

#: Internal column name for known ancestry labels
KNOWN_ANCESTRY_COL: str = "_known_ancestry"

# =============================================================================
# Variant Thresholds
# =============================================================================

#: Minimum number of shared variants required for analysis
MIN_SHARED_VARIANTS: int = 10_000

#: Warning threshold for shared variants (below triggers a warning)
WARN_SHARED_VARIANTS: int = 50_000

# =============================================================================
# 1000 Genomes Phase 3 Population Definitions
# =============================================================================

#: Super-population to sub-population mapping for 1000 Genomes Phase 3
#: Source: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
#:         integrated_call_samples_v3.20130502.ALL.panel
SUPERPOP_TO_SUBPOPS: Dict[str, Tuple[str, ...]] = {
    # African (AFR) - 7 populations
    "AFR": (
        "ACB",  # African Caribbean in Barbados
        "ASW",  # African Ancestry in Southwest US
        "ESN",  # Esan in Nigeria
        "GWD",  # Gambian in Western Division, The Gambia
        "LWK",  # Luhya in Webuye, Kenya
        "MSL",  # Mende in Sierra Leone
        "YRI",  # Yoruba in Ibadan, Nigeria
    ),
    # Ad Mixed American (AMR) - 4 populations
    "AMR": (
        "CLM",  # Colombian in Medellin, Colombia
        "MXL",  # Mexican Ancestry in Los Angeles, California
        "PEL",  # Peruvian in Lima, Peru
        "PUR",  # Puerto Rican in Puerto Rico
    ),
    # East Asian (EAS) - 5 populations
    "EAS": (
        "CDX",  # Chinese Dai in Xishuangbanna, China
        "CHB",  # Han Chinese in Beijing, China
        "CHS",  # Southern Han Chinese, China
        "JPT",  # Japanese in Tokyo, Japan
        "KHV",  # Kinh in Ho Chi Minh City, Vietnam
    ),
    # European (EUR) - 5 populations
    "EUR": (
        "CEU",  # Utah residents with Northern/Western European ancestry
        "FIN",  # Finnish in Finland
        "GBR",  # British in England and Scotland
        "IBS",  # Iberian populations in Spain
        "TSI",  # Toscani in Italia
    ),
    # South Asian (SAS) - 5 populations
    "SAS": (
        "BEB",  # Bengali in Bangladesh
        "GIH",  # Gujarati Indians in Houston, Texas
        "ITU",  # Indian Telugu in the UK
        "PJL",  # Punjabi in Lahore, Pakistan
        "STU",  # Sri Lankan Tamil in the UK
    ),
}

#: Reverse mapping: sub-population to super-population
SUBPOP_TO_SUPERPOP: Dict[str, str] = {
    subpop: superpop
    for superpop, subpops in SUPERPOP_TO_SUBPOPS.items()
    for subpop in subpops
}

#: All 1000 Genomes Phase 3 sub-population codes (26 populations)
ALL_1KG_SUBPOPS: Tuple[str, ...] = tuple(sorted(SUBPOP_TO_SUPERPOP.keys()))

#: All 1000 Genomes Phase 3 super-population codes (5 populations)
ALL_1KG_SUPERPOPS: Tuple[str, ...] = ("AFR", "AMR", "EAS", "EUR", "SAS")

#: Population full names for documentation and reports
POPULATION_NAMES: Dict[str, str] = {
    # Super-populations
    "AFR": "African",
    "AMR": "Ad Mixed American",
    "EAS": "East Asian",
    "EUR": "European",
    "SAS": "South Asian",
    # Sub-populations
    "ACB": "African Caribbean in Barbados",
    "ASW": "African Ancestry in Southwest US",
    "ESN": "Esan in Nigeria",
    "GWD": "Gambian in Western Division, The Gambia",
    "LWK": "Luhya in Webuye, Kenya",
    "MSL": "Mende in Sierra Leone",
    "YRI": "Yoruba in Ibadan, Nigeria",
    "CLM": "Colombian in Medellin, Colombia",
    "MXL": "Mexican Ancestry in Los Angeles, California",
    "PEL": "Peruvian in Lima, Peru",
    "PUR": "Puerto Rican in Puerto Rico",
    "CDX": "Chinese Dai in Xishuangbanna, China",
    "CHB": "Han Chinese in Beijing, China",
    "CHS": "Southern Han Chinese, China",
    "JPT": "Japanese in Tokyo, Japan",
    "KHV": "Kinh in Ho Chi Minh City, Vietnam",
    "CEU": "Utah residents with Northern/Western European ancestry",
    "FIN": "Finnish in Finland",
    "GBR": "British in England and Scotland",
    "IBS": "Iberian populations in Spain",
    "TSI": "Toscani in Italia",
    "BEB": "Bengali in Bangladesh",
    "GIH": "Gujarati Indians in Houston, Texas",
    "ITU": "Indian Telugu in the UK",
    "PJL": "Punjabi in Lahore, Pakistan",
    "STU": "Sri Lankan Tamil in the UK",
    # Special
    "unassigned": "Unassigned (below probability threshold)",
}

# =============================================================================
# Population Color Palette
# =============================================================================

#: Default color palette for super-population visualization
#: Based on 1000 Genomes Phase 3 super-populations
SUPERPOP_COLORS: Dict[str, str] = {
    "AFR": "#ff7f0e",  # Orange
    "AMR": "#9467bd",  # Purple
    "EAS": "#2ca02c",  # Green
    "EUR": "#1f77b4",  # Blue
    "SAS": "#d62728",  # Red
    # Special labels
    "unassigned": "#7f7f7f",  # Gray
}

#: Extended color palette including sub-populations
#: Sub-populations use shades of their super-population color
POPULATION_COLORS: Dict[str, str] = {
    # Super-populations (same as SUPERPOP_COLORS)
    "AFR": "#ff7f0e",  # Orange
    "AMR": "#9467bd",  # Purple
    "EAS": "#2ca02c",  # Green
    "EUR": "#1f77b4",  # Blue
    "SAS": "#d62728",  # Red
    # AFR sub-populations (orange shades)
    "ACB": "#ff9f40",
    "ASW": "#ffbf70",
    "ESN": "#e56b00",
    "GWD": "#cc5f00",
    "LWK": "#ff8c1a",
    "MSL": "#ffa64d",
    "YRI": "#b35300",
    # AMR sub-populations (purple shades)
    "CLM": "#a683c9",
    "MXL": "#b89fd6",
    "PEL": "#7b4fa2",
    "PUR": "#c9b8e3",
    # EAS sub-populations (green shades)
    "CDX": "#45b745",
    "CHB": "#1e8c1e",
    "CHS": "#5fd35f",
    "JPT": "#33a633",
    "KHV": "#79e079",
    # EUR sub-populations (blue shades)
    "CEU": "#3498db",
    "FIN": "#5faee3",
    "GBR": "#1a6aa5",
    "IBS": "#89c4ec",
    "TSI": "#0d4f78",
    # SAS sub-populations (red shades)
    "BEB": "#e94848",
    "GIH": "#a11f1f",
    "ITU": "#f47272",
    "PJL": "#c12929",
    "STU": "#ff9c9c",
    # Special labels
    "unassigned": "#7f7f7f",  # Gray
    "UNKNOWN": "#bcbd22",  # Olive
}

# =============================================================================
# Legacy Constants (for backward compatibility)
# =============================================================================

#: Alias for backward compatibility
KNOWN_1KG_SUPERPOPS: Tuple[str, ...] = ALL_1KG_SUPERPOPS

#: Known population codes from HapMap Phase 3
#: Note: CHD (Chinese in Metropolitan Denver) was in HapMap but not in 1KG Phase 3
KNOWN_HAPMAP_POPS: Tuple[str, ...] = (
    "ASW",  # African ancestry in Southwest USA
    "CEU",  # Utah residents with European ancestry
    "CHB",  # Han Chinese in Beijing
    "CHD",  # Chinese in Metropolitan Denver (HapMap only)
    "GIH",  # Gujarati Indians in Houston
    "JPT",  # Japanese in Tokyo
    "LWK",  # Luhya in Webuye, Kenya
    "MEX",  # Mexican ancestry in Los Angeles (HapMap code, MXL in 1KG)
    "MKK",  # Maasai in Kinyawa, Kenya (HapMap only)
    "TSI",  # Toscani in Italia
    "YRI",  # Yoruba in Ibadan, Nigeria
)
