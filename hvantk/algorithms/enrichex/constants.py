"""
Constants and default values for EnrichEx module.
"""

# Default thresholds for overlap enrichment
DEFAULT_ALPHA = 0.05
DEFAULT_CORRECTION_METHOD = "benjamini-hochberg"

# Default thresholds for burden testing
DEFAULT_MAX_AF = 0.01  # Maximum allele frequency (1%)
DEFAULT_MIN_SCORE = 20.0  # Minimum prediction score (e.g., CADD, REVEL)
DEFAULT_MIN_GQ = 20  # Minimum genotype quality
DEFAULT_MIN_DP = 10  # Minimum depth

# Default field names for variant filtering
DEFAULT_AF_FIELD = "gnomad_af"
DEFAULT_SCORE_FIELD = "cadd_phred"
DEFAULT_CONSEQUENCE_FIELD = "consequence"
DEFAULT_GENE_FIELD = "SYMBOL"

# Genotype aggregation methods
GENOTYPE_AGGREGATION_METHODS = ["hets", "homs", "multi_het", "homs_multi_het"]

# Default variant class presets for stratified burden analysis.
# Consequence terms should match the values in the user's MT consequence field.
# These use common VEP Sequence Ontology terms; adjust if your annotation
# pipeline uses different labels.
VARIANT_CLASS_PRESETS = {
    "lof": {
        "consequences": [
            "stop_gained",
            "frameshift_variant",
            "splice_donor_variant",
            "splice_acceptor_variant",
        ],
    },
    "missense_constrained": {
        "consequences": ["missense_variant"],
        "min_score": 25.0,
    },
    "synonymous": {
        "consequences": ["synonymous_variant"],
    },
}

# Competitive testing defaults
DEFAULT_N_PERMUTATIONS = 10000
DEFAULT_N_LENGTH_BINS = 10

# Multiple testing correction methods
CORRECTION_METHODS = ["bonferroni", "benjamini-hochberg", "none"]

# Phenotype types
PHENOTYPE_TYPES = ["binary", "continuous"]
