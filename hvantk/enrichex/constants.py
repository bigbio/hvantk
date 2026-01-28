"""
Constants and default values for EnrichEx module.
"""

# Default thresholds for overlap enrichment
DEFAULT_ALPHA = 0.05
DEFAULT_CORRECTION_METHOD = "benjamini-hochberg"

# Default thresholds for burden testing
DEFAULT_MAX_AF = 0.01  # Maximum allele frequency (1%)
DEFAULT_MIN_CADD = 20.0  # Minimum CADD score
DEFAULT_MIN_GQ = 20  # Minimum genotype quality
DEFAULT_MIN_DP = 10  # Minimum depth

# Default field names for variant filtering
DEFAULT_AF_FIELD = "gnomad_af"
DEFAULT_CADD_FIELD = "cadd_phred"
DEFAULT_CONSEQUENCE_FIELD = "consequence"
DEFAULT_GENE_FIELD = "SYMBOL"

# Genotype aggregation methods
GENOTYPE_AGGREGATION_METHODS = ["hets", "homs", "chets", "homs_chets"]

# Multiple testing correction methods
CORRECTION_METHODS = ["bonferroni", "benjamini-hochberg", "none"]

# Phenotype types
PHENOTYPE_TYPES = ["binary", "continuous"]
