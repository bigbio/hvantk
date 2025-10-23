"""
A set of constants used throughout the package.

"""

# Genome references
HG38_GENOME_REFERENCE = "GRCh38"
HG37_GENOME_REFERENCE = "GRCh37"

# File extensions
GVCF_EXTENSION = ".g.vcf.gz"
GVCF_EXTENSION_TBI = f"{GVCF_EXTENSION}.tbi"
VCF_EXTENSION = ".vcf.gz"
VCF_EXTENSION_TBI = f"{VCF_EXTENSION}.tbi"
VDS_EXTENSION = ".vds"

# Entry fields
GT_FIELD = "GT"
AD_FIELD = "AD"
DP_FIELD = "DP"
GQ_FIELD = "GQ"
PL_FIELD = "PL"
PID_FIELD = "PID"
SB_FIELD = "SB"
MIN_DP_FIELD = "MIN_DP"
ADJ_GT_FIELD = "adj"
