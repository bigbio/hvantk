"""Plugin-local constants for the clinvar skill.

Moved from hvantk/core/constants.py per issue #120 (clean-core principle).
These pathogenicity-label vocabularies describe ClinVar's CLNSIG
classification. Algorithms must NOT import them from here: the dependency
guards forbid ``algorithms/`` importing from ``skills/``, so the algorithm
modules (algorithms/ptm, algorithms/psroc, algorithms/annotation) keep
their own temporary local copies. Decoupling those copies via an injected
vocabulary is tracked by issue #133 (and follow-up #145).
"""

# ClinVar FTP base URLs
CLINVAR_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38"
CLINVAR_FTP_BASE_GRCh37 = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37"

# ClinVar CLNSIG classification labels.
CLINVAR_PATHOGENIC_LABELS = [
    "Pathogenic/Likely_pathogenic",
    "Likely_pathogenic",
    "Pathogenic",
]

CLINVAR_BENIGN_LABELS = [
    "Benign/Likely_benign",
    "Likely_benign",
    "Benign",
]
