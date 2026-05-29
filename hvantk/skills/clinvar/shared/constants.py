"""Plugin-local constants for the clinvar skill.

Moved from hvantk/core/constants.py per issue #120 (clean-core principle).
These pathogenicity-label vocabularies are also consumed by three
algorithms (algorithms/ptm, algorithms/psroc, algorithms/annotation)
which import them directly from this module; the dependency is
intentional (those algorithms interpret ClinVar's vocabulary by design).
"""

# ClinVar FTP base URLs
CLINVAR_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38"
CLINVAR_FTP_BASE_GRCh37 = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37"

# ClinVar classification labels (vocabulary used by algorithms/ptm,
# algorithms/psroc, algorithms/annotation).
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
