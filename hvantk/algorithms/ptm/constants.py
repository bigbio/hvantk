"""Algorithm-local constants for the PTM analysis algorithm.

Moved from hvantk/core/ptm_constants.py per issue #122 (clean-core principle).
"""

# Ensembl GTF -- release pinned once, in the ensembl-gene plugin. Re-exported here so
# existing importers of hvantk.algorithms.ptm.constants keep working unchanged.
from hvantk.resources.ensembl_release import (  # noqa: F401
    ENSEMBL_GTF_FILENAME,
    ENSEMBL_GTF_URL,
    ENSEMBL_RELEASE,
)

# PTM type categories (UniProt MOD_RES description prefixes)
PTM_TYPE_CATEGORIES = {
    "Phosphoserine": "phosphorylation",
    "Phosphothreonine": "phosphorylation",
    "Phosphotyrosine": "phosphorylation",
    "N6-acetyllysine": "acetylation",
    "N6-methyllysine": "methylation",
    "N6,N6-dimethyllysine": "methylation",
    "N6,N6,N6-trimethyllysine": "methylation",
    "Asymmetric dimethylarginine": "methylation",
    "Symmetric dimethylarginine": "methylation",
    "Omega-N-methylarginine": "methylation",
    "N-linked (GlcNAc...) asparagine": "n_glycosylation",
    "O-linked (GalNAc...) threonine": "o_glycosylation",
    "O-linked (GalNAc...) serine": "o_glycosylation",
    "Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in SUMO)": "sumoylation",
    "Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in SUMO2)": "sumoylation",
    "Glycyl lysine isopeptide (Lys-Gly)": "ubiquitination",
}

# Default flanking window sizes to test (in residues)
DEFAULT_FLANKING_WINDOWS = [3, 5, 7, 10]
DEFAULT_FLANKING_CODONS = 5

# Output TSV column names
PTM_OUTPUT_COLUMNS = [
    "chrom",
    "codon_start",
    "codon_end",
    "strand",
    "uniprot_id",
    "gene_symbol",
    "residue_pos",
    "amino_acid",
    "ptm_type",
    "ptm_category",
    "source_db",
    "evidence_type",
    "n_observations",
    "tissue_type",
]

# Transcript resolution method names (for logging/QC)
RESOLUTION_METHODS = ["xref_mane", "xref_any", "gene_mane"]

# ---------------------------------------------------------------------------
# Variant-annotation and LMM defaults
# (consumed by hvantk.algorithms.ptm.annotate and hvantk.algorithms.ptm.lmm).
# ---------------------------------------------------------------------------

# Proximal flanking window for gene-symbol variant annotation:
# 7 codons * 3 bp = 21 bp, applied one-sided to codon_start and codon_end.
PROXIMAL_BP = 21

# Default MAF bins (common / rare / ultra-rare) for callers that want a
# consistent set of cut-offs. Callers with source-specific thresholds
# should pass them explicitly.
DEFAULT_MAF_THRESHOLDS = {
    "common": 0.01,
    "rare": 0.001,
    "ultra_rare": 1e-5,
}

# Ordered categories for the binned-interaction LMM.
# "b0_none" is the reference (zero-expression) bin; Q1-Q4 are quartiles of
# log2(expr + 1) across positive-expression variants. Realised bin count
# may be smaller when pd.qcut collapses ties (duplicates="drop").
EXPRESSION_BIN_LABELS = ["b0_none", "b1_Q1", "b2_Q2", "b3_Q3", "b4_Q4"]

# Pseudocount for the log10(AF + eps) transformation used in LMM fits.
LOG_AF_EPSILON = 1e-8

# Minimum per-stratum sample counts for the constraint LMM
# (hvantk.algorithms.ptm.lmm.fit_constraint_lmm). Below these the per-gene estimate
# is unstable.
LMM_MIN_N_PTM = 30
LMM_MIN_N_NONPTM = 30
LMM_MIN_MIXED_GENES = 10

# Sparsity gates for the binned-interaction LMM
# (hvantk.algorithms.ptm.lmm.fit_binned_interaction_lmm).
LMM_BINNED_MIN_POS_EXPR = 100
LMM_BINNED_MIN_CELL_N = 5
