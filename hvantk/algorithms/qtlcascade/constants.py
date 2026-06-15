"""Algorithm-local constants for QTL cascade analysis.

Moved from hvantk/core/qtl_constants.py per issue #122 (clean-core principle).

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
    "discordant": "Discordant (eQTL ≠ pQTL direction)",
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
# GWAS × eQTL colocalization (gwas_coloc.py)
# ---------------------------------------------------------------------------
# Trait-specific prior effect-size variances W (Wakefield 2009). coloc's
# conventional defaults: case-control GWAS sd(logOR)=0.2 -> W=0.04;
# quantitative cis-eQTL sd=0.15 -> W=0.0225.
DEFAULT_GWAS_W_CC = 0.04
DEFAULT_EQTL_W_QUANT = 0.0225

# Minimum shared variants in a region for a gene to be coloc-tested.
DEFAULT_COLOC_MIN_SNPS = 20

# Public summary-statistic sources (remote-tabix; no bulk download).
# FinnGen R10 GWAS: contig has no 'chr'; cols chrom pos ref alt rsids
#   nearest_genes pval mlogp beta sebeta af*.
FINNGEN_R10_SUMSTATS_URL = (
    "https://storage.googleapis.com/finngen-public-data-r10/"
    "summary_stats/finngen_R10_{endpoint}.gz"
)
# eQTL Catalogue cis-eQTL all-variant sumstats; cols gene_id chrom pos ref alt
#   variant ma_samples maf pvalue beta se. Default study QTS000015 = GTEx.
EQTL_CATALOGUE_SUMSTATS_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"
    "{study}/{dataset}/{dataset}.all.tsv.gz"
)
EQTL_CATALOGUE_DEFAULT_STUDY = "QTS000015"

# ---------------------------------------------------------------------------
# Fine-mapping LD reference (finemap.py) — 1000 Genomes high-coverage GRCh38.
# Used to build a EUR LD matrix for SuSiE-RSS + coloc.susie. Optional layer.
# ---------------------------------------------------------------------------
KG_PHASED_VCF_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000G_2504_high_coverage/working/20201028_3202_phased/"
    "CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chrom}."
    "filtered.shapeit2-duohmm-phased.vcf.gz"
)
KG_PANEL_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt"
)
# Default fine-map super-population for the LD reference.
DEFAULT_FINEMAP_SUPERPOP = "EUR"

# ---------------------------------------------------------------------------
# SuSiE-RSS + coloc.susie (susie.py) — pure-Python fine-mapping (issue #193).
# Ports susieR::susie_rss(z, R, n, L) and coloc::coloc.susie; the priors below
# are the coloc package's coloc.susie defaults (p12 differs from the
# single-variant ABF p12 = 1e-5 above).
# ---------------------------------------------------------------------------
DEFAULT_SUSIE_L = 10              # max single effects (susie_rss L)
DEFAULT_SUSIE_COVERAGE = 0.95     # credible-set coverage
DEFAULT_SUSIE_MIN_ABS_CORR = 0.5  # credible-set purity filter (min |r|)
DEFAULT_SUSIE_MAX_ITER = 100      # IBSS iterations
DEFAULT_COLOC_SUSIE_P1 = 1e-4     # coloc.susie P(causal for trait 1 only)
DEFAULT_COLOC_SUSIE_P2 = 1e-4     # coloc.susie P(causal for trait 2 only)
DEFAULT_COLOC_SUSIE_P12 = 5e-6    # coloc.susie P(shared causal variant)
