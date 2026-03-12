"""
Synthetic cohort generator for burden testing validation.

This module provides functions to generate synthetic MatrixTables with
planted burden signal, useful for testing and validating the burden
analysis pipeline.

Key functions:
- generate_synthetic_burden_cohort: Create a synthetic MT with planted signal
- check_type_i_error: Verify p-value calibration under the null
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

try:  # Optional dependency
    import hail as hl
except ModuleNotFoundError as exc:  # pragma: no cover - depends on env
    hl = None  # type: ignore
    _HAIL_IMPORT_ERROR = exc
else:
    _HAIL_IMPORT_ERROR = None

logger = logging.getLogger(__name__)


def _require_hail() -> None:
    if hl is None:  # pragma: no cover - depends on env
        raise ImportError(
            "Hail is required for synthetic cohort generation. "
            "Install hvantk with the 'hail' extra."
        ) from _HAIL_IMPORT_ERROR


# ---------------------------------------------------------------------------
# Consequence categories
# ---------------------------------------------------------------------------
_CONSEQUENCE_TYPES = ["stop_gained", "missense_variant", "synonymous_variant"]
_DAMAGING_CONSEQUENCES = {"stop_gained", "missense_variant"}


def generate_synthetic_burden_cohort(
    n_cases: int = 500,
    n_controls: int = 500,
    n_genes: int = 200,
    variants_per_gene: int = 10,
    signal_gene_sets: Optional[Dict[str, List[str]]] = None,
    effect_sizes: Optional[Dict[str, float]] = None,
    baseline_carrier_rate: float = 0.05,
    seed: Optional[int] = None,
) -> Tuple["hl.MatrixTable", "hl.Table", Dict[str, List[str]]]:
    """Generate a synthetic cohort with planted burden signal.

    Creates a MatrixTable with:
    - n_cases + n_controls samples
    - n_genes * variants_per_gene variants
    - Binary genotypes (0/0 or 0/1)
    - Gene annotation (SYMBOL field)
    - Consequence annotation (stop_gained, missense_variant, synonymous_variant)

    Signal is planted by increasing the carrier rate for variants in
    signal gene sets in cases by the specified effect size (odds ratio).
    Synonymous variants have NO planted signal (calibration control).

    Parameters
    ----------
    n_cases : int
        Number of case samples.
    n_controls : int
        Number of control samples.
    n_genes : int
        Number of genes to simulate.
    variants_per_gene : int
        Number of rare variants per gene.
    signal_gene_sets : Dict[str, List[str]], optional
        Gene sets with planted signal. Keys are set names, values are lists
        of gene names. If None, creates two default sets from the first
        10 and 20 genes.
    effect_sizes : Dict[str, float], optional
        Odds ratio for each signal gene set. Keys must match signal_gene_sets.
        If None, uses OR=2.0 for all signal sets.
    baseline_carrier_rate : float
        Baseline per-variant carrier rate in controls.
    seed : int, optional
        Random seed for reproducibility.

    Returns
    -------
    Tuple[hl.MatrixTable, hl.Table, Dict[str, List[str]]]
        - mt: Synthetic cohort MatrixTable
        - phenotype_ht: Phenotype table with 'is_case' field
        - gene_sets: The gene set dict suitable for run_burden_analysis()
    """
    _require_hail()
    rng = np.random.default_rng(seed)

    n_samples = n_cases + n_controls
    n_variants = n_genes * variants_per_gene

    # ------------------------------------------------------------------
    # 1. Generate gene names
    # ------------------------------------------------------------------
    gene_names = [f"GENE_{i}" for i in range(n_genes)]

    # ------------------------------------------------------------------
    # 2. Default signal gene sets if not provided
    # ------------------------------------------------------------------
    if signal_gene_sets is None:
        signal_gene_sets = {
            "signal_set_small": gene_names[:10],
            "signal_set_large": gene_names[:20],
        }
    if effect_sizes is None:
        effect_sizes = {gs_name: 2.0 for gs_name in signal_gene_sets}

    # Validate that signal gene names exist
    gene_name_set = set(gene_names)
    for gs_name, genes in signal_gene_sets.items():
        unknown = set(genes) - gene_name_set
        if unknown:
            raise ValueError(
                f"Signal gene set '{gs_name}' contains unknown genes: {unknown}. "
                f"Valid gene names are GENE_0 .. GENE_{n_genes - 1}."
            )

    # Build lookup: gene_name -> gene_index for fast access
    gene_name_to_idx = {g: i for i, g in enumerate(gene_names)}

    # ------------------------------------------------------------------
    # 3. Assign consequences per variant (deterministic, cycled)
    # ------------------------------------------------------------------
    consequences: List[str] = []
    for v in range(n_variants):
        consequences.append(_CONSEQUENCE_TYPES[v % len(_CONSEQUENCE_TYPES)])

    # ------------------------------------------------------------------
    # 4. Build per-variant SYMBOL list
    # ------------------------------------------------------------------
    symbols: List[str] = []
    for g_idx in range(n_genes):
        symbols.extend([gene_names[g_idx]] * variants_per_gene)

    # ------------------------------------------------------------------
    # 5. Generate genotype matrix
    # ------------------------------------------------------------------
    # Start with baseline carrier probability for everyone
    carrier_probs = np.full((n_variants, n_samples), baseline_carrier_rate)

    # Increase carrier rate for damaging variants in signal genes in cases
    for gs_name, genes in signal_gene_sets.items():
        or_val = effect_sizes.get(gs_name, 2.0)
        for gene in genes:
            gene_idx = gene_name_to_idx[gene]
            var_start = gene_idx * variants_per_gene
            for v in range(variants_per_gene):
                var_idx = var_start + v
                consequence = consequences[var_idx]
                if consequence in _DAMAGING_CONSEQUENCES:
                    # Increase rate for cases (first n_cases columns)
                    carrier_probs[var_idx, :n_cases] *= or_val

    # Cap at a reasonable maximum
    carrier_probs = np.clip(carrier_probs, 0, 0.5)

    # Draw carrier status
    is_carrier = (rng.random((n_variants, n_samples)) < carrier_probs).astype(int)

    logger.info(
        "Generated synthetic genotype matrix: %d variants x %d samples",
        n_variants,
        n_samples,
    )

    # ------------------------------------------------------------------
    # 6. Build Hail MatrixTable
    # ------------------------------------------------------------------
    mt = hl.utils.range_matrix_table(n_variants, n_samples)
    mt = mt.key_cols_by(s=hl.str(mt.col_idx))

    # Row annotations
    symbols_literal = hl.literal(symbols)
    consequences_literal = hl.literal(consequences)
    mt = mt.annotate_rows(
        SYMBOL=symbols_literal[mt.row_idx],
        consequence=consequences_literal[mt.row_idx],
        gnomad_af=hl.literal(baseline_carrier_rate),
    )

    # Entry annotations — genotype from the pre-computed numpy matrix
    gt_list = is_carrier.tolist()
    gt_literal = hl.literal(gt_list)
    mt = mt.annotate_entries(
        GT=hl.call(0, gt_literal[mt.row_idx][mt.col_idx]),
    )

    # ------------------------------------------------------------------
    # 7. Build phenotype Table
    # ------------------------------------------------------------------
    pheno_rows = [hl.struct(s=str(i), is_case=(i < n_cases)) for i in range(n_samples)]
    phenotype_ht = hl.Table.parallelize(
        pheno_rows,
        schema=hl.tstruct(s=hl.tstr, is_case=hl.tbool),
    ).key_by("s")

    logger.info(
        "Synthetic cohort: %d cases, %d controls, %d genes, " "%d signal gene sets",
        n_cases,
        n_controls,
        n_genes,
        len(signal_gene_sets),
    )

    return mt, phenotype_ht, signal_gene_sets


# ---------------------------------------------------------------------------
# Type-I error calibration check
# ---------------------------------------------------------------------------


def check_type_i_error(
    p_values: List[float],
    alpha: float = 0.05,
    tolerance: float = 2.0,
) -> Dict[str, Any]:
    """Check if p-values are uniformly distributed (type-I error calibration).

    Under the null hypothesis, p-values should be uniformly distributed on
    [0, 1].  This function checks the rejection rate at the given *alpha*
    and runs a Kolmogorov-Smirnov test for uniformity.

    Parameters
    ----------
    p_values : List[float]
        P-values from burden tests under the null hypothesis.
    alpha : float
        Significance level.
    tolerance : float
        Allowed fold-deviation from expected rejection rate.  The test is
        considered calibrated if the observed rejection rate is within
        ``[expected / tolerance, expected * tolerance]``.

    Returns
    -------
    Dict[str, Any]
        Dictionary with keys:

        - ``n_tests`` : int -- number of p-values
        - ``n_rejected`` : int -- number of p-values <= alpha
        - ``rejection_rate`` : float -- fraction rejected
        - ``expected_rate`` : float -- alpha
        - ``calibrated`` : bool -- whether rejection rate is within tolerance
        - ``ks_pvalue`` : float -- Kolmogorov-Smirnov test p-value
          (null = uniform)
    """
    from scipy import stats as sp_stats

    p_arr = np.asarray(p_values, dtype=float)
    n_tests = len(p_arr)

    if n_tests == 0:
        return {
            "n_tests": 0,
            "n_rejected": 0,
            "rejection_rate": 0.0,
            "expected_rate": alpha,
            "calibrated": False,
            "ks_pvalue": 0.0,
        }

    n_rejected = int(np.sum(p_arr <= alpha))
    rejection_rate = n_rejected / n_tests
    expected_rate = alpha

    # Check if rejection rate is within tolerance band
    lower_bound = expected_rate / tolerance
    upper_bound = expected_rate * tolerance
    rate_ok = lower_bound <= rejection_rate <= upper_bound

    # Kolmogorov-Smirnov test against uniform(0, 1)
    ks_stat, ks_pvalue = sp_stats.kstest(p_arr, "uniform")

    # Calibrated if both the rejection rate is within tolerance AND
    # the KS test does not reject uniformity at a generous threshold
    calibrated = rate_ok and ks_pvalue > 0.001

    return {
        "n_tests": n_tests,
        "n_rejected": n_rejected,
        "rejection_rate": rejection_rate,
        "expected_rate": expected_rate,
        "calibrated": calibrated,
        "ks_pvalue": float(ks_pvalue),
    }
