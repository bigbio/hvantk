"""
Colocalization analysis using Approximate Bayes Factors (ABF).

Tests whether eQTL and pQTL association signals in a genomic region share
the same causal variant, distinguishing true signal propagation
(variant -> mRNA -> protein) from LD confounding.

Five hypotheses per region:
    H0 — no association with either trait
    H1 — association with eQTL only
    H2 — association with pQTL only
    H3 — both associated, **different** causal variants (LD artifact)
    H4 — both associated, **shared** causal variant (true cascade)

Implementation follows the R coloc package
(github.com/chr1swallace/coloc, core logic in ``R/claudia.R``).

References
----------
- Giambartolomei et al. (2014) PLoS Genet 10(5):e1004383 — original method
- Wakefield (2009) Am J Hum Genet 84(1):60-71 — ABF formula
- Pullin & Wallace (2025) PLoS Genet 21(5):e1011697 — v6 extension
"""

import logging
from typing import Optional

import numpy as np
import pandas as pd

from hvantk.qtlcascade.constants import (
    DEFAULT_COLOC_H4_THRESHOLD,
    DEFAULT_COLOC_P1,
    DEFAULT_COLOC_P2,
    DEFAULT_COLOC_P12,
    DEFAULT_COLOC_W,
    DEFAULT_COLOC_WINDOW_KB,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Core ABF computation
# ---------------------------------------------------------------------------


def compute_log_abf(
    beta: np.ndarray,
    se: np.ndarray,
    W: float = DEFAULT_COLOC_W,
) -> np.ndarray:
    """Per-variant log Approximate Bayes Factor (Wakefield 2009, Eq. 2).

    Parameters
    ----------
    beta : array
        Effect sizes (regression coefficients).
    se : array
        Standard errors of the effect sizes.
    W : float
        Prior variance on the true effect size.

    Returns
    -------
    array
        Log ABF for each variant.
    """
    r = W / (W + se**2)
    z = beta / se
    return 0.5 * (np.log(1.0 - r) + r * z**2)


def _logsumexp(x: np.ndarray) -> float:
    """Numerically stable log-sum-exp."""
    m = np.max(x)
    return float(m + np.log(np.sum(np.exp(x - m))))


def coloc_abf(
    eqtl_beta: np.ndarray,
    eqtl_se: np.ndarray,
    pqtl_beta: np.ndarray,
    pqtl_se: np.ndarray,
    p1: float = DEFAULT_COLOC_P1,
    p2: float = DEFAULT_COLOC_P2,
    p12: float = DEFAULT_COLOC_P12,
    W: float = DEFAULT_COLOC_W,
) -> dict:
    """Run coloc ABF for a single genomic region.

    Parameters
    ----------
    eqtl_beta, eqtl_se : array
        Summary statistics for the eQTL trait.
    pqtl_beta, pqtl_se : array
        Summary statistics for the pQTL trait.
    p1, p2, p12 : float
        Prior probabilities (see module docstring).
    W : float
        Prior variance on effect size.

    Returns
    -------
    dict
        H0–H4 posterior probabilities, ``n_variants``, and
        ``lead_snp_h4_idx`` (variant contributing most to H4).
    """
    n = len(eqtl_beta)
    if n == 0:
        return {
            "H0": 1.0, "H1": 0.0, "H2": 0.0, "H3": 0.0, "H4": 0.0,
            "n_variants": 0, "lead_snp_h4_idx": -1,
        }

    log_abf1 = compute_log_abf(eqtl_beta, eqtl_se, W)
    log_abf2 = compute_log_abf(pqtl_beta, pqtl_se, W)

    sum_abf1 = _logsumexp(log_abf1)
    sum_abf2 = _logsumexp(log_abf2)

    # Hypothesis log-likelihoods (unnormalised)
    log_h0 = 0.0
    log_h1 = np.log(p1) + sum_abf1
    log_h2 = np.log(p2) + sum_abf2
    log_h3 = np.log(p1) + np.log(p2) + sum_abf1 + sum_abf2
    log_abf_both = log_abf1 + log_abf2
    log_h4 = np.log(p12) + _logsumexp(log_abf_both)

    # Posterior via softmax
    all_log = np.array([log_h0, log_h1, log_h2, log_h3, log_h4])
    posteriors = np.exp(all_log - np.max(all_log))
    posteriors /= posteriors.sum()

    return {
        "H0": float(posteriors[0]),
        "H1": float(posteriors[1]),
        "H2": float(posteriors[2]),
        "H3": float(posteriors[3]),
        "H4": float(posteriors[4]),
        "n_variants": n,
        "lead_snp_h4_idx": int(np.argmax(log_abf_both)),
    }


# ---------------------------------------------------------------------------
# Per-gene coloc driver (Hail + NumPy hybrid)
# ---------------------------------------------------------------------------


def run_coloc_per_gene(
    eqtl_allpairs_ht_path: str,
    pqtl_allpairs_ht_path: str,
    cascade_genes: list,
    tissue: Optional[str] = None,
    window_kb: int = DEFAULT_COLOC_WINDOW_KB,
    p1: float = DEFAULT_COLOC_P1,
    p2: float = DEFAULT_COLOC_P2,
    p12: float = DEFAULT_COLOC_P12,
    W: float = DEFAULT_COLOC_W,
) -> pd.DataFrame:
    """Run coloc for all cascade genes.

    Uses Hail for bulk data extraction (one Spark job) and NumPy for
    per-gene ABF computation.

    Parameters
    ----------
    eqtl_allpairs_ht_path : str
        Path to allpairs eQTL Hail Table.
    pqtl_allpairs_ht_path : str
        Path to allpairs pQTL Hail Table.
    cascade_genes : list[str]
        Gene IDs with both eQTL and pQTL evidence.
    tissue : str, optional
        Filter allpairs tables to this tissue.
    window_kb : int
        Window (±kb) around lead variant for regional extraction.
    p1, p2, p12, W : float
        Coloc prior parameters.

    Returns
    -------
    pd.DataFrame
        Columns: gene_id, tissue, H0–H4, n_variants.
    """
    import hail as hl

    result_cols = ["gene_id", "tissue", "H0", "H1", "H2", "H3", "H4", "n_variants"]
    empty = pd.DataFrame(columns=result_cols)

    if not cascade_genes:
        return empty

    eqtl_ht = hl.read_table(eqtl_allpairs_ht_path)
    pqtl_ht = hl.read_table(pqtl_allpairs_ht_path)

    # Filter to cascade genes (single Spark filter)
    gene_set = hl.literal(set(cascade_genes))
    eqtl_ht = eqtl_ht.filter(gene_set.contains(eqtl_ht.gene_id))
    pqtl_ht = pqtl_ht.filter(gene_set.contains(pqtl_ht.gene_id))

    if tissue:
        if "tissue" in list(eqtl_ht.row):
            eqtl_ht = eqtl_ht.filter(eqtl_ht.tissue == tissue)
        if "tissue" in list(pqtl_ht.row):
            pqtl_ht = pqtl_ht.filter(pqtl_ht.tissue == tissue)

    # Select fields for join; annotate position for windowing
    eqtl_sel = eqtl_ht.select(
        eqtl_beta=eqtl_ht.beta,
        eqtl_se=eqtl_ht.se,
        eqtl_p=eqtl_ht.p_value,
        position=eqtl_ht.locus.position,
    )
    pqtl_sel = pqtl_ht.select(
        pqtl_beta=pqtl_ht.beta,
        pqtl_se=pqtl_ht.se,
    )

    # Inner join on (locus, alleles, gene_id) — single Spark job
    joined = eqtl_sel.join(pqtl_sel, how="inner")

    # Flatten for pandas export
    joined = joined.key_by()
    joined = joined.annotate(
        contig=joined.locus.contig,
        pos=joined.locus.position,
    )
    joined = joined.select(
        "gene_id", "pos", "eqtl_beta", "eqtl_se", "eqtl_p",
        "pqtl_beta", "pqtl_se",
    )

    logger.info(
        "Exporting joined allpairs for coloc (%d cascade genes)", len(cascade_genes)
    )
    df = joined.to_pandas()

    if df.empty:
        logger.warning("No overlapping variants found between allpairs tables")
        return empty

    # Per-gene coloc with regional windowing (pure Python)
    window_bp = window_kb * 1000
    results = []

    for gene_id, group in df.groupby("gene_id"):
        # Window around lead eQTL variant
        lead_pos = group.loc[group["eqtl_p"].idxmin(), "pos"]
        region = group[np.abs(group["pos"] - lead_pos) <= window_bp]

        if len(region) < 2:
            continue

        row = coloc_abf(
            eqtl_beta=region["eqtl_beta"].values,
            eqtl_se=region["eqtl_se"].values,
            pqtl_beta=region["pqtl_beta"].values,
            pqtl_se=region["pqtl_se"].values,
            p1=p1, p2=p2, p12=p12, W=W,
        )
        row["gene_id"] = gene_id
        row["tissue"] = tissue or "unknown"
        results.append(row)

    if not results:
        logger.warning("Coloc produced no results (check variant overlap)")
        return empty

    logger.info("Coloc completed for %d / %d genes", len(results), len(cascade_genes))
    return pd.DataFrame(results)[result_cols]
