"""
Wilcoxon rank-sum marker gene detection for scRNA-seq data.

Implements vectorised one-vs-rest Wilcoxon rank-sum tests (Mann–Whitney U)
with tie correction and multiple-testing correction.  Designed as the
statistical backend for ``hvantk expression markers --method wilcoxon``.

The implementation follows the two-phase approach used by Seurat's
``FindAllMarkers``:

1. **Pre-filter** candidate genes (fold-change + fraction-expressed cutoffs).
2. **Rank-sum test** all candidates in a single vectorised pass.

All functions in this module are pure numpy/scipy — no Hail dependency.

References
----------
- Pullin JM, McCarthy DJ. "A comparison of marker gene selection methods for
  single-cell RNA sequencing data". *Genome Biology* 25, 56 (2024).
- Luecken MD, Theis FJ. "Current best practices in single-cell RNA-seq
  analysis: a tutorial". *Mol Syst Biol* 15, e8746 (2019).
- Squair JW et al. "Confronting false discoveries in single-cell
  differential expression". *Nat Commun* 12, 5692 (2021).
- Stuart T et al. "Comprehensive Integration of Single-Cell Data". *Cell*
  177, 1888–1902 (2019). [Seurat v3]
- Wolf FA et al. "SCANPY: large-scale single-cell gene expression data
  analysis". *Genome Biol* 19, 15 (2018).
- Benjamini Y, Hochberg Y. "Controlling the false discovery rate: a
  practical and powerful approach to multiple testing". *JRSS B* 57,
  289–300 (1995).
- Mann HB, Whitney DR. "On a test of whether one of two random variables is
  stochastically larger than the other". *Ann Math Stat* 18, 50–60 (1947).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Dict, Optional, Set, Tuple

import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata

from hvantk.enrichex.correction import apply_correction
from hvantk.utils.gene_sets import GeneSet, GeneSetCollection

logger = logging.getLogger(__name__)

__all__ = [
    "WilcoxonParams",
    "rank_genes_groups",
    "results_to_gene_set_collection",
]


# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------


@dataclass
class WilcoxonParams:
    """Parameters for Wilcoxon rank-sum marker gene detection.

    Attributes
    ----------
    min_fold_change : float
        Pre-filter: minimum fold-change (group mean / rest mean) for a gene
        to be tested.
    min_fraction_expressed : float
        Pre-filter: minimum fraction of cells in the focal group expressing
        the gene (expression > 0).
    max_candidates : int
        Maximum number of candidate genes to carry into Phase 2 (rank-sum).
        Genes are prioritised by fold-change.
    top_n : int
        Maximum significant markers to return per group.
    correction_method : str
        Multiple-testing correction method passed to
        :func:`~hvantk.enrichex.correction.apply_correction`.
    alpha : float
        Significance threshold on adjusted p-values.
    tie_correction : bool
        Apply tie-correction to the variance of the U statistic.
    """

    min_fold_change: float = 1.5
    min_fraction_expressed: float = 0.1
    max_candidates: int = 5000
    top_n: int = 200
    correction_method: str = "benjamini-hochberg"
    alpha: float = 0.05
    tie_correction: bool = True


# ---------------------------------------------------------------------------
# Core statistical functions
# ---------------------------------------------------------------------------


def _compute_rank_matrix(expression: np.ndarray) -> np.ndarray:
    """Rank expression values per gene across all cells.

    Parameters
    ----------
    expression : np.ndarray
        Dense expression matrix of shape ``(n_cells, n_genes)``.

    Returns
    -------
    np.ndarray
        Rank matrix of the same shape, using average tie-breaking.
    """
    n_cells, n_genes = expression.shape
    ranks = np.empty_like(expression, dtype=np.float64)
    for j in range(n_genes):
        ranks[:, j] = rankdata(expression[:, j], method="average")
    return ranks


def _compute_tie_correction(expression: np.ndarray, n_total: int) -> np.ndarray:
    """Compute per-gene tie-correction factors.

    The correction factor is ``1 - sum(t**3 - t) / (N**3 - N)`` where *t*
    is the number of observations sharing a rank and *N* is the total number
    of observations.  A value of 1.0 means no ties; 0.0 means all values
    are identical.

    Parameters
    ----------
    expression : np.ndarray
        Dense expression matrix ``(n_cells, n_genes)``.
    n_total : int
        Total number of cells (rows).

    Returns
    -------
    np.ndarray
        1-D array of length ``n_genes`` with tie-correction factors ∈ [0, 1].
    """
    n_genes = expression.shape[1]
    corrections = np.ones(n_genes, dtype=np.float64)
    denom = float(n_total) ** 3 - float(n_total)
    if denom == 0:
        return corrections  # n_total <= 1 — no correction possible

    for j in range(n_genes):
        _, counts = np.unique(expression[:, j], return_counts=True)
        tie_sum = np.sum(counts.astype(np.float64) ** 3 - counts.astype(np.float64))
        corrections[j] = 1.0 - tie_sum / denom

    return corrections


def _wilcoxon_one_vs_rest(
    ranks: np.ndarray,
    group_mask: np.ndarray,
    n_total: int,
    tie_correction: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Vectorised Wilcoxon rank-sum test for one group vs the rest.

    Parameters
    ----------
    ranks : np.ndarray
        Rank matrix ``(n_cells, n_genes)`` (pre-computed).
    group_mask : np.ndarray
        Boolean mask ``(n_cells,)`` selecting the focal group.
    n_total : int
        Total number of cells.
    tie_correction : np.ndarray
        Per-gene tie-correction factors ``(n_genes,)``.

    Returns
    -------
    u_stats : np.ndarray
        Mann–Whitney U statistics ``(n_genes,)``.
    z_scores : np.ndarray
        Normal-approximation z-scores ``(n_genes,)``.
    p_values : np.ndarray
        Two-sided p-values ``(n_genes,)``.
    """
    n1 = int(group_mask.sum())
    n2 = n_total - n1

    if n1 == 0 or n2 == 0:
        n_genes = ranks.shape[1]
        return (
            np.full(n_genes, np.nan),
            np.full(n_genes, 0.0),
            np.full(n_genes, 1.0),
        )

    # Rank sums for the focal group (vectorised across genes)
    R1 = ranks[group_mask].sum(axis=0)

    # U statistic
    U = R1 - n1 * (n1 + 1) / 2.0

    # Expected value and tie-corrected variance under H0
    mu = n1 * n2 / 2.0
    # sigma^2 = n1 * n2 / 12 * ((N + 1) * tie_correction)
    # where tie_correction absorbs the tie term
    variance = (n1 * n2 / 12.0) * (n_total + 1) * tie_correction

    # Avoid division by zero for constant columns
    sigma = np.sqrt(np.maximum(variance, 1e-300))
    z = (U - mu) / sigma

    # Two-sided p-value from normal approximation
    p = 2.0 * norm.sf(np.abs(z))

    return U, z, p


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def rank_genes_groups(
    expression: np.ndarray,
    group_labels: np.ndarray,
    gene_ids: np.ndarray,
    gene_names: Optional[np.ndarray] = None,
    params: Optional[WilcoxonParams] = None,
) -> pd.DataFrame:
    """Wilcoxon rank-sum marker gene detection (one-vs-rest).

    Parameters
    ----------
    expression : np.ndarray
        Dense expression matrix ``(n_cells, n_genes)``.
    group_labels : np.ndarray
        Group labels ``(n_cells,)``, one per cell.
    gene_ids : np.ndarray
        Gene identifiers ``(n_genes,)``.
    gene_names : np.ndarray, optional
        Gene names ``(n_genes,)``.  Carried through to output.
    params : WilcoxonParams, optional
        Test parameters.  Defaults to ``WilcoxonParams()``.

    Returns
    -------
    pd.DataFrame
        Results with columns: ``group``, ``gene_id``, ``gene_name`` (if
        provided), ``u_statistic``, ``z_score``, ``pvalue``,
        ``pvalue_adj``, ``fold_change``, ``fraction_expressed``,
        ``fraction_expressed_rest``, ``log2_fold_change``.
        Rows are sorted by ``(group, pvalue_adj, -fold_change)``.
    """
    if params is None:
        params = WilcoxonParams()

    n_cells, n_genes = expression.shape
    groups = np.unique(group_labels)

    logger.info(
        "Wilcoxon rank-sum: %d cells, %d genes, %d groups",
        n_cells, n_genes, len(groups),
    )

    # --- Compute ranks and tie correction once ---
    ranks = _compute_rank_matrix(expression)
    tc = (
        _compute_tie_correction(expression, n_cells)
        if params.tie_correction
        else np.ones(n_genes, dtype=np.float64)
    )

    # --- Run per-group tests ---
    all_results = []

    for group in groups:
        mask = group_labels == group
        n_in = int(mask.sum())
        n_out = n_cells - n_in

        if n_in < 2 or n_out < 2:
            logger.warning(
                "Group '%s' has %d cells (rest=%d) — skipping.",
                group, n_in, n_out,
            )
            continue

        # Per-gene fold change and fraction expressed (pre-filter)
        mean_in = expression[mask].mean(axis=0)
        mean_out = expression[~mask].mean(axis=0)
        frac_in = (expression[mask] > 0).mean(axis=0)
        frac_out = (expression[~mask] > 0).mean(axis=0)
        denom = np.where(mean_out > 0, mean_out, 1e-10)
        fc = mean_in / denom

        # Pre-filter candidates
        candidate_mask = (fc >= params.min_fold_change) & (
            frac_in >= params.min_fraction_expressed
        )
        candidate_idx = np.where(candidate_mask)[0]

        if len(candidate_idx) == 0:
            logger.info(
                "Group '%s': no candidates passed pre-filter.", group
            )
            continue

        # Limit candidates by fold change rank
        if len(candidate_idx) > params.max_candidates:
            top_order = np.argsort(-fc[candidate_idx])[: params.max_candidates]
            candidate_idx = candidate_idx[top_order]

        # Subset ranks and tie-correction for candidates
        ranks_sub = ranks[:, candidate_idx]
        tc_sub = tc[candidate_idx]

        U, z, p = _wilcoxon_one_vs_rest(ranks_sub, mask, n_cells, tc_sub)

        # Multiple testing correction
        p_adj = np.array(
            apply_correction(p.tolist(), method=params.correction_method)
        )

        # Build results for this group
        n_cand = len(candidate_idx)
        gids = gene_ids[candidate_idx]
        log2fc = np.log2(np.maximum(fc[candidate_idx], 1e-300))

        group_df = pd.DataFrame(
            {
                "group": [str(group)] * n_cand,
                "gene_id": gids,
                "u_statistic": U,
                "z_score": z,
                "pvalue": p,
                "pvalue_adj": p_adj,
                "fold_change": fc[candidate_idx],
                "log2_fold_change": log2fc,
                "fraction_expressed": frac_in[candidate_idx],
                "fraction_expressed_rest": frac_out[candidate_idx],
            }
        )

        if gene_names is not None:
            group_df.insert(2, "gene_name", gene_names[candidate_idx])

        all_results.append(group_df)

    if not all_results:
        cols = [
            "group", "gene_id", "u_statistic", "z_score", "pvalue",
            "pvalue_adj", "fold_change", "log2_fold_change",
            "fraction_expressed", "fraction_expressed_rest",
        ]
        if gene_names is not None:
            cols.insert(2, "gene_name")
        return pd.DataFrame(columns=cols)

    results = pd.concat(all_results, ignore_index=True)

    # Sort: group, then by significance (ascending adj-p, descending fc)
    results = results.sort_values(
        ["group", "pvalue_adj", "fold_change"],
        ascending=[True, True, False],
    ).reset_index(drop=True)

    logger.info(
        "Wilcoxon results: %d total rows, %d significant (alpha=%.3f)",
        len(results),
        int((results["pvalue_adj"] <= params.alpha).sum()),
        params.alpha,
    )

    return results


# ---------------------------------------------------------------------------
# Conversion to GeneSetCollection
# ---------------------------------------------------------------------------


def results_to_gene_set_collection(
    results_df: pd.DataFrame,
    background_genes: Set[str],
    top_n: int = 200,
    alpha: float = 0.05,
    gene_col: str = "gene_name",
) -> GeneSetCollection:
    """Convert Wilcoxon results DataFrame to a GeneSetCollection.

    For each group, takes the top *top_n* significant genes (adj-p ≤ alpha),
    sorted by fold-change descending.

    Parameters
    ----------
    results_df : pd.DataFrame
        Output of :func:`rank_genes_groups`.
    background_genes : Set[str]
        Background gene universe.
    top_n : int
        Maximum markers per group.
    alpha : float
        Adjusted p-value threshold.
    gene_col : str
        Column to use for gene identifiers in the GeneSet.  Defaults to
        ``"gene_name"``; falls back to ``"gene_id"`` if absent.

    Returns
    -------
    GeneSetCollection
    """
    if gene_col not in results_df.columns:
        gene_col = "gene_id"

    gene_sets: Dict[str, GeneSet] = {}

    for group, gdf in results_df.groupby("group"):
        sig = gdf[gdf["pvalue_adj"] <= alpha].copy()
        sig = sig.sort_values("fold_change", ascending=False).head(top_n)

        if sig.empty:
            logger.info("Group '%s': no significant markers at alpha=%.3f", group, alpha)
            continue

        genes = set(sig[gene_col].tolist())
        scores = {
            row[gene_col]: round(row["fold_change"], 4)
            for _, row in sig.iterrows()
        }
        pvals = {
            row[gene_col]: float(row["pvalue_adj"])
            for _, row in sig.iterrows()
        }

        gene_sets[str(group)] = GeneSet(
            name=str(group),
            genes=genes,
            source=f"wilcoxon:{group}",
            metadata={
                "fold_changes": scores,
                "adjusted_pvalues": pvals,
                "method": "wilcoxon",
            },
        )

    logger.info(
        "Built GeneSetCollection: %d groups with significant markers",
        len(gene_sets),
    )

    return GeneSetCollection(
        gene_sets=gene_sets,
        background_genes=background_genes,
        source_description=(
            f"Wilcoxon rank-sum markers (top_n={top_n}, alpha={alpha})"
        ),
        metadata={"method": "wilcoxon", "top_n": top_n, "alpha": alpha},
    )
