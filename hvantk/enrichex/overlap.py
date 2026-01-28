"""
Overlap enrichment analysis using Fisher's exact test.

This module provides functions for testing whether a query gene list is
enriched in gene sets using Fisher's exact test via Hail.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

try:  # Optional dependency: Hail is only needed for actual enrichment tests
    import hail as hl
except ModuleNotFoundError as exc:  # pragma: no cover - depends on env
    hl = None  # type: ignore
    _HAIL_IMPORT_ERROR = exc
else:
    _HAIL_IMPORT_ERROR = None

from hvantk.enrichex.correction import apply_correction
from hvantk.enrichex.gene_sets import GeneSetCollection

logger = logging.getLogger(__name__)


@dataclass
class OverlapResult:
    """Result of Fisher's exact test for one gene set.

    Attributes
    ----------
    gene_set_name : str
        Name of the tested gene set
    n_query : int
        Number of query genes in background universe
    n_gene_set : int
        Size of the gene set (in background)
    n_overlap : int
        Number of genes in overlap
    n_background : int
        Size of background universe
    p_value : float
        Raw p-value from Fisher's exact test
    odds_ratio : float
        Odds ratio (enrichment factor)
    ci_lower : float
        Lower bound of 95% confidence interval
    ci_upper : float
        Upper bound of 95% confidence interval
    overlap_genes : List[str]
        List of genes in the overlap
    p_adjusted : Optional[float]
        Adjusted p-value (after multiple testing correction)
    """

    gene_set_name: str
    n_query: int
    n_gene_set: int
    n_overlap: int
    n_background: int
    p_value: float
    odds_ratio: float
    ci_lower: float
    ci_upper: float
    overlap_genes: List[str]
    p_adjusted: Optional[float] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization.

        Returns
        -------
        Dict[str, Any]
            Dictionary representation
        """
        return {
            "gene_set_name": self.gene_set_name,
            "n_query": self.n_query,
            "n_gene_set": self.n_gene_set,
            "n_overlap": self.n_overlap,
            "n_background": self.n_background,
            "p_value": self.p_value,
            "odds_ratio": self.odds_ratio,
            "ci_lower": self.ci_lower,
            "ci_upper": self.ci_upper,
            "overlap_genes": self.overlap_genes,
            "p_adjusted": self.p_adjusted,
        }


def compute_overlap_enrichment(
    query_genes: List[str],
    gene_set_collection: GeneSetCollection,
    correction_method: str = "benjamini-hochberg",
) -> List[OverlapResult]:
    """Test enrichment of query genes in each gene set using Fisher's exact test.

    This function uses Hail's fisher_exact_test() to compute enrichment
    statistics for each gene set in the collection.

    The test uses a 2x2 contingency table:
                      | In Gene Set | Not in Gene Set |
        In Query      |     a       |       b         |
        Not in Query  |     c       |       d         |

    where:
    - a = number of query genes in gene set
    - b = number of query genes not in gene set
    - c = number of non-query genes in gene set
    - d = number of non-query genes not in gene set

    Parameters
    ----------
    query_genes : List[str]
        User's gene list to test for enrichment
    gene_set_collection : GeneSetCollection
        Collection of gene sets to test against
    correction_method : str
        Multiple testing correction method: "bonferroni", "benjamini-hochberg", or "none"

    Returns
    -------
    List[OverlapResult]
        Results for each gene set, sorted by p-value (ascending)

    Examples
    --------
    >>> from hvantk.enrichex import GeneSetCollection, compute_overlap_enrichment
    >>> query = ["BRCA1", "TP53", "EGFR"]
    >>> gene_sets = GeneSetCollection.load("cell_type_markers.json")
    >>> results = compute_overlap_enrichment(query, gene_sets)
    >>> for r in results[:5]:
    ...     print(f"{r.gene_set_name}: p={r.p_adjusted:.2e}, OR={r.odds_ratio:.2f}")
    """
    _require_hail()
    query_set = set(query_genes)
    background = gene_set_collection.background_genes

    # Intersect query with background
    query_in_background = query_set & background
    n_query = len(query_in_background)
    n_background = len(background)

    logger.info(f"Query: {len(query_set)} genes, {n_query} in background")
    logger.info(f"Background: {n_background} genes")
    logger.info(f"Testing {len(gene_set_collection)} gene sets")

    if n_query == 0:
        logger.warning("No query genes found in background universe!")
        return []

    results = []

    for gene_set in gene_set_collection:
        # Get gene set genes that are in background
        gs_genes = gene_set.genes & background
        n_gene_set = len(gs_genes)

        # Skip empty gene sets
        if n_gene_set == 0:
            logger.warning(
                f"Gene set '{gene_set.name}' has no genes in background, skipping"
            )
            continue

        # Compute overlap
        overlap = query_in_background & gs_genes
        n_overlap = len(overlap)

        # Build 2x2 contingency table
        # |                  | In Gene Set | Not in Gene Set |
        # | In Query         |     a       |       b         |
        # | Not in Query     |     c       |       d         |
        a = n_overlap
        b = n_query - n_overlap
        c = n_gene_set - n_overlap
        d = n_background - n_query - n_gene_set + n_overlap

        # Sanity check
        if a < 0 or b < 0 or c < 0 or d < 0:
            logger.error(
                f"Invalid contingency table for {gene_set.name}: "
                f"a={a}, b={b}, c={c}, d={d}"
            )
            continue

        # Fisher's exact test using Hail
        try:
            fisher_result = hl.eval(hl.fisher_exact_test(a, b, c, d))

            results.append(
                OverlapResult(
                    gene_set_name=gene_set.name,
                    n_query=n_query,
                    n_gene_set=n_gene_set,
                    n_overlap=n_overlap,
                    n_background=n_background,
                    p_value=fisher_result.p_value,
                    odds_ratio=fisher_result.odds_ratio,
                    ci_lower=fisher_result.ci_95_lower,
                    ci_upper=fisher_result.ci_95_upper,
                    overlap_genes=sorted(overlap),
                )
            )
        except Exception as e:
            logger.error(f"Fisher's exact test failed for {gene_set.name}: {e}")
            continue

    if not results:
        logger.warning("No valid results generated")
        return []

    # Sort by p-value
    results.sort(key=lambda r: r.p_value)

    # Apply multiple testing correction
    logger.info(f"Applying {correction_method} correction to {len(results)} results")
    p_values = [r.p_value for r in results]
    p_adjusted = apply_correction(p_values, method=correction_method)

    for i, r in enumerate(results):
        r.p_adjusted = p_adjusted[i]

    return results


def compute_overlap_enrichment_pandas(
    query_genes: List[str],
    gene_set_collection: GeneSetCollection,
    correction_method: str = "benjamini-hochberg",
):
    """Compute enrichment and return as pandas DataFrame.

    This is a convenience wrapper around compute_overlap_enrichment()
    that returns results as a pandas DataFrame for easy analysis.

    Parameters
    ----------
    query_genes : List[str]
        Query gene list
    gene_set_collection : GeneSetCollection
        Gene set collection
    correction_method : str
        Multiple testing correction method

    Returns
    -------
    pd.DataFrame
        Results DataFrame with columns:
        - gene_set_name
        - n_query, n_gene_set, n_overlap, n_background
        - p_value, odds_ratio, ci_lower, ci_upper
        - p_adjusted
        - overlap_genes (comma-separated string)
        - significant (bool, p_adjusted < 0.05)

    Examples
    --------
    >>> df = compute_overlap_enrichment_pandas(query_genes, gene_sets)
    >>> df[df['significant']].sort_values('p_adjusted')
    """
    _require_hail()
    import pandas as pd

    results = compute_overlap_enrichment(
        query_genes, gene_set_collection, correction_method
    )

    if not results:
        # Return empty DataFrame with expected columns
        return pd.DataFrame(
            columns=[
                "gene_set_name",
                "n_query",
                "n_gene_set",
                "n_overlap",
                "n_background",
                "p_value",
                "odds_ratio",
                "ci_lower",
                "ci_upper",
                "p_adjusted",
                "overlap_genes",
                "significant",
            ]
        )

    df = pd.DataFrame([r.to_dict() for r in results])
    df["overlap_genes"] = df["overlap_genes"].apply(lambda x: ",".join(x))
    df["significant"] = df["p_adjusted"] < 0.05

    return df


def _require_hail() -> None:
    if hl is None:  # pragma: no cover - depends on env
        raise ImportError(
            "Hail is required for overlap enrichment. "
            "Install hvantk with the 'hail' extra."
        ) from _HAIL_IMPORT_ERROR
