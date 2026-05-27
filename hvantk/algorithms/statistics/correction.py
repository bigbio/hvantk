"""
Multiple testing correction methods.

Provides Bonferroni and Benjamini-Hochberg (FDR) correction methods.
"""

import logging
from typing import List, Literal

import numpy as np

logger = logging.getLogger(__name__)


def apply_correction(
    p_values: List[float],
    method: Literal["bonferroni", "benjamini-hochberg", "none"] = "benjamini-hochberg",
    n_total: int | None = None,
) -> List[float]:
    """Apply multiple testing correction to p-values.

    Parameters
    ----------
    p_values : List[float]
        Raw p-values to correct.
    method : str
        Correction method:
        - "bonferroni": Bonferroni correction (conservative)
        - "benjamini-hochberg": Benjamini-Hochberg FDR correction
        - "none": No correction (return original p-values)
    n_total : int, optional
        Total number of hypotheses tested.  When provided, the correction
        uses this value instead of ``len(p_values)``.  This is useful when
        a pre-filter reduces the set of p-values but the correction should
        account for the original number of tests (Seurat-style correction;
        see ``p.adjust(p, method, n)`` in R).

    Returns
    -------
    List[float]
        Adjusted p-values

    Examples
    --------
    >>> p_values = [0.001, 0.01, 0.05, 0.1]
    >>> apply_correction(p_values, method="bonferroni")
    [0.004, 0.04, 0.2, 0.4]

    >>> apply_correction(p_values, method="benjamini-hochberg")
    [0.004, 0.02, 0.0667, 0.1]
    """
    if not p_values:
        return []

    n_tests = n_total if n_total is not None else len(p_values)
    if n_tests < len(p_values):
        raise ValueError(
            f"n_total ({n_tests}) must be >= len(p_values) ({len(p_values)})"
        )

    if method == "none":
        logger.debug("No correction applied")
        return list(p_values)

    elif method == "bonferroni":
        logger.debug(f"Applying Bonferroni correction (n={n_tests})")
        adjusted = [min(1.0, p * n_tests) for p in p_values]
        return adjusted

    elif method == "benjamini-hochberg":
        logger.debug(f"Applying Benjamini-Hochberg correction (n={n_tests})")

        # Convert to numpy for efficient computation
        p_array = np.array(p_values)
        n_p = len(p_values)

        # Sort indices by p-value
        sorted_indices = np.argsort(p_array)
        sorted_p = p_array[sorted_indices]

        # BH adjustment: p_adj[i] = p[i] * n_total / rank[i]
        # Iterate over the actual p-values (len(p_values)) but use
        # n_tests (which may be n_total > len(p_values)) in the formula.
        # This matches R's p.adjust(p, method="BH", n=n_total).
        adjusted = np.zeros(n_p)

        # Compute adjusted p-values in reverse order
        for i in range(n_p - 1, -1, -1):
            rank = i + 1
            adjusted[sorted_indices[i]] = sorted_p[i] * n_tests / rank

        # Ensure monotonicity (cumulative minimum from end)
        # This ensures that if p[i] < p[j], then adj_p[i] <= adj_p[j]
        for i in range(n_p - 2, -1, -1):
            if adjusted[sorted_indices[i]] > adjusted[sorted_indices[i + 1]]:
                adjusted[sorted_indices[i]] = adjusted[sorted_indices[i + 1]]

        # Cap at 1.0
        adjusted = np.minimum(adjusted, 1.0)

        return adjusted.tolist()

    else:
        raise ValueError(
            f"Unknown correction method: {method}. "
            f"Must be one of: bonferroni, benjamini-hochberg, none"
        )


def fdr_threshold(
    p_values: List[float],
    alpha: float = 0.05,
) -> float:
    """Calculate FDR threshold using Benjamini-Hochberg procedure.

    This returns the largest p-value that would be considered significant
    at the given FDR level.

    Parameters
    ----------
    p_values : List[float]
        Raw p-values
    alpha : float
        Desired FDR level (default: 0.05)

    Returns
    -------
    float
        FDR threshold p-value

    Examples
    --------
    >>> p_values = [0.001, 0.01, 0.05, 0.1, 0.2]
    >>> fdr_threshold(p_values, alpha=0.05)
    0.05
    """
    if not p_values:
        return 0.0

    n_tests = len(p_values)
    sorted_p = sorted(p_values)

    # Find largest i where p[i] <= (i/n) * alpha
    threshold = 0.0
    for i, p in enumerate(sorted_p, start=1):
        if p <= (i / n_tests) * alpha:
            threshold = p

    return threshold
