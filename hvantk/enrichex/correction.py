"""
Backward-compatible re-export of correction utilities.

The canonical location is :mod:`hvantk.algorithms.statistics.correction`.
"""

from hvantk.algorithms.statistics.correction import apply_correction, fdr_threshold  # noqa: F401

__all__ = ["apply_correction", "fdr_threshold"]
