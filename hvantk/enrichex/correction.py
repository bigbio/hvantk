"""
Backward-compatible re-export of correction utilities.

The canonical location is :mod:`hvantk.utils.correction`.
"""

from hvantk.utils.correction import apply_correction, fdr_threshold  # noqa: F401

__all__ = ["apply_correction", "fdr_threshold"]
