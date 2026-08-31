"""Thin wrapper around :mod:`tspex` tissue-specificity metrics.

Delegates computation to the maintained `tspex <https://pypi.org/project/tspex/>`_
package rather than reimplementing Yanai's τ or related indices inline. A single
entry point returns a ``pd.Series`` of per-gene specificity values given a
gene x group expression matrix.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd

from hvantk.core.models.backends import algorithm, Backend

if TYPE_CHECKING:
    # AnnotationTable referenced in string annotations on
    # compute_specificity_artifact; imported under TYPE_CHECKING so the linter
    # sees the name without forcing a runtime import at module load.
    from hvantk.core.models import AnnotationTable  # noqa: F401

logger = logging.getLogger(__name__)

__all__ = ["compute_specificity", "compute_specificity_artifact"]

_SUPPORTED_METHODS = (
    "tau",
    "tsi",
    "gini",
    "shannon_specificity",
    "roku_specificity",
    "counts",
    "spm",
)

Method = Literal[
    "tau",
    "tsi",
    "gini",
    "shannon_specificity",
    "roku_specificity",
    "counts",
    "spm",
]


@algorithm(name="tissue_specificity", backends=[Backend.PANDAS])
def compute_specificity(
    gene_x_group_df: pd.DataFrame,
    method: Method = "tau",
    log: bool = False,
) -> pd.Series:
    """Compute per-gene tissue specificity via ``tspex``.

    Parameters
    ----------
    gene_x_group_df
        Rows = genes, columns = groups (tissues/cell types), values = aggregate
        expression (e.g. median TPM or log-relative abundance). Must be numeric,
        non-negative, and free of NaN rows.
    method
        Specificity metric. ``"tau"`` (Yanai 2005) is the default. Any method
        supported by :class:`tspex.TissueSpecificity` is accepted.
    log
        If ``True``, let tspex log2-transform the input internally. Pass
        already-log-transformed values with ``log=False`` (common for
        proteomics log-relative matrices).

    Returns
    -------
    pandas.Series
        One value per gene, indexed by input row index.

    Raises
    ------
    ImportError
        If ``tspex`` is not installed.
    ValueError
        If the input contains NaNs, negative values, or an unsupported method.
    """
    if method not in _SUPPORTED_METHODS:
        raise ValueError(
            f"Unsupported method '{method}'. Expected one of {_SUPPORTED_METHODS}."
        )

    try:
        import tspex
    except ImportError as exc:
        raise ImportError(
            "tspex is required for tissue specificity metrics. "
            "Install with 'pip install tspex' or the 'constraint' extra."
        ) from exc

    if gene_x_group_df.empty:
        raise ValueError("gene_x_group_df is empty.")
    if gene_x_group_df.shape[1] < 2:
        raise ValueError(
            "Need at least 2 groups (columns) to compute specificity; got "
            f"{gene_x_group_df.shape[1]}."
        )

    numeric = gene_x_group_df.apply(pd.to_numeric, errors="coerce")
    if numeric.isna().any().any():
        n_bad = int(numeric.isna().any(axis=1).sum())
        raise ValueError(
            f"Input has NaN values in {n_bad} row(s); drop or impute before calling."
        )
    if (numeric.values < 0).any():
        raise ValueError(
            "Input contains negative values; tspex expects non-negative expression."
        )

    logger.info(
        "Computing tspex %s on matrix %d genes x %d groups (log=%s)",
        method,
        numeric.shape[0],
        numeric.shape[1],
        log,
    )

    ts = tspex.TissueSpecificity(numeric, method=method, log=log)
    values = ts.tissue_specificity
    values.name = method

    if values.isna().any():
        logger.warning(
            "tspex returned NaN for %d gene(s) (likely all-zero rows).",
            int(values.isna().sum()),
        )

    return values


@algorithm(
    name="tissue_specificity_artifact",
    backends=[Backend.PANDAS],
    inputs={"ann": "AnnotationTable"},
    outputs={"specificity": "AnnotationTable"},
)
def compute_specificity_artifact(
    ann: "AnnotationTable",
    method: "Method" = "tau",
    log: bool = False,
) -> "AnnotationTable":
    """Phase P artifact-typed wrapper for compute_specificity.

    Accepts an AnnotationTable of gene-by-tissue/cell-type expression
    (rows = genes, columns = groups). Delegates to compute_specificity
    via .to_pandas(); wraps the resulting pd.Series back into an
    AnnotationTable carrying chained provenance from the input.

    The legacy pd.DataFrame-typed compute_specificity stays as the
    canonical implementation; this wrapper is the canonical entry point
    for callers consuming the artifact contract.
    """
    from hvantk.core.models import AnnotationTable

    df = ann.to_pandas()
    # The gene identifier column is conventionally the index; restore it.
    if len(df.columns) > 0 and df.columns[0] in {
        "gene_id", "gene_symbol", "gene", "ensembl_id"
    }:
        df = df.set_index(df.columns[0])

    result = compute_specificity(df, method=method, log=log)
    # result is a pd.Series — promote to a 2-column DataFrame: gene + specificity
    out_df = pd.DataFrame({
        "gene_id": result.index,
        "specificity": result.values,
    })
    return AnnotationTable.from_pandas(out_df, provenance=ann.provenance)


def tau_yanai_reference(
    gene_x_group_df: pd.DataFrame, eps: float = 1e-9
) -> pd.Series:
    """Reference Yanai τ used only for validation against ``compute_specificity``.

    Not part of the public API — exported for the smoke test so we can confirm
    tspex's τ matches the textbook definition on a small fixture.
    """
    x = gene_x_group_df.to_numpy(dtype=float)
    n = x.shape[1]
    if n < 2:
        raise ValueError("Need >= 2 columns.")

    row_max = x.max(axis=1)
    safe_max = np.where(row_max > eps, row_max, np.nan)
    normalised = x / safe_max[:, None]
    tau = (1.0 - normalised).sum(axis=1) / (n - 1)
    return pd.Series(tau, index=gene_x_group_df.index, name="tau")
