"""Reduce an expression summary (groups x genes) to a symbol-keyed per-gene feature table.

Layering: pandas/numpy/anndata + stdlib only, plus sibling ``hvantk.algorithms`` modules; never
``hvantk.skills``/``hvantk.tools``. The atlas arrives as data (an AnnData handed in by the caller).
The reduced table then rides the existing P2c-2 ``key: symbol`` path onto the spine.

Specificity is the EWCE fraction ``mean_group / sum_groups(mean)`` (the metric the manuscript used;
the only one that reproduces the locked expression numbers).
"""
from __future__ import annotations

import logging
import re

logger = logging.getLogger(__name__)


def _san(label: str) -> str:
    """Column-safe token from a group label: lowercase, non-alnum -> single underscore."""
    return re.sub(
        r"_+", "_", re.sub(r"[^0-9a-zA-Z]+", "_", label.strip().lower())
    ).strip("_")


def ewce_specificity(mean_gg):
    """Row-normalized mean expression: spec[g, k] = mean[g, k] / sum_k mean[g, :].

    ``mean_gg`` is a genes x groups DataFrame. A gene with zero total mean (expressed nowhere)
    maps to 0 in every group (not NaN), so it reads as non-specific rather than missing.
    """
    row_sums = mean_gg.sum(axis=1)
    return mean_gg.div(row_sums.where(row_sums > 0), axis=0).fillna(0.0)


def _layer_genes_by_group(adata, name, keep_idx, groups):
    """Extract a layer as a genes x groups DataFrame (the summary stores groups x genes)."""
    import numpy as np
    import pandas as pd

    m = adata.layers[name]
    m = m.toarray() if hasattr(m, "toarray") else np.asarray(m)
    return pd.DataFrame(
        m[keep_idx, :].T,
        index=list(adata.var.index),
        columns=[groups[i] for i in keep_idx],
    )


def reduce_matrix_to_gene(adata, mspec):
    """Reduce a summary AnnData to a symbol-keyed per-gene feature table.

    Returns a DataFrame with a ``symbol`` column and the feature columns, one row per gene.
    """
    import pandas as pd

    groups = [str(g).strip() for g in adata.obs.index]
    drop = {str(g).strip() for g in mspec.drop_groups}
    keep_idx = [i for i, g in enumerate(groups) if g not in drop]

    mean_gg = _layer_genes_by_group(adata, "mean", keep_idx, groups)
    cols = {}

    for stat in mspec.stats:
        layer = "fraction_expressed" if stat in ("fraction_expressed", "frac") else stat
        df = (
            mean_gg
            if layer == "mean"
            else _layer_genes_by_group(adata, layer, keep_idx, groups)
        )
        stat_tag = "frac" if layer == "fraction_expressed" else layer
        for g in df.columns:
            cols[f"{mspec.atlas}_{_san(g)}_{stat_tag}"] = df[g]

    if mspec.specificity is not None:
        spec_gg = ewce_specificity(mean_gg)
        targets = [t.strip() for t in mspec.specificity.targets]
        present = [c for c in spec_gg.columns if c in targets]
        if not present:
            raise ValueError(
                f"specificity targets {targets} not found among groups {list(spec_gg.columns)}"
            )
        tgt = spec_gg[present]
        combined = (
            tgt.max(axis=1) if mspec.specificity.combine == "max" else tgt.mean(axis=1)
        )
        cols[f"{mspec.atlas}_{mspec.specificity.name}"] = combined

    out = pd.DataFrame(cols)
    out.index.name = "symbol"
    out = out.groupby(
        level=0
    ).max()  # one row per gene; duplicate symbols collapse to max
    logger.info(
        "reduce_matrix_to_gene: %d genes, %d feature columns", len(out), out.shape[1]
    )
    return out.reset_index()
