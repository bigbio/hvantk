"""The pandas path onto a cohort manifest -- ``attach()``'s non-Hail sibling.

``rerank`` (``hvantk/algorithms/rerank/``) is pandas + scikit-learn with zero Hail, so
a :class:`~hvantk.algorithms.cohort.spec.CohortManifest` cannot reach it through
:func:`hvantk.algorithms.cohort.attach.attach`, which joins onto a Hail Layer-1 table.
This module is the loader that does reach it: it reads the same on-disk cohort table
``attach`` and ``hvantk cohort validate`` read, and must never disagree with them about
delimiter, gzip handling, or which rows/columns are valid -- so it reuses
:mod:`hvantk.algorithms.cohort.checks` rather than re-deriving any of that logic.

Layering: stdlib + pandas + sibling ``hvantk.algorithms.cohort`` modules only. This
module must never import Hail, ``hvantk.skills``, ``hvantk.tools``, or ``click``; every
function here raises plain ``ValueError``.
"""
from __future__ import annotations

import pandas as pd

from hvantk.algorithms.cohort.checks import (
    check_declared_columns_exist,
    check_key_column_exists,
    detect_delimiter,
    is_compressed,
    read_header,
)


def _read_table(manifest) -> pd.DataFrame:
    """Read the cohort's on-disk table with the delimiter/gzip rules ``checks``
    already applies, so validation and loading can never disagree about either."""
    delimiter = detect_delimiter(manifest.table)
    compression = "gzip" if is_compressed(manifest.table) else None
    return pd.read_csv(manifest.table, sep=delimiter, compression=compression)


def _check_no_duplicate_genes(genes: pd.Series, manifest) -> None:
    """Fail loud on more than one row for the same gene -- the same rule
    :func:`hvantk.algorithms.cohort.attach.attach` enforces on its Hail side.

    A cohort table is a spine-mappable gene key plus a cohort-derived prior
    statistic: more than one row per gene means a downstream consumer would have to
    pick one arbitrarily and silently drop the rest.
    """
    counts = genes.value_counts()
    duplicated = sorted(str(g) for g in counts[counts > 1].index)
    if not duplicated:
        return
    n_dup_rows = int((counts[counts > 1] - 1).sum())
    raise ValueError(
        f"cohort {manifest.name!r} table has {n_dup_rows} duplicate row(s) for gene "
        f"key(s) {', '.join(duplicated[:10])}; the cohort table must carry exactly "
        "one row per gene -- aggregate the table before loading"
    )


def load_cohort_frame(manifest) -> pd.DataFrame:
    """Load a cohort manifest's table as a pandas frame: one row per gene.

    Returns a ``gene`` column (renamed from ``manifest.key_column``) plus every column
    in ``manifest.declared_columns()`` -- the prior column and every cohort axis
    column, in declaration order. Enforces the same contract
    ``hvantk cohort validate``/``attach`` enforce: the key column and every declared
    column must exist in the table's header, and the table must carry exactly one row
    per gene.

    Raises
    ------
    ValueError
        If the table is missing or empty (:func:`~hvantk.algorithms.cohort.checks.read_header`),
        if the key column or a declared column is missing from the header, or if the
        table has more than one row for the same gene.
    """
    header = read_header(manifest.table)
    check_key_column_exists(manifest, header)
    check_declared_columns_exist(manifest, header)

    df = _read_table(manifest)
    _check_no_duplicate_genes(df[manifest.key_column], manifest)

    declared = [c for c in manifest.declared_columns() if c != manifest.key_column]
    frame = df[[manifest.key_column] + declared].rename(
        columns={manifest.key_column: "gene"}
    )
    return frame.reset_index(drop=True)


def load_prior_frame(manifest) -> pd.DataFrame:
    """Load a cohort manifest's prior statistic in ``PriorSpec.load()``'s shape.

    Returns exactly two columns, ``unit`` (the gene key) and ``prior_stat`` (the
    manifest's declared prior column) -- the same column names
    ``hvantk.algorithms.rerank.config.PriorSpec.load()`` returns today, so a
    ``CohortManifest`` can stand in for a ``PriorSpec`` without anything downstream of
    this call noticing the difference.
    """
    frame = load_cohort_frame(manifest)
    return frame[["gene", manifest.prior.column]].rename(
        columns={"gene": "unit", manifest.prior.column: "prior_stat"}
    )
