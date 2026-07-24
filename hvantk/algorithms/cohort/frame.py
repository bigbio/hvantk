"""The pandas path onto a cohort manifest -- ``attach()``'s non-Hail sibling.

``rerank`` (``hvantk/algorithms/rerank/``) is pandas + scikit-learn with zero Hail, so
a :class:`~hvantk.algorithms.cohort.spec.CohortManifest` cannot reach it through
:func:`hvantk.algorithms.cohort.attach.attach`, which joins onto a Hail Layer-1 table.
This module is the loader that does reach it: it reads the same on-disk cohort table
``attach`` and ``hvantk cohort validate`` read, and must never disagree with them about
delimiter, gzip handling, or which *columns are present in the header* -- so it reuses
:mod:`hvantk.algorithms.cohort.checks` rather than re-deriving any of that logic.

``validate`` is header-only by design: it never opens a single data row. Everything
about the data rows themselves -- one row per gene, a non-null gene key on every row,
whitespace-padded cell values that would otherwise silently miss a join -- is enforced
only here, at load time. A table can pass ``hvantk cohort validate`` cleanly and still
raise from :func:`load_cohort_frame`; that is by design, not a bug to "fix" by teaching
``validate`` to read rows.

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


def _check_no_duplicate_column_labels(columns: pd.Index, manifest) -> None:
    """Fail loud if stripping whitespace from column names collapsed two distinct
    headers onto the same label.

    A header like ``"gene\\tgene \\tminp"`` has two distinct raw names, but
    :func:`_read_table` strips both to ``"gene"`` (matching
    :func:`~hvantk.algorithms.cohort.checks.read_header`, which strips too) so the
    literal ``"gene, minp"``-with-a-space style validates and loads the same way.
    That stripping can also manufacture a duplicate that was not there before: without
    this check, ``df[manifest.key_column]`` would silently return a two-column
    ``DataFrame`` instead of a ``Series``, and the first thing to notice would be an
    ``AttributeError`` several calls downstream that names neither the manifest nor
    the column.
    """
    dupes = sorted(set(columns[columns.duplicated()]))
    if not dupes:
        return
    names = ", ".join(repr(d) for d in dupes)
    raise ValueError(
        f"cohort {manifest.name!r} table header has column name(s) {names} that "
        "collapse onto a duplicate label once whitespace is stripped from column "
        "names; rename the source column(s) in the table so they remain distinct "
        "after stripping"
    )


def _read_table(manifest, *, raw: bool = False) -> pd.DataFrame:
    """Read the cohort's on-disk table with the delimiter/gzip rules ``checks``
    already applies, so validation and loading can never disagree about either.

    Column *names* are stripped of surrounding whitespace, matching
    :func:`~hvantk.algorithms.cohort.checks.read_header` (which strips too): a header
    like ``"gene, minp"`` -- a common comma-delimited style -- has header names that
    validate and load agree on.

    That alone is not enough: a ``", "``-delimited table also pads every *value*
    after the first field with a leading space (``"g0"`` arrives as ``" g0"``). If
    only names were stripped, such a table would load "successfully" with
    space-padded keys that then silently fail to join anywhere downstream --
    strictly worse than the pre-stripping behaviour of failing loud with a
    ``KeyError``. So column *values* are stripped too, for every object-dtype
    (string) column, not only the key column: any declared column could equally be
    used as a join or grouping key by a downstream consumer, and there is no
    legitimate cohort-table use for leading/trailing whitespace in a cell.

    ``raw=True`` skips both the NA-token coercion pandas normally applies and the
    value-stripping above, returning literal, unparsed text for every cell (used by
    :func:`_raw_key_value_at` to show a user the exact raw cell behind a rejected
    row -- showing the *stripped* or *NaN-coerced* value there would defeat the
    point of naming a "raw" value in the first place).
    """
    delimiter = detect_delimiter(manifest.table)
    compression = "gzip" if is_compressed(manifest.table) else None
    read_kwargs = {"dtype": str, "keep_default_na": False} if raw else {}
    df = pd.read_csv(
        manifest.table, sep=delimiter, compression=compression, **read_kwargs
    )
    df.columns = df.columns.str.strip()
    _check_no_duplicate_column_labels(df.columns, manifest)
    if not raw:
        for col in df.select_dtypes(include="object").columns:
            df[col] = df[col].str.strip()
    return df


def _raw_key_value_at(manifest, row_index: int) -> str:
    """The literal, unparsed text of the key column at ``row_index`` (0-based, data
    rows only, header excluded), bypassing both pandas' default NA-token coercion
    and the whitespace-stripping :func:`_read_table` normally applies.

    Used only to build a helpful error message: showing the user ``NaN`` (what the
    cell became) instead of ``"NA"`` (what was actually typed in the file) makes the
    row much harder to find with a text search.
    """
    raw = _read_table(manifest, raw=True)
    return raw[manifest.key_column].iloc[row_index]


def _check_no_null_genes(genes: pd.Series, manifest) -> None:
    """Fail loud on an empty/blank/null gene key instead of letting it enter the frame
    as ``NaN``.

    ``pd.Series.value_counts()`` (used by :func:`_check_no_duplicate_genes`) drops NaN
    by default, so a table with several blank key cells would otherwise sail straight
    through the duplicate-row check, enter the frame as ``NaN`` genes, and then vanish
    silently in a downstream left-merge -- neither rejected nor reported, in violation
    of the "exactly one row per gene" contract.

    A cell counts as blank here either because it is genuinely empty/whitespace-only,
    or because ``pandas.read_csv`` parsed it as one of its default NA tokens (``"NA"``,
    ``"N/A"``, ``"NULL"``, ``"None"``, ``"nan"``, ``"null"``, ...) -- a gene literally
    named one of those strings is indistinguishable, post-parse, from a truly missing
    cell, so it is rejected the same way. (``"0"`` is not among pandas' NA tokens and
    is accepted normally.) The message below names the raw text so a cell that
    *looks* non-empty is not reported as "empty".
    """
    blank = genes.isna() | (genes.astype(str).str.strip() == "")
    n_blank = int(blank.sum())
    if not n_blank:
        return
    first_row = int(blank[blank].index[0])
    raw_value = _raw_key_value_at(manifest, first_row)
    raise ValueError(
        f"cohort {manifest.name!r} table has {n_blank} row(s) whose "
        f"{manifest.key_column!r} (key_column) value is blank, or was parsed by "
        "pandas as an NA token (pandas.read_csv's default NA strings include '', "
        "'NA', 'N/A', 'NULL', 'None', 'nan', 'null', among others -- '0' is not one "
        f"of them and is accepted); e.g. data row {first_row + 1} has raw "
        f"{manifest.key_column!r} value {raw_value!r}. The cohort table must carry "
        "a valid gene key on every row -- fix or drop these rows before loading"
    )


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


def load_cohort_frame(manifest, *, include_prior: bool = True) -> pd.DataFrame:
    """Load a cohort manifest's table as a pandas frame: one row per gene.

    Returns a ``gene`` column (renamed from ``manifest.key_column``) plus every column
    in ``manifest.declared_columns()`` -- the prior column and every cohort axis
    column, in declaration order (or just the axis columns; see ``include_prior``
    below).

    Column-presence enforcement matches ``hvantk cohort validate``/``attach``
    exactly: the key column and every declared column must exist in the table's
    header (both sides agree on whitespace-stripped names). Everything else this
    function enforces -- one row per gene, a non-null gene key on every row, no
    post-strip duplicate column labels, whitespace-padded values that would silently
    miss a join -- is checked here ONLY. ``validate`` is header-only by design and
    never reads a single data row, so a table can pass ``hvantk cohort validate``
    cleanly and still raise here.

    ``include_prior`` (default ``True``) controls whether the returned frame carries
    the manifest's prior column (``manifest.prior.column``) alongside the axis
    columns: ``True`` returns ``manifest.declared_columns()``, ``False`` returns only
    ``manifest.axis_columns()``. Pass ``False`` when the caller has already consumed
    the prior under its own name (e.g. as ``prior_stat``) and only wants the cohort
    axis columns -- most notably ``engine.rerank()``'s audit merge, where
    re-including the raw prior column would (a) add nothing not already carried as
    ``prior_stat`` and (b) make the prior column a spurious collision candidate
    against a feature axis that legitimately reuses the same column name (e.g. a
    "burden" axis whose model feature is the same p-value the cohort declares as its
    prior).

    Raises
    ------
    ValueError
        If the table is missing or empty (:func:`~hvantk.algorithms.cohort.checks.read_header`),
        if the key column or a declared column is missing from the header, if
        stripping whitespace from the header collapses two column names onto the
        same label, if the table has more than one row for the same gene, or if the
        gene key is blank or a pandas NA token on any row.
    """
    header = read_header(manifest.table)
    check_key_column_exists(manifest, header)
    check_declared_columns_exist(manifest, header)

    df = _read_table(manifest)
    _check_no_null_genes(df[manifest.key_column], manifest)
    _check_no_duplicate_genes(df[manifest.key_column], manifest)

    cols = manifest.declared_columns() if include_prior else manifest.axis_columns()
    declared = [c for c in cols if c != manifest.key_column]
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
