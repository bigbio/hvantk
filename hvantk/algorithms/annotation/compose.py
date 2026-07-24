"""Compose prepared Layer-1 sources onto the gene spine (Stage 2 of the annotation build).

`compose` is Stage 2 of the annotation build (design §3): it left-joins every axis's
prepared table (Stage 1, ``prepare.py``) onto the ``gene_id``-keyed spine, driven
entirely by a :class:`~hvantk.algorithms.annotation.spec.FeatureSpec`. It is deliberately
generic -- axis names, output columns, and any cohort-specific assumptions come only from
the spec and the ``prepared_by_axis`` mapping the caller hands in; nothing about a
specific axis or cohort is hard-coded here.

A gene absent from an axis's prepared table stays missing (Hail missing) in that axis's
declared columns -- compose never fills in a default. Each axis additionally gets a
``{axis}_present`` boolean column, true iff the gene had a row in that axis's prepared
table, so downstream consumers can tell "missing because absent from source" apart from
any other reason a value might be null.

Layering: this module imports Hail and sibling ``hvantk.algorithms.annotation`` modules
only -- it must never import ``hvantk.skills`` or ``hvantk.tools`` (the source is data,
already reduced to a Hail Table by Stage 1; compose has no business reading raw sources).
"""
from __future__ import annotations

import logging

logger = logging.getLogger(__name__)


def compose(spine, prepared_by_axis, spec):
    """Left-join every ``spec.layer1`` axis's prepared table onto the spine.

    Parameters
    ----------
    spine : hail.Table
        The ``gene_id``-keyed spine (Stage 0).
    prepared_by_axis : dict[str, hail.Table]
        Axis label -> that axis's prepared, ``gene_id``-keyed table (Stage 1 output,
        e.g. :func:`hvantk.algorithms.annotation.prepare.prepare_source`), one entry
        per axis in ``spec.layer1``.
    spec : hvantk.algorithms.annotation.spec.FeatureSpec

    Returns
    -------
    (hail.Table, dict)
        The composed, ``gene_id``-keyed table -- spine columns, every axis's declared
        columns, and a ``{axis}_present`` bool per axis -- and a manifest reporting,
        per axis and column, the non-null rate over the spine and the positive rate
        (fraction > 0) among non-null values.

    Raises
    ------
    ValueError
        If two ``spec.layer1`` axes declare the same output column name. Checked in
        pure Python before any Hail call.
    KeyError
        If ``prepared_by_axis`` is missing an axis that ``spec.layer1`` declares.
    """
    import hail as hl

    _check_no_column_collisions(spec)

    ht = spine
    for entry in spec.layer1:
        axis = entry.axis
        if axis not in prepared_by_axis:
            raise KeyError(
                f"prepared_by_axis is missing axis {axis!r}, declared in spec "
                f"{spec.name!r}"
            )
        prepared = prepared_by_axis[axis]
        # `prepared[ht.gene_id]` is Hail's index-join idiom (Table.__getitem__ ->
        # Table.index): a StructExpression of prepared's non-key fields, missing where
        # ht.gene_id has no match in prepared -- exactly the left join we want, and the
        # struct being missing (rather than any individual field) is what makes the
        # presence flag correct regardless of whether a matched row's own column values
        # happen to be null.
        joined = prepared[ht.gene_id]
        updates = {col: joined[col] for col in entry.columns}
        updates[f"{axis}_present"] = hl.is_defined(joined)
        ht = ht.annotate(**updates)

    manifest = _build_manifest(ht, spec)
    return ht, manifest


def _check_no_column_collisions(spec) -> None:
    """Raise if two ``spec.layer1`` axes declare the same output column name.

    Pure Python, no Hail -- runs before any Hail call so the check (and the tests that
    exercise it) do not need a Hail session.
    """
    owner: dict[str, str] = {}
    for entry in spec.layer1:
        for col in entry.columns:
            if col in owner:
                raise ValueError(
                    f"duplicate output column {col!r}: declared by axis {owner[col]!r} "
                    f"and axis {entry.axis!r}; name columns uniquely across the spec's "
                    "layer1 axes"
                )
            owner[col] = entry.axis


def _build_manifest(ht, spec) -> dict:
    """Per axis+column non-null rate (over the spine) and positive rate (among non-null).

    Computed with a single ``ht.aggregate`` call over every declared column, per the
    design's manifest contract.
    """
    import hail as hl

    columns = [col for entry in spec.layer1 for col in entry.columns]
    n_genes = ht.count()

    if not columns:
        return {"n_genes": n_genes, "axes": {entry.axis: {} for entry in spec.layer1}}

    stats = ht.aggregate(
        hl.struct(
            **{
                col: hl.struct(
                    nn=hl.agg.count_where(hl.is_defined(ht[col])),
                    pos=hl.agg.count_where(ht[col] > 0),
                )
                for col in columns
            }
        )
    )

    axes: dict = {}
    for entry in spec.layer1:
        axis_manifest = {}
        for col in entry.columns:
            col_stats = stats[col]
            nn = col_stats["nn"]
            pos = col_stats["pos"]
            axis_manifest[col] = {
                "non_null_rate": (nn / n_genes) if n_genes else 0.0,
                "positive_rate": (pos / nn) if nn else None,
            }
        axes[entry.axis] = axis_manifest

    return {"n_genes": n_genes, "axes": axes}
