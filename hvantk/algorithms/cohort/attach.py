"""Attach a cohort onto the Layer-1 gene feature matrix.

``attach`` maps the cohort table's gene key onto ``gene_id`` (reusing Stage-1's
identifier reconciliation and its mapping-rate report), then left-joins the cohort's
declared columns onto the Layer-1 matrix over the FULL spine universe.

Two properties matter and are the reason this is not a plain merge:

* **Nothing is imputed.** A gene absent from the cohort table keeps missing values in
  every cohort column. It is not zero, and it is not "tested and null".
* **Absence is recorded.** ``cohort_tested`` is true exactly for genes that had a row in
  the cohort table, so "never tested by this cohort" stays distinguishable from "tested,
  result was null". The analysis universe is then a downstream choice rather than an
  accident of which rows the input file happened to contain.

Layering: imports Hail, sibling ``hvantk.algorithms`` modules and ``core`` only -- never
``hvantk.skills`` or ``hvantk.tools``. The HGNC catalog is passed in by the caller.
"""
from __future__ import annotations

import logging

logger = logging.getLogger(__name__)

#: The full universe of column names ``attach`` itself can introduce on the output
#: table. Not all of these are reserved for every manifest -- ``cohort_tested`` always
#: is, but ``label`` is only written (and therefore only reserved) when the manifest
#: declares labels; see :func:`_reserved_columns` for the per-manifest active set that
#: :func:`check_no_layer1_collisions` actually checks against.
RESERVED_COLUMNS = ("cohort_tested", "label")


def as_source_entry(manifest):
    """Express the cohort manifest as a Stage-1 ``SourceEntry``.

    A cohort table is just another source that has to be reconciled onto the spine, so
    it reuses ``prepare_source`` rather than duplicating alias resolution and rate
    reporting. Every declared column travels as one entry because a cohort supplies a
    single table for all of its axes.
    """
    from hvantk.algorithms.annotation.spec import SourceEntry

    return SourceEntry(
        axis=manifest.name,
        source=manifest.name,
        key=manifest.key,
        columns=manifest.declared_columns(),
        min_mapping_rate=manifest.min_mapping_rate,
    )


def _expose_key_column(cohort_ht, manifest):
    """Expose the cohort table's identifier column under the name ``prepare_source`` expects.

    ``prepare_source`` (Stage 1, ``hvantk.algorithms.annotation.prepare``) resolves the
    identifier column by plain field-name lookup on ``entry.key`` (``source_ht[key]``)
    -- it never consults ``source_ht``'s actual Hail key. ``manifest.key`` names the
    identifier SPACE (``gene_id`` / ``hgnc_id`` / ``symbol``), which is what
    ``prepare_source`` dispatches the mapper on; ``manifest.key_column`` names the
    column that actually holds those values in the cohort's own table, and the two are
    often different (a real cohort's symbol column is commonly named ``gene``). When
    they differ, rename the column so the field ``prepare_source`` looks up exists; when
    they match -- every manifest written before ``key_column`` existed -- this is a
    no-op, which is what keeps old manifests working unchanged.
    """
    if manifest.key_column == manifest.key:
        return cohort_ht
    return cohort_ht.rename({manifest.key_column: manifest.key})


def _reserved_columns(manifest) -> set:
    """The reserved output names ``attach()`` will actually write for this manifest.

    ``cohort_tested`` is unconditional -- every call writes it. ``label`` is written
    only when ``manifest.labels is not None`` (see the ``if manifest.labels is not
    None`` branch in :func:`attach`); reserving it unconditionally would reject a
    cohort's unrelated ``label`` column (a common name in a case/control table) for a
    collision that can never actually happen for that manifest.
    """
    reserved = {"cohort_tested"}
    if manifest.labels is not None:
        reserved.add("label")
    return reserved


def check_no_layer1_collisions(layer1, manifest) -> None:
    """Raise if a declared cohort column, or a reserved output name, would overwrite
    an existing column -- on either side of the join.

    ``Table.annotate`` silently replaces an existing field, so without this a cohort
    column named like a Layer-1 feature would quietly shadow it and every downstream
    number would be computed on the wrong values. Schema introspection only -- reading
    field names starts no Hail job -- so this runs before any computation.

    Three checks, in order:

    1. A declared cohort column named ``cohort_tested`` or (when the manifest
       declares labels) ``label`` -- the reserved names ``attach`` assigns itself
       after the left-join. Without this, the cohort's own data would be silently
       overwritten by the synthetic value with no exception and no warning.
    2. A reserved name that ALREADY exists in the Layer-1 schema. Without this,
       ``layer1.annotate(cohort_tested=..., label=...)`` replaces the existing
       Layer-1 field outright -- e.g. a prior `attach` run's ``cohort_tested`` column,
       or a Layer-1 source that happens to carry a ``label`` column of its own -- with
       no error, exit 0. This is the same failure class as (1), just checked against
       a schema ``attach`` itself is about to write rather than one the manifest
       declares.
    3. A declared cohort column that already exists in the Layer-1 matrix.
    """
    declared = set(manifest.declared_columns())
    reserved = _reserved_columns(manifest)

    reserved_clashing = sorted(declared & reserved)
    if reserved_clashing:
        raise ValueError(
            f"cohort {manifest.name!r} declares column(s) {', '.join(reserved_clashing)} "
            f"that collide with attach()'s reserved output column name(s) "
            f"({', '.join(sorted(reserved))}); rename the cohort column(s) -- attach "
            "would silently overwrite the cohort's data with its own computed value"
        )

    existing = set(layer1.row)

    reserved_in_layer1 = sorted(reserved & existing)
    if reserved_in_layer1:
        raise ValueError(
            f"the Layer-1 matrix already has column(s) {', '.join(reserved_in_layer1)}, "
            f"which cohort {manifest.name!r} would overwrite with attach()'s own "
            f"computed value ({', '.join(sorted(reserved))} are reserved output "
            "names); rename or drop the existing Layer-1 column(s) first -- attaching "
            "would silently replace them, including a prior attach run's own output"
        )

    clashing = sorted(declared & existing)
    if clashing:
        raise ValueError(
            f"cohort {manifest.name!r} declares column(s) {', '.join(clashing)} "
            "that already exist in the Layer-1 matrix; rename the cohort column(s) "
            "-- attaching would silently overwrite the Layer-1 values"
        )


def _check_no_duplicate_cohort_rows(cohort_ht, manifest) -> None:
    """Raise if the cohort table has more than one row for the same gene key.

    ``prepare_source``'s ``key == "gene_id"`` branch (Stage-1 machinery, shared with
    other pipelines and never modified here) does not enforce one row per gene -- it
    only filters onto the spine and selects columns. Without this check, the later
    index-join (``prepared[layer1.gene_id]``) picks one of several matching rows
    arbitrarily and silently drops the rest, in violation of the contract's R3 rule
    ("one row per tested gene"). Checked here, on the cohort side, precisely because
    ``prepare_source`` must not change behaviour for its other callers.

    Cheap by construction -- a row count and a distinct-key count, no collect -- with
    the actual duplicated keys materialised only once those two counts have already
    proven a duplicate exists.
    """
    import hail as hl

    key = manifest.key
    n_rows = cohort_ht.count()
    n_keys = cohort_ht.key_by(key).distinct().count()
    if n_rows == n_keys:
        return

    counts = cohort_ht.aggregate(hl.agg.counter(cohort_ht[key]))
    dup_keys = sorted(k for k, c in counts.items() if c > 1)
    raise ValueError(
        f"cohort {manifest.name!r} table has {n_rows - n_keys} duplicate row(s) for "
        f"gene key(s) {', '.join(dup_keys[:10])}; the cohort table must carry "
        "exactly one row per gene -- aggregate the table before attaching"
    )


def attach(layer1, cohort_ht, manifest, *, hgnc=None):
    """Left-join a cohort's declared columns onto the Layer-1 matrix.

    Parameters
    ----------
    layer1 : hail.Table
        The composed Layer-1 matrix, keyed by ``gene_id``. It already spans the whole
        spine (``compose`` emits one row per spine gene), so no separate spine argument
        is needed -- and passing one would admit a mismatch.
    cohort_ht : hail.Table
        The cohort's gene-level table, keyed by whatever ``manifest.key_column`` names
        (defaults to ``manifest.key`` when the manifest omits it).
    manifest : hvantk.algorithms.cohort.spec.CohortManifest
    hgnc : optional
        An HGNC gene-catalog streamer, required when ``manifest.key`` is ``hgnc_id`` or
        ``symbol``. Constructed by the caller (the CLI) so this module never imports
        the skills layer.

    Returns
    -------
    (hail.Table, dict)
        The attached table -- Layer-1 columns, the cohort's declared columns, and a
        ``cohort_tested`` bool -- and a report carrying the mapping rate, the tested
        count, the prior's declared direction, and per-axis per-column rates.

    Raises
    ------
    ValueError
        If a declared cohort column collides with a Layer-1 column or with a
        reserved output column name (``RESERVED_COLUMNS``), if the cohort table
        carries more than one row for the same gene key, or if the manifest
        declares labels but no ``hgnc`` was passed.
    hvantk.algorithms.annotation.mapping.MappingRateError
        If the gene-key mapping rate falls below ``manifest.min_mapping_rate``.
    """
    import hail as hl

    from hvantk.algorithms.annotation.mapping import enforce_rate
    from hvantk.algorithms.annotation.prepare import prepare_source

    check_no_layer1_collisions(layer1, manifest)

    cohort_ht = _expose_key_column(cohort_ht, manifest)
    _check_no_duplicate_cohort_rows(cohort_ht, manifest)

    spine_gene_ids = layer1.gene_id.collect()
    prepared, report = prepare_source(
        cohort_ht, spine_gene_ids, as_source_entry(manifest), hgnc=hgnc
    )
    enforce_rate(report, manifest.min_mapping_rate)
    logger.info("cohort %s: %s", manifest.name, report.summary())

    # Index-join: a missing struct means the gene had no row in the cohort table, which
    # is what makes `cohort_tested` correct even when a matched row's own values are null.
    joined = prepared[layer1.gene_id]
    updates = {col: joined[col] for col in manifest.declared_columns()}
    updates["cohort_tested"] = hl.is_defined(joined)

    label_report = None
    if manifest.labels is not None:
        label_ids, label_report = _resolve_labels(manifest, spine_gene_ids, hgnc=hgnc)
        updates["label"] = hl.literal(label_ids, dtype=hl.tset(hl.tstr)).contains(
            layer1.gene_id
        )

    attached = layer1.annotate(**updates)
    return attached, _build_report(attached, manifest, report, label_report)


def _resolve_labels(manifest, spine_gene_ids, *, hgnc):
    """Map the declared label gene set onto spine ``gene_id``s.

    Label sets are symbol-keyed -- ``hvantk genesets prepare`` rejects Ensembl and
    Entrez identifiers outright -- so they always need reconciliation before they can
    meet a ``gene_id``-keyed matrix. The rate is enforced because silently dropping a
    slice of a label set corrupts every downstream metric while still producing
    plausible-looking output. An informally curated clinical panel is exactly the input
    where stale symbols are likely.
    """
    from hvantk.algorithms.annotation.mapping import GeneIdMapper, enforce_rate

    from hvantk.algorithms.cohort.labels import load_label_genes

    if hgnc is None:
        raise ValueError(
            f"cohort {manifest.name!r} declares labels, which are always "
            "symbol-keyed and need the HGNC streamer; pass hgnc=<HGNCGeneCatalogStreamer>"
        )

    genes = load_label_genes(manifest.labels)
    mapper = GeneIdMapper(hgnc, spine_gene_ids)
    mapping, report = mapper.from_symbols(genes, source=f"{manifest.name}:labels")
    enforce_rate(report, manifest.labels.min_mapping_rate)
    logger.info("cohort %s labels: %s", manifest.name, report.summary())

    return {g for g in mapping.values() if g is not None}, report


def _build_report(ht, manifest, mapping_report, label_report=None) -> dict:
    """Coverage report: tested count, mapping rate, and per-column rates.

    Mirrors the shape ``compose`` emits so both layers' manifests read the same way.

    ``positive_rate`` needs an ordered comparison (``column > 0``), so it is computed
    only for declared columns Hail's own schema reports as numeric
    (:func:`hail.is_numeric` on ``ht[col].dtype``). A declared column that
    ``hl.import_table(..., impute=True)`` typed as ``str`` -- a single dirty row like
    ``<0.001``, ``.``, ``Inf``, or a stray space is enough -- gets
    ``positive_rate: None`` and a logged warning naming the column and its type,
    rather than crashing here with a bare ``TypeError`` after mapping, joining and
    label resolution have already finished.
    """
    import hail as hl

    columns = list(manifest.declared_columns())
    n_genes = ht.count()

    numeric_columns = {col for col in columns if hl.is_numeric(ht[col].dtype)}
    for col in columns:
        if col not in numeric_columns:
            logger.warning(
                "cohort %s: declared column %r has non-numeric type %s; "
                "positive_rate will be reported as null for this column",
                manifest.name,
                col,
                ht[col].dtype,
            )

    stats = ht.aggregate(
        hl.struct(
            n_tested=hl.agg.count_where(ht.cohort_tested),
            **{
                col: hl.struct(
                    nn=hl.agg.count_where(hl.is_defined(ht[col])),
                    **(
                        {"pos": hl.agg.count_where(ht[col] > 0)}
                        if col in numeric_columns
                        else {}
                    ),
                )
                for col in columns
            },
        )
    )

    def _rates(col):
        nn = stats[col]["nn"]
        non_null_rate = (nn / n_genes) if n_genes else 0.0
        if col not in numeric_columns:
            return {"non_null_rate": non_null_rate, "positive_rate": None}
        pos = stats[col]["pos"]
        return {
            "non_null_rate": non_null_rate,
            "positive_rate": (pos / nn) if nn else None,
        }

    report = {
        "cohort": manifest.name,
        "n_genes": n_genes,
        "n_tested": stats["n_tested"],
        "key": manifest.key,
        "key_column": manifest.key_column,
        "mapping_rate": mapping_report.rate,
        "prior": {
            "column": manifest.prior.column,
            "direction": manifest.prior.direction,
            **_rates(manifest.prior.column),
        },
        "axes": {
            entry.axis: {col: _rates(col) for col in entry.columns}
            for entry in manifest.cohort_axes
        },
    }

    if label_report is not None:
        report["labels"] = {
            "gene_set": manifest.labels.gene_set,
            "set_name": manifest.labels.set_name,
            "n_in": label_report.n_in,
            "n_mapped": label_report.n_mapped,
            "mapping_rate": label_report.rate,
        }

    return report
