"""Prepare a single annotation source into a per-gene, gene_id-keyed artifact.

`prepare` is Stage 1 of the annotation build (design §3): it reads one already-built
source, maps its key onto the spine's ``gene_id``, keeps only the columns a spec entry
declares, restricts to rows that land on the spine, and returns the prepared table plus a
:class:`MappingReport`. Stage 2 (compose, P3) left-joins every prepared artifact to the
spine. Sources keyed on ``gene_id``, ``hgnc_id``, or ``symbol`` are supported; non-``gene_id``
keys are re-keyed onto ``gene_id`` via the HGNC-backed mapper before column selection.

Layering: this module reads sources as Hail Tables handed in by the caller (the CLI loads
them by path). It must never import a ``hvantk.skills`` module -- the source is data, not
code.
"""
from __future__ import annotations

import logging

from hvantk.algorithms.annotation.mapping import GeneIdMapper, MappingRateError

logger = logging.getLogger(__name__)

# Gene-level key spaces a prepared source may declare. Single source of truth for the
# guard below; `_resolve_to_gene_id` dispatches the same set, and the feature-spec schema
# advertises it. Keep the three in step -- they drifted once already.
_GENE_KEY_SPACES = ("gene_id", "hgnc_id", "symbol", "uniprot_id")


def _collapse_reducers():
    """Reductions available to a spec's ``collapse:`` key, built lazily (needs Hail).

    Mapping a source key onto ``gene_id`` can be legitimately many-to-one -- one gene
    routinely has several UniProt accessions (isoforms, historical entries), so both
    ``insider:interfaces`` and ``uniprot-ptm:sites`` collapse. Without a declared policy
    that is an error, because silently dropping or double-counting rows is worse than
    failing. ``max`` is the safe pick for accession collapse: the rows describe the SAME
    protein, so it never inflates the way ``sum`` would.
    """
    import hail as hl

    return {"max": hl.agg.max, "min": hl.agg.min, "mean": hl.agg.mean,
            "sum": hl.agg.sum}


class _LazyCollapseReducers(dict):
    """Populates from Hail on first use so importing this module stays Hail-free."""

    def _ensure(self):
        if not super().__len__():
            self.update(_collapse_reducers())

    def __contains__(self, key):
        self._ensure()
        return super().__contains__(key)

    def __getitem__(self, key):
        self._ensure()
        return super().__getitem__(key)

    def __iter__(self):
        self._ensure()
        return super().__iter__()

    def __len__(self):
        self._ensure()
        return super().__len__()


COLLAPSE_REDUCERS = _LazyCollapseReducers()


def prepare_source(source_ht, spine_gene_ids, entry, *, hgnc=None):
    """Map ``source_ht`` onto the spine and select ``entry.columns``.

    Parameters
    ----------
    source_ht : hail.Table
        The built source, keyed on ``entry.key`` (``gene_id``, ``hgnc_id``, or ``symbol``).
    spine_gene_ids : Iterable[str]
        The spine's ``gene_id`` values.
    entry : hvantk.algorithms.annotation.spec.SourceEntry
    hgnc : optional
        An ``HGNCGeneCatalogStreamer``; required for non-``gene_id`` keys.

    Returns
    -------
    (hail.Table, MappingReport)
        A ``gene_id``-keyed table with exactly ``entry.columns``, restricted to the spine,
        and the mapping report for the source.
    """
    import hail as hl

    key = entry.key
    # Keep in step with the schema's key enum and `_resolve_to_gene_id`'s dispatch --
    # this listed the accepted spaces a third time and so silently excluded uniprot_id,
    # which both the schema and the mapper already supported.
    if key not in _GENE_KEY_SPACES:
        raise ValueError(
            f"prepare_source supports {', '.join(_GENE_KEY_SPACES)} keys; entry "
            f"{entry.source!r} declares key {key!r}"
        )
    if key != "gene_id" and hgnc is None:
        raise ValueError(
            f"entry {entry.source!r} has key {key!r}, which needs the HGNC streamer; "
            f"pass hgnc=<HGNCGeneCatalogStreamer>"
        )

    spine = set(spine_gene_ids)
    mapper = GeneIdMapper(hgnc, spine)

    source_keys = source_ht[key].collect()
    resolved, report = _resolve_to_gene_id(source_keys, key, entry.source, mapper)

    if key == "gene_id":
        mapped_ids = hl.literal(set(resolved.values()))
        prepared = source_ht.filter(mapped_ids.contains(source_ht.gene_id))
        prepared = prepared.select(*entry.columns)
    else:
        prepared = _rekey_onto_gene_id(
            source_ht, key, resolved, entry.columns, entry.source,
            collapse=entry.collapse,
        )

    logger.info(report.summary())
    return prepared, report


def _resolve_to_gene_id(source_keys, to_space, source_label, mapper):
    """Dispatch the mapper by id-space and drop unmapped; raise if nothing maps."""
    if to_space == "gene_id":
        mapping, report = mapper.from_gene_ids(source_keys, source=source_label)
    elif to_space == "hgnc_id":
        mapping, report = mapper.from_hgnc_ids(source_keys, source=source_label)
    elif to_space == "symbol":
        mapping, report = mapper.from_symbols(source_keys, source=source_label)
    elif to_space == "uniprot_id":
        mapping, report = mapper.from_uniprot_ids(source_keys, source=source_label)
    else:
        raise ValueError(f"unsupported id-space {to_space!r}")
    resolved = {k: v for k, v in mapping.items() if v is not None}
    if not resolved:
        raise MappingRateError(
            f"{source_label} mapped 0 of {report.n_in} identifiers onto the spine; "
            f"nothing to prepare. First unmapped: {', '.join(report.unmapped[:10])}"
        )
    return resolved, report


def _rekey_onto_gene_id(source_ht, key_col, resolved, columns, source_label, *,
                        collapse=None):
    """Attach gene_id from ``resolved``, drop unresolved, key on gene_id, select columns.

    Enforces one row per gene. A many-to-one mapping raises unless the spec entry
    declares a ``collapse`` reduction -- see :data:`COLLAPSE_REDUCERS`.
    """
    import hail as hl

    # Unkey first: when `aggregate.by` is itself `gene_id` (pQTL and GTEx eQTL are both
    # gene_id-keyed), group_by has already keyed the table on gene_id and annotating it
    # raises "cannot overwrite key field". Dropping the key is safe -- the very next step
    # re-keys on gene_id anyway.
    source_ht = source_ht.key_by()

    lut = hl.literal(resolved)
    prepared = source_ht.annotate(gene_id=lut.get(source_ht[key_col]))
    prepared = prepared.filter(hl.is_defined(prepared.gene_id))
    prepared = prepared.key_by("gene_id").select(*columns)
    n_rows = prepared.count()
    n_genes = prepared.distinct().count()
    if n_rows == n_genes:
        return prepared

    if collapse is None:
        raise ValueError(
            f"{source_label}: {n_rows - n_genes} source rows resolve to the same gene_id "
            f"(many-to-one); aggregate the source before prepare, or declare a "
            f"`collapse:` reduction ({', '.join(sorted(COLLAPSE_REDUCERS))})"
        )
    if collapse not in COLLAPSE_REDUCERS:
        raise ValueError(
            f"{source_label}: unknown collapse reduction {collapse!r}; "
            f"expected one of {', '.join(sorted(COLLAPSE_REDUCERS))}"
        )
    agg = COLLAPSE_REDUCERS[collapse]
    logger.info(
        "%s: collapsing %d many-to-one rows onto gene_id with %r",
        source_label, n_rows - n_genes, collapse,
    )
    return prepared.group_by(prepared.gene_id).aggregate(
        **{c: agg(prepared[c]) for c in columns}
    )


def prepare_matrix_source(matrix_ad, spine_gene_ids, entry, *, hgnc=None):
    """Reduce an expression-matrix AnnData to per-gene features, then reconcile onto the spine.

    ``entry.matrix`` declares the reduction (group axis, per-group stats, EWCE specificity). The
    reduced table is symbol-keyed; it rides the existing symbol resolve step onto gene_id, then
    collapses symbol-collisions (two symbols -> one gene) by max rather than raising, because a
    matrix is a matrix->gene reduction (see :func:`_collapse_matrix_onto_gene_id`).
    """
    import hail as hl

    from hvantk.algorithms.annotation import matrix as matrix_mod

    if entry.matrix is None:
        raise ValueError(
            f"entry {entry.source!r} prepared as matrix but has no matrix block"
        )
    if hgnc is None:
        raise ValueError(
            f"entry {entry.source!r} is a symbol-keyed matrix source; needs the HGNC streamer"
        )

    df = matrix_mod.reduce_matrix_to_gene(
        matrix_ad, entry.matrix
    )  # 'symbol' + feature cols
    # A matrix source's column set is decided at RUNTIME by the atlas's groups -- with
    # specificity emit="vector" (the default) there is one column per surviving group,
    # whose names come from _san() over labels the spec author cannot enumerate ahead of
    # time. So carry everything the reducer produced, rather than the static
    # entry.columns: restricting to the declared list here silently dropped the whole
    # vector on the floor, leaving emit="vector" a no-op through the real pipeline.
    produced = [c for c in df.columns if c != "symbol"]
    # entry.columns stays meaningful as an assertion: whatever the spec explicitly names
    # MUST be produced. That catches a typo, a renamed group, or an atlas swapped under
    # the spec -- all of which would otherwise surface as a silently absent feature.
    missing = [c for c in entry.columns if c not in produced]
    if missing:
        raise ValueError(
            f"entry {entry.source!r} declares column(s) {missing} that the matrix "
            f"reducer did not produce; it emitted {produced}. Check the spec's "
            f"specificity.name / stats against the atlas's group labels"
        )
    grouped = hl.Table.from_pandas(df, key=["symbol"])
    mapper = GeneIdMapper(hgnc, set(spine_gene_ids))
    resolved, report = _resolve_to_gene_id(
        grouped.symbol.collect(), "symbol", entry.source, mapper
    )
    prepared = _collapse_matrix_onto_gene_id(
        grouped, "symbol", resolved, produced, entry.source
    )
    logger.info(report.summary())
    return prepared, report


def _collapse_matrix_onto_gene_id(source_ht, key_col, resolved, columns, source_label):
    """Attach gene_id, drop unresolved, and collapse symbol-collisions onto one row per gene.

    Two distinct symbols can resolve to one ``gene_id`` (HGNC aliases / previous symbols). For a
    plain gene table that is an error -- there is no combine rule -- so :func:`_rekey_onto_gene_id`
    raises. An expression matrix, however, is a matrix->gene reduction: colliding symbols are
    combined by ``max``, matching the reducer's own duplicate-symbol collapse
    (:func:`hvantk.algorithms.annotation.matrix.reduce_matrix_to_gene`). All P2c-4 matrix features
    (EWCE specificity, fraction expressed, mean expression) are "larger = more signal", so taking
    the max keeps the strongest evidence when one gene carries two symbol rows.
    """
    import hail as hl

    lut = hl.literal(resolved)
    attached = source_ht.annotate(gene_id=lut.get(source_ht[key_col]))
    attached = attached.filter(hl.is_defined(attached.gene_id))
    collapsed = attached.group_by(attached.gene_id).aggregate(
        **{c: hl.agg.max(attached[c]) for c in columns}
    )
    return collapsed


def prepare_variant_source(source_ht, spine_gene_ids, entry, *, hgnc=None):
    """Aggregate a variant source to per-gene stats, then reconcile onto the spine.

    ``entry.aggregate`` declares the group column, id-space, filter, and per-score stats. The
    aggregated table is keyed on ``aggregate.by`` (an id in ``aggregate.to`` space); it is mapped
    onto gene_id exactly like a gene-keyed source.
    """
    from hvantk.algorithms.annotation import transforms

    agg = entry.aggregate
    if agg is None:
        raise ValueError(
            f"entry {entry.source!r} has key 'variant' but no aggregate block"
        )
    if agg.to != "gene_id" and hgnc is None:
        raise ValueError(
            f"entry {entry.source!r} aggregates to {agg.to!r}, which needs the HGNC streamer"
        )

    grouped = transforms.aggregate_to_gene(source_ht, agg)  # keyed on agg.by
    mapper = GeneIdMapper(hgnc, set(spine_gene_ids))
    source_keys = grouped[agg.by].collect()
    resolved, report = _resolve_to_gene_id(source_keys, agg.to, entry.source, mapper)
    prepared = _rekey_onto_gene_id(
        grouped, agg.by, resolved, entry.columns, entry.source, collapse=entry.collapse
    )
    logger.info(report.summary())
    return prepared, report
