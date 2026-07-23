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
    if key not in ("gene_id", "hgnc_id", "symbol"):
        raise ValueError(
            f"prepare_source supports gene_id, hgnc_id, symbol keys; entry "
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
            source_ht, key, resolved, entry.columns, entry.source
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
    else:
        raise ValueError(f"unsupported id-space {to_space!r}")
    resolved = {k: v for k, v in mapping.items() if v is not None}
    if not resolved:
        raise MappingRateError(
            f"{source_label} mapped 0 of {report.n_in} identifiers onto the spine; "
            f"nothing to prepare. First unmapped: {', '.join(report.unmapped[:10])}"
        )
    return resolved, report


def _rekey_onto_gene_id(source_ht, key_col, resolved, columns, source_label):
    """Attach gene_id from ``resolved``, drop unresolved, key on gene_id, select columns.

    Enforces one row per gene: a many-to-one mapping (or an un-aggregated source) raises.
    """
    import hail as hl

    lut = hl.literal(resolved)
    prepared = source_ht.annotate(gene_id=lut.get(source_ht[key_col]))
    prepared = prepared.filter(hl.is_defined(prepared.gene_id))
    prepared = prepared.key_by("gene_id").select(*columns)
    n_rows = prepared.count()
    n_genes = prepared.distinct().count()
    if n_rows != n_genes:
        raise ValueError(
            f"{source_label}: {n_rows - n_genes} source rows resolve to the same gene_id "
            f"(many-to-one); this source needs aggregation before prepare"
        )
    return prepared


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
        grouped, agg.by, resolved, entry.columns, entry.source
    )
    logger.info(report.summary())
    return prepared, report
