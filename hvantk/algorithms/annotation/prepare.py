"""Prepare a single annotation source into a per-gene, gene_id-keyed artifact.

`prepare` is Stage 1 of the annotation build (design §3): it reads one already-built
source, maps its key onto the spine's ``gene_id``, keeps only the columns a spec entry
declares, restricts to rows that land on the spine, and returns the prepared table plus a
:class:`MappingReport`. Stage 2 (compose, P3) left-joins every prepared artifact to the
spine. P2a supports ``gene_id``-keyed sources only; other key types raise until P2c.

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
        The built source, keyed on ``entry.key``.
    spine_gene_ids : Iterable[str]
        The spine's ``gene_id`` values.
    entry : hvantk.algorithms.annotation.spec.SourceEntry
    hgnc : optional
        An ``HGNCGeneCatalogStreamer``; required only for non-``gene_id`` keys (P2c).

    Returns
    -------
    (hail.Table, MappingReport)
        A ``gene_id``-keyed table with exactly ``entry.columns``, restricted to the spine,
        and the mapping report for the source.
    """
    import hail as hl

    if entry.key != "gene_id":
        raise ValueError(
            f"prepare_source (P2a) supports only gene_id-keyed sources; entry "
            f"{entry.source!r} declares key {entry.key!r}"
        )

    spine = set(spine_gene_ids)
    mapper = GeneIdMapper(hgnc, spine)

    source_ids = source_ht.gene_id.collect()
    mapping, report = mapper.from_gene_ids(source_ids, source=entry.source)

    # Keep only rows whose gene_id is on the spine, then narrow to the declared columns.
    mapped = {g for g, v in mapping.items() if v is not None}
    if not mapped:
        # hl.literal(set()) cannot infer an element type, so an empty mapping would crash
        # opaquely here. A source that reaches the spine with nothing is a hard failure:
        # surface it as a clear MappingRateError before any table is built.
        raise MappingRateError(
            f"{entry.source} mapped 0 of {report.n_in} identifiers onto the spine; "
            f"nothing to prepare. First unmapped: {', '.join(report.unmapped[:10])}"
        )

    mapped_ids = hl.literal(mapped)
    prepared = source_ht.filter(mapped_ids.contains(source_ht.gene_id))
    prepared = prepared.select(*entry.columns)

    logger.info(report.summary())
    return prepared, report
