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
    if key == "gene_id":
        mapping, report = mapper.from_gene_ids(source_keys, source=entry.source)
    elif key == "hgnc_id":
        mapping, report = mapper.from_hgnc_ids(source_keys, source=entry.source)
    else:  # symbol
        mapping, report = mapper.from_symbols(source_keys, source=entry.source)

    resolved = {k: v for k, v in mapping.items() if v is not None}
    if not resolved:
        raise MappingRateError(
            f"{entry.source} mapped 0 of {report.n_in} identifiers onto the spine; "
            f"nothing to prepare. First unmapped: {', '.join(report.unmapped[:10])}"
        )

    if key == "gene_id":
        # Already gene_id-keyed: filter to the mapped ids and select columns.
        mapped_ids = hl.literal(set(resolved.values()))
        prepared = source_ht.filter(mapped_ids.contains(source_ht.gene_id))
        prepared = prepared.select(*entry.columns)
    else:
        # Re-key: attach gene_id from the mapping, drop rows that did not resolve, key on gene_id.
        lut = hl.literal(resolved)  # source_key -> gene_id
        prepared = source_ht.annotate(gene_id=lut.get(source_ht[key]))
        prepared = prepared.filter(hl.is_defined(prepared.gene_id))
        prepared = prepared.key_by("gene_id").select(*entry.columns)
        # One row per gene is a hard contract; a many-to-one mapping would duplicate a gene.
        n_rows = prepared.count()
        n_genes = prepared.distinct().count()
        if n_rows != n_genes:
            raise ValueError(
                f"{entry.source}: {n_rows - n_genes} source rows resolve to the same "
                f"gene_id (many-to-one); this source needs aggregation before prepare, not "
                f"available until a later increment"
            )

    logger.info(report.summary())
    return prepared, report
