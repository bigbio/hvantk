"""Declarative variant->gene aggregation (design decision 5's transform vocabulary).

Layering: this module imports hail + stdlib ONLY. It never imports a ``hvantk.skills`` or
``hvantk.tools`` module -- the source arrives as a Hail Table handed in by the caller (same
contract as ``prepare.py``). Identifier reconciliation onto the spine happens afterwards, in
``prepare.py``, via the existing GeneIdMapper; this module only reshapes a variant table into
one row per ``agg.by`` value.
"""
from __future__ import annotations

import logging

logger = logging.getLogger(__name__)

FRAC_GT_PREFIX = "frac_gt_"


def _missense_predicate(ht):
    """A dbNSFP row is missense iff aaref/aaalt are two distinct standard residues.

    dbNSFP's ``variant_all`` file enumerates all possible nonsynonymous SNVs; there is no VEP
    Consequence column. aaalt == 'X' is a stop-gain, aaref == 'X' a stop-loss, '.' is missing;
    all are excluded, leaving substitutions between standard amino acids.
    """
    import hail as hl

    return (
        (ht.aaref != ht.aaalt)
        & (ht.aaref != "X")
        & (ht.aaalt != "X")
        & (ht.aaref != ".")
        & (ht.aaalt != ".")
    )


# Named row predicates a spec ``aggregate.filter`` may reference. A closed registry (not an
# eval'd expression) keeps declarative config free of an injection surface; extend by adding a
# named predicate, not by growing an expression grammar.
PREDICATES = {"missense": _missense_predicate}


def _max_over_dict(dict_expr):
    """Reduce a per-transcript ``dict<transcript_id, float>`` to one scalar per variant.

    ``max`` over the dict's non-missing values. REVEL/CADD are broadcast-identical across a
    variant's transcripts (so max == the single value); MPC can differ mildly and max is a
    defensible single pick. ``hl.max`` skips missing elements (filter_missing=True default) and
    returns missing for an all-missing/empty dict.
    """
    import hail as hl

    return hl.max(dict_expr.values())


# Named transcript-dict reducers (``aggregate.reduce``). Default ``max``.
REDUCERS = {"max": _max_over_dict}


def _agg_for(token, col):
    """Map a stat token to a Hail aggregator over a scalar column expression."""
    import hail as hl

    if token == "mean":
        return hl.agg.mean(col)
    if token == "max":
        return hl.agg.max(col)
    if token == "min":
        return hl.agg.min(col)
    if token == "count":
        return hl.agg.count_where(hl.is_defined(col))
    if token.startswith(FRAC_GT_PREFIX):
        thr = float(token[len(FRAC_GT_PREFIX) :])
        return hl.agg.fraction(hl.is_defined(col) & (col > thr))
    raise ValueError(f"unknown stat token {token!r}")


def output_name(score, token):
    """Canonical feature-column name for a (score, stat) pair, e.g. revel_frac_gt_0.5."""
    return f"{score}_{token}"


def aggregate_to_gene(source_ht, agg):
    """Collapse a variant table to one row per ``agg.by`` value.

    Steps: filter (named predicate) -> reduce each score's transcript-dict to a per-variant
    scalar -> explode the ';'-delimited ``agg.by`` gene string to distinct genes -> group_by that
    gene -> the declared stats + ``n_possible_missense``. Returns a Table keyed on ``agg.by``,
    carrying exactly the ``output_name(score, token)`` columns plus ``n_possible_missense``. No
    spine reconciliation here -- the caller (prepare.py) maps ``agg.by`` onto gene_id.
    """
    import hail as hl

    ht = source_ht
    if agg.filter:
        ht = ht.filter(PREDICATES[agg.filter](ht))

    reducer = REDUCERS[agg.reduce]
    ht = ht.annotate(**{s.name: reducer(ht[s.column]) for s in agg.scores})

    # dbNSFP's Ensembl_geneid is a ';'-string (aligned with the transcript array); a variant can
    # be missense in several genes. Explode to distinct genes so each real gene gets the variant.
    ht = ht.annotate(_gene=hl.set(ht[agg.by].split(";")))
    ht = ht.explode(ht._gene)

    aggregations = {}
    for s in agg.scores:
        for token in s.stats:
            aggregations[output_name(s.name, token)] = _agg_for(token, ht[s.name])
    aggregations["n_possible_missense"] = hl.agg.count()

    grouped = ht.group_by(**{agg.by: ht._gene}).aggregate(**aggregations)
    logger.info("aggregate_to_gene: %d genes", grouped.count())
    return grouped
