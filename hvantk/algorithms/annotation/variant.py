"""Variant-level transforms for ML training sets.

Selecting per-variant columns from an allele-keyed source (e.g. dbNSFP rankscores and
allele frequencies), turning protein-site tables into locus sets (PTM codons), and
choosing one gene for a variant that a source assigns to several.

Layering: Hail plus ``hvantk.algorithms`` siblings only; sources arrive as tables handed in
by the caller, never read here.
"""

from __future__ import annotations

from typing import Sequence

import hail as hl

from hvantk.algorithms.annotation.spec import VariantColumn


def field_expr(ht: hl.Table, path: str) -> hl.Expression:
    """Resolve a dotted field path (``gnomAD.gnomAD_exomes_AF``) on ``ht``'s rows.

    Raises
    ------
    KeyError
        Naming the first segment that does not exist and the fields available there.
    """
    expr = ht.row
    walked = []
    for segment in path.split("."):
        fields = list(expr.dtype) if isinstance(expr.dtype, hl.tstruct) else []
        if segment not in fields:
            where = ".".join(walked) or "the row"
            raise KeyError(
                f"field {path!r}: {segment!r} is not a field of {where} "
                f"(available: {', '.join(fields[:20])}{', ...' if len(fields) > 20 else ''})"
            )
        expr = expr[segment]
        walked.append(segment)
    return expr


def _cast(expr: hl.Expression, cast: str | None) -> hl.Expression:
    if cast is None:
        return expr
    if cast != "float64":
        raise ValueError(f"unsupported cast {cast!r}; expected 'float64'")
    # dbNSFP ships many numbers as strings, with '.' for "no value".
    return hl.parse_float64(expr) if expr.dtype == hl.tstr else hl.float64(expr)


def column_expr(ht: hl.Table, column: VariantColumn) -> hl.Expression:
    """The value of one allele-join column on ``ht``'s rows."""
    if column.fields:
        if column.reduce != "max":
            raise ValueError(f"column {column.name!r}: unsupported reduce {column.reduce!r}")
        exprs = [_cast(field_expr(ht, f), column.cast) for f in column.fields]
        return hl.max(*exprs, filter_missing=True)
    return _cast(field_expr(ht, column.field), column.cast)


def allele_columns(source: hl.Table, columns: Sequence[VariantColumn]) -> hl.Table:
    """``source`` reduced to one field per output column, keyed as ``source`` is."""
    return source.select(**{c.name: column_expr(source, c) for c in columns})


def codon_sites(
    sites: hl.Table,
    *,
    chrom: str = "chrom",
    start: str = "codon_start",
    end: str = "codon_end",
    reference_genome: str = "GRCh38",
) -> hl.Table:
    """The loci covered by each site's codon, as a distinct ``locus``-keyed table.

    Contigs are written with a ``chr`` prefix and the mitochondrion (``MT``/``M``) as
    ``chrM``, as GRCh38 names them. A codon whose ends are two bases apart contributes all
    three positions; any other span is a codon split by an intron, which contributes its two
    ends only, never the intron between them. Positions that are not valid loci are dropped.
    """
    t = sites.key_by()
    lo = hl.min(t[start], t[end])
    hi = hl.max(t[start], t[end])
    c = t[chrom]
    contig = (
        hl.case()
        .when(c.startswith("chr"), c)
        .when((c == "MT") | (c == "M"), "chrM")
        .default("chr" + c)
    )
    t = t.select(
        _contig=contig,
        _pos=hl.if_else(hi - lo == 2, hl.range(lo, hi + 1), hl.array([lo, hi])),
    )
    t = t.explode("_pos")
    t = t.filter(hl.is_valid_locus(t._contig, t._pos, reference_genome=reference_genome))
    t = t.select(locus=hl.locus(t._contig, t._pos, reference_genome=reference_genome))
    return t.key_by("locus").distinct()


def choose_candidate(
    candidates: hl.ArrayExpression, preferred: hl.SetExpression
) -> hl.StructExpression:
    """Pick one of a variant's candidate genes: ``struct(index: int32, rule: str)``.

    ``candidates`` is the variant's gene ids in the source's listed order (an element may
    be missing when an id did not map); ``preferred`` is the set of gene ids a second
    source assigns to the variant (missing = none). ``rule`` is

    - ``"none"``: no candidates (``index`` missing);
    - ``"single"``: one candidate (index 0), whatever ``preferred`` holds;
    - ``"preferred"``: several candidates; the first defined one in ``preferred``;
    - ``"first"``: several candidates, none in ``preferred``; index 0.
    """
    preferred = hl.or_else(preferred, hl.empty_set(hl.tstr))
    n = hl.len(candidates)
    hits = hl.range(n).filter(
        lambda i: hl.is_defined(candidates[i]) & preferred.contains(candidates[i])
    )
    return (
        hl.case()
        .when(n == 0, hl.struct(index=hl.missing(hl.tint32), rule="none"))
        .when(n == 1, hl.struct(index=0, rule="single"))
        .when(hl.len(hits) > 0, hl.struct(index=hits[0], rule="preferred"))
        .default(hl.struct(index=0, rule="first"))
    )
