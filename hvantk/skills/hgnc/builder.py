"""Hail Table builder for HGNC gene nomenclature data.

Owns the ``build_hgnc_gene_lookup`` builder. Turns the HGNC
complete-set TSV into an ``AnnotationTable`` keyed by ``hgnc_id``.
"""

from __future__ import annotations

import logging

import hail as hl

from hvantk.skills.hgnc.shared.constants import (
    HGNC_GENE_FIELDS,
    HGNC_PIPE_SEPARATED_FIELDS,
)
from hvantk.core.utils.table_utils import get_row_fields

logger = logging.getLogger(__name__)


def _missing_declared_fields(row_fields) -> list[str]:
    """Declared HGNC_GENE_FIELDS keys absent from the input header.

    Pure helper (no Hail) so a declared field upstream silently drops --
    e.g. location_sortable in #381, caught only by the drift probe -- shows
    up as a warning instead of being tolerated forever by the rename map.
    """
    return [k for k in HGNC_GENE_FIELDS if k not in row_fields]


def _split_pipe_tokens(expr: hl.expr.StringExpression) -> hl.expr.ArrayExpression:
    """Tokenize one HGNC pipe-separated cell into its de-quoted, stripped parts.

    hgnc_complete_set.txt QUOTES its multi-value fields (prev_symbol is written
    `"H1F4|HIST1H1E"`), so the quotes must come off BEFORE the split -- otherwise
    the first and last elements keep a stray `"` and every alias/prev-symbol
    lookup for the clean symbol misses. No HGNC value legitimately contains a
    double quote, so stripping them all is safe and simpler than anchoring to
    the ends. Empty/missing input yields an empty array rather than `[""]`.

    Shared by the plain `array<str>` fields in HGNC_PIPE_SEPARATED_FIELDS and
    by lsdb, which tokenizes the same way but then re-pairs the tokens into
    `array<struct{name, url}>` instead of keeping them flat (#391).
    """
    return hl.if_else(
        hl.is_defined(expr) & (expr != ""),
        expr.replace('"', "")
        .split("\\|")
        .map(lambda x: x.strip())
        .filter(lambda x: x != ""),
        hl.empty_array(hl.tstr),
    )


def build_hgnc_gene_lookup(
    parsed_input,
    ctx,
    *,
    include_withdrawn: bool = False,
    fields=None,
    strict_fields: bool = False,
):
    """Plugin builder — returns an AnnotationTable.

    Parameters
    ----------
    parsed_input : str | Path
        Path to the HGNC complete-set TSV.
    ctx : hvantk.core.models.BuildContext
        Platform-provided context; supplies provenance.
    include_withdrawn : bool, optional
        If True, include withdrawn genes (default False).
    fields : list of str, optional
        Optional field subset; if None, all fields are kept.
    strict_fields : bool, optional
        If True, raise ValueError when a field declared in HGNC_GENE_FIELDS is
        absent from the input header, instead of only logging a warning
        (default False).

    Returns
    -------
    hvantk.core.models.AnnotationTable
        The HGNC table wrapped with Provenance.

    Raises
    ------
    ValueError
        If ``strict_fields`` is True and the input header is missing one or
        more fields declared in HGNC_GENE_FIELDS.
    """
    from hvantk.core.models import AnnotationTable

    ht = hl.import_table(
        paths=str(parsed_input),
        impute=False,
        min_partitions=10,
        missing="",
    )

    # Rename fields to standardized names
    row_fields = get_row_fields(ht)
    missing_fields = _missing_declared_fields(row_fields)
    if missing_fields:
        if strict_fields:
            raise ValueError(
                "HGNC input is missing declared field(s): " + ", ".join(missing_fields)
            )
        logger.warning(
            "HGNC input is missing declared field(s): %s", ", ".join(missing_fields)
        )
    rename_map = {k: v for k, v in HGNC_GENE_FIELDS.items() if k in row_fields}
    ht = ht.rename(rename_map)

    if not include_withdrawn:
        ht = ht.filter(ht.status == "Approved")

    # Parse pipe-separated fields into arrays. See _split_pipe_tokens for the
    # quote-strip/split/strip/filter-empty rationale.
    row_fields = get_row_fields(ht)
    for field in HGNC_PIPE_SEPARATED_FIELDS:
        if field in row_fields:
            ht = ht.annotate(**{field: _split_pipe_tokens(ht[field])})

    # lsdb's tokens are (database name, URL) PAIRS, not independent values
    # (e.g. `"Breast Cancer|http://.../LRG_292|http://..."`), so it is kept out
    # of HGNC_PIPE_SEPARATED_FIELDS (see the comment there) and re-paired here
    # into array<struct{name, url}> instead of one flat array<str> that would
    # interleave names and URLs (#391). An odd token count (a trailing name
    # with no URL) gets url=missing rather than erroring.
    if "lsdb" in row_fields:
        tokens = _split_pipe_tokens(ht["lsdb"])
        n_tokens = hl.len(tokens)
        ht = ht.annotate(
            lsdb=hl.range(0, (n_tokens + 1) // 2).map(
                lambda i: hl.struct(
                    name=tokens[2 * i],
                    url=hl.or_missing(2 * i + 1 < n_tokens, tokens[2 * i + 1]),
                )
            )
        )

    if fields is not None:
        ht = ht.select(*fields)

    ht = ht.key_by("hgnc_id")
    return AnnotationTable.from_hail(
        ht, provenance=ctx.provenance(schema_id="hgnc-lookup-v1")
    )
