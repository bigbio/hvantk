"""Hail Table builder for HGNC gene nomenclature data.

Owns the Phase B ``build_hgnc_gene_lookup`` builder. Turns the HGNC
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


def build_hgnc_gene_lookup(
    parsed_input,
    ctx,
    *,
    include_withdrawn: bool = False,
    fields=None,
):
    """Phase B builder — returns an AnnotationTable.

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

    Returns
    -------
    hvantk.core.models.AnnotationTable
        The HGNC table wrapped with Provenance.
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
    rename_map = {k: v for k, v in HGNC_GENE_FIELDS.items() if k in row_fields}
    ht = ht.rename(rename_map)

    if not include_withdrawn:
        ht = ht.filter(ht.status == "Approved")

    # Parse pipe-separated fields into arrays.
    #
    # hgnc_complete_set.txt QUOTES its multi-value fields (prev_symbol is written
    # `"H1F4|HIST1H1E"`), so the quotes must come off BEFORE the split -- otherwise the
    # first and last elements keep a stray `"` and every alias/prev-symbol lookup for the
    # clean symbol misses. No HGNC value legitimately contains a double quote, so
    # stripping them all is safe and simpler than anchoring to the ends.
    row_fields = get_row_fields(ht)
    for field in HGNC_PIPE_SEPARATED_FIELDS:
        if field in row_fields:
            ht = ht.annotate(
                **{
                    field: hl.if_else(
                        hl.is_defined(ht[field]) & (ht[field] != ""),
                        ht[field]
                        .replace('"', "")
                        .split("\\|")
                        .map(lambda x: x.strip())
                        .filter(lambda x: x != ""),
                        hl.empty_array(hl.tstr),
                    )
                }
            )

    if fields is not None:
        ht = ht.select(*fields)

    ht = ht.key_by("hgnc_id")
    return AnnotationTable.from_hail(
        ht, provenance=ctx.provenance(schema_id="hgnc-lookup-v1")
    )
