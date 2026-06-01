"""Hail Table builder for MSigDB gene-set GMT files.

Owns the Phase B ``build_msigdb_genesets`` builder. Turns an MSigDB GMT file
(e.g., C2 Canonical Pathways) into an ``AnnotationTable`` keyed by
``set_name``.

The GMT format is tab-separated with variable-width rows: each row is one
gene set, column 1 is the set name, column 2 is a description (a
gsea-msigdb URL in MSigDB-issued files), and columns 3..N are gene members.
Imported via ``hl.import_lines`` because ``hl.import_table`` rejects
variable column counts.
"""

from __future__ import annotations

import logging

import hail as hl

logger = logging.getLogger(__name__)


def build_msigdb_genesets(
    parsed_input,
    ctx,
    **params,
):
    """Phase B builder — returns an AnnotationTable.

    Each row is one MSigDB gene set keyed by ``set_name``. The ``genes`` column
    is ``array<str>``. msigdb may be promoted to a GeneSet collection
    artifact in a future phase; for Phase B it stays as an AnnotationTable.
    """
    from hvantk.core.models import AnnotationTable

    ht = hl.import_lines(paths=str(parsed_input), min_partitions=4)
    parts = ht.text.split("\t")
    ht = ht.select(
        set_name=parts[0],
        source_url=parts[1],
        genes=parts[2:],
    )
    ht = ht.filter(ht.set_name != "")
    ht = ht.key_by("set_name")

    return AnnotationTable.from_hail(
        ht, provenance=ctx.provenance(schema_id="msigdb-genesets-v1")
    )
