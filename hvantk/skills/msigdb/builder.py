"""Hail Table builder for MSigDB gene-set GMT files.

Owns the ``build_msigdb_genesets`` builder. Turns an MSigDB GMT file
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

from hvantk.core.utils.raw_inputs import resolve_raw_file

logger = logging.getLogger(__name__)

# msigdb is bring-your-own (credentialed): a user downloads and unzips a collection's
# GMT by hand, so `hvantk reprocess` hands this builder the raw directory it was
# placed in (see resolve_raw_file). Any unzipped MSigDB GMT matches; the committed
# fixture is also a .gmt.
_RAW_FILE_PATTERNS = ("*.gmt",)


def build_msigdb_genesets(
    parsed_input,
    ctx,
    **params,
):
    """Plugin builder — returns an AnnotationTable.

    Each row is one MSigDB gene set keyed by ``set_name``. The ``genes`` column
    is ``array<str>``. msigdb may be promoted to a GeneSet collection
    artifact later; for now it stays an AnnotationTable.
    """
    from hvantk.core.models import AnnotationTable

    ht = hl.import_lines(
        paths=resolve_raw_file(parsed_input, patterns=_RAW_FILE_PATTERNS),
        min_partitions=4,
    )
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
