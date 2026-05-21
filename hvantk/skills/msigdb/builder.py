"""Hail Table builder for MSigDB gene-set GMT files.

This module owns ``create_msigdb_tb``, the canonical builder that turns an
MSigDB GMT file (e.g., C2 Canonical Pathways) into a Hail Table keyed by
``set_name``. It was migrated out of :mod:`hvantk.core.builders.table` so
that everything MSigDB-specific (builder, drift probe, tests, fixtures,
SKILL) lives under the plugin folder at :mod:`hvantk.skills.msigdb`.

The shared helper ``create_table_base`` intentionally stays in
``hvantk.core.builders.table`` because it is reused by every builder.
"""

from __future__ import annotations

import logging

import hail as hl

from hvantk.core.builders.table import create_table_base

logger = logging.getLogger(__name__)


def create_msigdb_tb(
    input_path: str,
    output_path: str,
    overwrite: bool = False,
    export_tsv: bool = False,
) -> "hl.Table":
    """Create a Hail Table from an MSigDB GMT gene-set file keyed by set_name.

    Implements the contract in ``hvantk/skills/msigdb/SKILL.md``. The GMT
    format is tab-separated with variable-width rows: each row is one gene
    set, column 1 is the set name, column 2 is a description (a gsea-msigdb
    URL in MSigDB-issued files), and columns 3..N are gene members.

    Imported via ``hl.import_lines`` (one row per line, single ``text``
    field) because ``hl.import_table`` rejects variable column counts.
    The transform splits the line on ``\\t`` and slices the gene members
    into an ``array<str>``.

    Parameters
    ----------
    input_path : str
        Path to the GMT file (e.g., ``c2.cp.v2026.1.Hs.symbols.gmt``).
    output_path : str
        Path to write the output Hail Table (``.ht`` directory).
    overwrite : bool, optional
        Overwrite the output if present (default: False).
    export_tsv : bool, optional
        If True, also export a TSV alongside the HT (default: False).

    Returns
    -------
    hl.Table
        Hail Table keyed by ``set_name`` with fields:
          - ``set_name: str``
          - ``source_url: str`` (the GMT description column, verbatim)
          - ``genes: array<str>`` (gene members; symbols for ``.Hs.symbols.gmt``)
    """

    def transform(ht: hl.Table) -> hl.Table:
        # hl.import_lines yields rows with `file: str` and `text: str`.
        # Split on tab; slice [2:] for the variable-width gene-member tail.
        parts = ht.text.split("\t")
        ht = ht.select(
            set_name=parts[0],
            source_url=parts[1],
            genes=parts[2:],
        )
        # Defensive: drop blank lines (parts would be a 1-element array).
        ht = ht.filter(ht.set_name != "")
        ht = ht.key_by("set_name")
        return ht

    return create_table_base(
        source_name="MSigDB gene sets",
        input_path=input_path,
        output_path=output_path,
        import_func=lambda: hl.import_lines(paths=input_path, min_partitions=4),
        transform_func=transform,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )


def build_msigdb_genesets(
    parsed_input,
    ctx,
):
    """Phase B builder — returns an AnnotationTable.

    Each row is one MSigDB gene set keyed by set_name. The genes member
    column is an array<str>. msigdb may be promoted to a GeneSet collection
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
