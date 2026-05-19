"""Hail Table builder for HGNC gene nomenclature data.

This module owns ``create_hgnc_gene_tb``, the canonical builder that turns the
HGNC complete-set TSV into a Hail Table keyed by ``hgnc_id``. It was migrated
out of :mod:`hvantk.core.builders.table` so that everything HGNC-specific
(builder, downloader, drift probe, tests, fixtures, SKILL) lives under the
plugin folder at :mod:`hvantk.skills.hgnc`.

The shared helper ``_create_table_base`` and the field-mapping constants
(``HGNC_GENE_FIELDS``, ``HGNC_PIPE_SEPARATED_FIELDS``) intentionally stay in
their existing modules because they are reused by other builders.
"""

from __future__ import annotations

import logging
from typing import List, Optional

import hail as hl

from hvantk.core.constants import HGNC_GENE_FIELDS, HGNC_PIPE_SEPARATED_FIELDS
from hvantk.core.builders.table import _create_table_base
from hvantk.core.utils.table_utils import get_row_fields

logger = logging.getLogger(__name__)


def create_hgnc_gene_tb(
    input_path: str,
    output_path: str,
    include_withdrawn: bool = False,
    fields: Optional[List[str]] = None,
    overwrite: bool = False,
    export_tsv: bool = False,
) -> "hl.Table":
    """
    Create a Hail Table from HGNC gene nomenclature data keyed by hgnc_id.

    HGNC (HUGO Gene Nomenclature Committee) provides the authoritative source
    for human gene symbols and cross-references to other databases.

    Example usage:
        ht = create_hgnc_gene_tb(
            input_path="/data/hgnc_complete_set.txt",
            output_path="/tables/hgnc.ht"
        )

    Parameters
    ----------
    input_path : str
        Path to the HGNC complete set TSV file (hgnc_complete_set.txt).
    output_path : str
        Path to write the output Hail Table.
    include_withdrawn : bool, optional
        If True, include withdrawn/non-approved genes (default: False).
    fields : list of str, optional
        List of fields to select from the table (default: None, keeps all).
    overwrite : bool, optional
        Whether to overwrite the output file if it exists (default: False).
    export_tsv : bool, optional
        If True, also export a TSV version (default: False).

    Returns
    -------
    hl.Table
        Hail Table keyed by hgnc_id with gene nomenclature and cross-references.

    Notes
    -----
    The table includes:
    - Core identifiers: hgnc_id, gene_symbol, gene_name, status
    - Symbol history: alias_symbols, prev_symbols (as arrays)
    - Cross-references: ensembl_gene_id, entrez_id, uniprot_ids, etc.
    - Gene classification: locus_group, locus_type, gene_group
    - Disease/clinical links: omim_id, orphanet_id, gencc, mane_select

    Pipe-separated fields (alias_symbol, prev_symbol, uniprot_ids, etc.) are
    automatically parsed into arrays.
    """

    def transform(ht: hl.Table) -> hl.Table:
        # Rename fields to standardized names
        logger.info("Renaming HGNC fields to standardized names")
        row_fields = get_row_fields(ht)
        rename_map = {k: v for k, v in HGNC_GENE_FIELDS.items() if k in row_fields}
        ht = ht.rename(rename_map)

        # Filter to approved genes unless include_withdrawn is True
        if not include_withdrawn:
            logger.info("Filtering to approved genes only")
            ht = ht.filter(ht.status == "Approved")

        # Parse pipe-separated fields into arrays
        logger.info("Parsing pipe-separated fields into arrays")
        row_fields = get_row_fields(ht)
        for field in HGNC_PIPE_SEPARATED_FIELDS:
            if field in row_fields:
                # Split on pipe, filter empty strings
                ht = ht.annotate(
                    **{
                        field: hl.if_else(
                            hl.is_defined(ht[field]) & (ht[field] != ""),
                            ht[field].split("\\|").filter(lambda x: x != ""),
                            hl.empty_array(hl.tstr),
                        )
                    }
                )

        # Key by hgnc_id
        logger.info("Keying table by hgnc_id")
        ht = ht.key_by("hgnc_id")

        return ht

    return _create_table_base(
        source_name="HGNC gene nomenclature",
        input_path=input_path,
        output_path=output_path,
        import_func=lambda: hl.import_table(
            paths=input_path,
            impute=False,
            min_partitions=10,
            missing="",
        ),
        transform_func=transform,
        fields=fields,
        overwrite=overwrite,
        export_tsv=export_tsv,
    )
