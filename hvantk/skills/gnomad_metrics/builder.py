"""Hail Table builder for the gnomAD constraint gene metrics resource.

Phase K plugin promotion — inlines the import logic from the legacy
create_gnomad_constraint_gene_metrics_tb function from
hvantk.core.builders.table with the Phase B contract. The legacy function
stays in place for backward compatibility.
"""
from __future__ import annotations

import logging
from typing import Any

import hail as hl

logger = logging.getLogger(__name__)


def build_gnomad_metrics_metrics(
    parsed_input,
    ctx,
    **params,
):
    """Phase B builder — returns an AnnotationTable.

    Imports the gnomAD constraint gene metrics TSV (keyed by gene_id) and
    wraps it with Provenance. Accepts **params for compatibility (fields, etc.).

    Parameters
    ----------
    parsed_input : str | Path
        Path to the gnomAD lof_metrics TSV/BGZ file.
    ctx : hvantk.core.models.BuildContext
        Platform-provided context; supplies provenance.
    **params
        Optional: fields (list of str) to select from the table.
    """
    from hvantk.core.models import AnnotationTable

    fields = params.get("fields", None)

    # 1. Import (inline from create_gnomad_constraint_gene_metrics_tb's import_func)
    ht = hl.import_table(
        paths=str(parsed_input),
        impute=True,
        min_partitions=100,
        key="gene_id",
    )

    # 2. Optional field selection
    if fields is not None:
        logger.info("Selecting fields: %s", fields)
        ht = ht.select(*fields)

    # 3. Wrap with provenance
    return AnnotationTable.from_hail(
        ht, provenance=ctx.provenance(schema_id="gnomad-metrics-v1")
    )
