"""Hail Table builder for the gnomAD constraint gene metrics resource.

Owns the Phase B ``build_gnomad_metrics_metrics`` builder. Imports the
gnomAD lof_metrics TSV keyed by ``gene_id`` and wraps with Provenance.
"""
from __future__ import annotations

import logging

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
