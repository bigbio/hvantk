"""Hail Table builder for the gnomAD constraint gene metrics resource.

Owns the ``build_gnomad_metrics_metrics`` builder. Imports the
gnomAD lof_metrics TSV keyed by ``gene_id`` and wraps with Provenance.
"""

from __future__ import annotations

import logging

import hail as hl

from hvantk.core.utils.raw_inputs import resolve_raw_file
from hvantk.core.utils.table_utils import select_fields

logger = logging.getLogger(__name__)

# gnomad-metrics declares no lifecycle.parse, so `hvantk reprocess` hands the builder
# the raw download directory (see resolve_raw_file). The downloader writes a .bgz or
# .tsv depending on version/table; the committed fixtures are also .bgz.
_RAW_FILE_PATTERNS = ("*.bgz", "*.tsv")


def build_gnomad_metrics_metrics(
    parsed_input,
    ctx,
    **params,
):
    """Plugin builder — returns an AnnotationTable.

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
    # Key column: v2.1.1 by_gene has "gene_id"; v4.0 constraint_metrics has no
    # gene_id column (per-transcript rows), so callers pass e.g. key="transcript".
    key = params.get("key", "gene_id")

    ht = hl.import_table(
        paths=resolve_raw_file(parsed_input, patterns=_RAW_FILE_PATTERNS),
        impute=True,
        min_partitions=100,
        key=key,
    )

    # 2. Optional field selection
    if fields is not None:
        logger.info("Selecting fields: %s", fields)
        ht = select_fields(ht, fields)

    # 3. Wrap with provenance
    return AnnotationTable.from_hail(
        ht, provenance=ctx.provenance(schema_id="gnomad-metrics-v1")
    )
