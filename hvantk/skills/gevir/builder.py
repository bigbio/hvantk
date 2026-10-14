"""Hail Table builder for the GeVIR (Gene Vulnerability and Intolerance Rank) resource.

Owns the ``build_gevir_metrics`` builder. Imports the GeVIR metrics
TSV keyed by ``gene_id`` and wraps with Provenance.
"""

from __future__ import annotations

import logging

import hail as hl

from hvantk.core.utils.raw_inputs import resolve_raw_file

logger = logging.getLogger(__name__)

# gevir declares no lifecycle.parse, so `hvantk reprocess` hands the builder the raw
# download directory (see resolve_raw_file). The downloader writes
# gevir_metrics_pmid31873297.tsv.bgz; the committed fixture directory also holds a
# README.md alongside its (differently-named but same-extension) .tsv.bgz fixture, so
# the pattern must pick the data file out and leave the README behind.
_RAW_FILE_PATTERNS = ("*.tsv.bgz",)


def build_gevir_metrics(
    parsed_input,
    ctx,
    **params,
):
    """Plugin builder — returns an AnnotationTable.

    Imports the GeVIR TSV (keyed by gene_id) and wraps it with Provenance.
    Accepts **params for compatibility (fields, etc.).

    Parameters
    ----------
    parsed_input : str | Path
        Path to the GeVIR metrics TSV/BGZ file.
    ctx : hvantk.core.models.BuildContext
        Platform-provided context; supplies provenance.
    **params
        Optional: fields (list of str) to select from the table.
    """
    from hvantk.core.models import AnnotationTable

    fields = params.get("fields", None)

    ht = hl.import_table(
        paths=resolve_raw_file(parsed_input, patterns=_RAW_FILE_PATTERNS),
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
        ht, provenance=ctx.provenance(schema_id="gevir-metrics-v1")
    )
