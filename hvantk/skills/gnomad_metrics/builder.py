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
    from pathlib import Path

    from hvantk.core.models import AnnotationTable

    fields = params.get("fields", None)
    # Key column: v2.1.1 by_gene has "gene_id"; v4.0 constraint_metrics has no
    # gene_id column (per-transcript rows), so callers pass e.g. key="transcript".
    key = params.get("key", "gene_id")

    # `hvantk reprocess` hands a download-only plugin's builder the raw_dir (no
    # parse stage); resolve the constraint file inside it. An explicit file path
    # (run_builder_for_spec / tests) is used as-is.
    src = Path(parsed_input)
    if src.is_dir():
        candidates = sorted(src.glob("*.bgz")) + sorted(src.glob("*.tsv"))
        if not candidates:
            raise FileNotFoundError(
                f"No gnomAD constraint file (*.bgz/*.tsv) found in {src}"
            )
        if len(candidates) > 1:
            # Do not silently pick the first: a raw dir holding more than one constraint
            # file (e.g. v2.1.1 + v4.0, or by_gene + by_transcript) would build the wrong
            # table. Mirror the sibling gevir builder and fail loud. (PR #222 review.)
            raise ValueError(
                f"expected exactly one gnomAD constraint file in {src}, found "
                f"{len(candidates)}: {[c.name for c in candidates]}. Point --raw-dir at a "
                f"directory holding a single version's constraint file, or pass an explicit "
                f"file path."
            )
        src = candidates[0]
        logger.info("Resolved gnomAD constraint file: %s", src)

    ht = hl.import_table(
        paths=str(src),
        impute=True,
        min_partitions=100,
        key=key,
    )

    # 2. Optional field selection
    if fields is not None:
        logger.info("Selecting fields: %s", fields)
        ht = ht.select(*fields)

    # 3. Wrap with provenance
    return AnnotationTable.from_hail(
        ht, provenance=ctx.provenance(schema_id="gnomad-metrics-v1")
    )
