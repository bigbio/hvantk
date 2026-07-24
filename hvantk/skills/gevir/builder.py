"""Hail Table builder for the GeVIR (Gene Vulnerability and Intolerance Rank) resource.

Owns the Phase B ``build_gevir_metrics`` builder. Imports the GeVIR metrics
TSV keyed by ``gene_id`` and wraps with Provenance.
"""
from __future__ import annotations

import glob
import logging
import os

import hail as hl

logger = logging.getLogger(__name__)


def _resolve_gevir_path(parsed_input) -> str:
    """Resolve the builder's input to the GeVIR table file itself.

    gevir declares no ``lifecycle.parse``, so ``hvantk reprocess`` hands the builder the raw
    *directory* (reprocess forwards ``raw_dir`` verbatim as ``parsed_input`` for no-parse
    datasets), not the file inside it. ``hl.import_table`` cannot read a directory, so resolve
    a directory to the single raw file it contains. A path that already points at a file is
    returned unchanged, which keeps the direct ``build(<file>)`` calls (tests, snapshots)
    working. The GeVIR catalog does not reliably name the file, so glob for the one file
    rather than joining a fixed name (this is why it differs from structure's
    ``_resolve_gtf_path``).
    """
    path = str(parsed_input)
    if not os.path.isdir(path):
        return path
    candidates = sorted(
        f
        for f in glob.glob(os.path.join(path, "*"))
        if os.path.isfile(f) and not os.path.basename(f).startswith(".")
    )
    if len(candidates) != 1:
        raise ValueError(
            f"expected exactly one raw file in directory {path!r}, found "
            f"{len(candidates)}: {[os.path.basename(c) for c in candidates]}"
        )
    return candidates[0]


def build_gevir_metrics(
    parsed_input,
    ctx,
    **params,
):
    """Phase B builder — returns an AnnotationTable.

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
        paths=_resolve_gevir_path(parsed_input),
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
