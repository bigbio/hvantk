"""Hail Table builder for AlphaGenome variant effect predictions.

Owns the full pipeline: drives the local AlphaGenomePipeline to call the
external API for each variant, persists predictions and checkpoint files
to a temp directory, then re-emits a Hail Table keyed by (locus, alleles)
under the Phase B contract.
"""
from __future__ import annotations

import logging
import os
import shutil
import tempfile
from typing import Any

import hail as hl

from hvantk.core.utils.file_utils import resolve_compression

logger = logging.getLogger(__name__)


def _run_alphagenome_pipeline(
    input_path: str,
    output_path: str,
    config_path: str,
    no_resume: bool = False,
    overwrite: bool = False,
) -> "hl.Table":
    """Run AlphaGenome predictions and return a checkpointed Hail Table.

    Drives AlphaGenomePipeline for the API calls, then builds a minimal
    Hail Table keyed by (locus, alleles) from the input variants and
    checkpoints it under output_path/alphagenome_variants.ht.
    """
    from hvantk.skills.alphagenome.pipelines import AlphaGenomePipeline

    if overwrite and os.path.isdir(output_path):
        shutil.rmtree(output_path)

    pipeline = AlphaGenomePipeline(
        input_path=input_path,
        output_dir=output_path,
        config_path=config_path,
        no_resume=no_resume,
    )
    pipeline.setup()
    try:
        for _batch in pipeline.stream():
            pass  # checkpointing handled internally
    finally:
        pipeline.teardown()

    logger.info("Creating AlphaGenome variants table from %s", input_path)

    if input_path.endswith(".ht"):
        ht = hl.read_table(input_path)
        if "locus" not in ht.row or "alleles" not in ht.row:
            raise ValueError(
                "Input Hail Table must contain 'locus' and 'alleles' fields for "
                "AlphaGenome table builder output."
            )
    else:
        resolved_path, force_bgz = resolve_compression(input_path)
        import_kwargs: dict[str, Any] = {"impute": True}
        if force_bgz:
            import_kwargs["force_bgz"] = True
        ht = hl.import_table(resolved_path, **import_kwargs)
        required = {"chrom", "pos", "ref", "alt"}
        missing = required.difference(set(ht.row))
        if missing:
            raise ValueError(
                "TSV input missing required columns for AlphaGenome: "
                f"{', '.join(sorted(missing))}"
            )
        ht = ht.annotate(
            locus=hl.locus(ht.chrom, hl.int(ht.pos), reference_genome="GRCh38"),
            alleles=[ht.ref, ht.alt],
        ).key_by("locus", "alleles")

    if os.path.isdir(output_path):
        table_output_path = os.path.join(output_path, "alphagenome_variants.ht")
    else:
        table_output_path = output_path
    logger.info("Checkpointing AlphaGenome variants table to %s", table_output_path)
    return ht.checkpoint(output=table_output_path, overwrite=overwrite)


def build_alphagenome_predictions(parsed_input, ctx, **params):
    """Phase B builder — returns an AnnotationTable.

    Parameters
    ----------
    parsed_input : str | Path
        Path to a Hail Table (.ht) or TSV with chrom/pos/ref/alt columns.
    ctx : hvantk.core.models.BuildContext
        Platform-provided context; supplies provenance.
    **params
        Required: config_path (str) — path to AlphaGenome YAML config.
        Optional: no_resume (bool, default False).
    """
    from hvantk.core.models import AnnotationTable

    safe_params = {
        k: v
        for k, v in params.items()
        if k not in ("output_path", "overwrite")
    }

    with tempfile.TemporaryDirectory() as td:
        tmp_out = os.path.join(td, "alphagenome_out")
        os.makedirs(tmp_out, exist_ok=True)
        _run_alphagenome_pipeline(
            input_path=str(parsed_input),
            output_path=tmp_out,
            overwrite=True,
            **safe_params,
        )
        ht_path = os.path.join(tmp_out, "alphagenome_variants.ht")
        ht = hl.read_table(ht_path)
        persistent_tmp = hl.utils.new_temp_file(
            prefix="alphagenome_", extension=".ht"
        )
        ht = ht.checkpoint(persistent_tmp, overwrite=True)

    return AnnotationTable.from_hail(
        ht, provenance=ctx.provenance(schema_id="alphagenome-v1")
    )
