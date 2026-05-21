"""Hail Table builder for AlphaGenome variant effect predictions.

Phase K plugin promotion — delegates to the legacy create_alphagenome_tb
function from hvantk.core.builders.table with the Phase B contract. The
legacy function involves an external API streamer (AlphaGenomeStreamer) with
complex setup/teardown logic that is cleanly encapsulated in the legacy
builder. Using the delegation-stub pattern (per Phase K spec).

The legacy function stays in place for backward compatibility.
"""
from __future__ import annotations

import logging
import os
import tempfile
from typing import Any

import hail as hl

logger = logging.getLogger(__name__)


def build_alphagenome_predictions(
    parsed_input,
    ctx,
    **params,
):
    """Phase B builder — returns an AnnotationTable.

    Delegates to create_alphagenome_tb (delegation-stub pattern per Phase K
    spec) because the AlphaGenome builder wraps an external streaming API
    (AlphaGenomeStreamer) with complex setup/teardown that is cleanly
    encapsulated in the legacy builder.

    Parameters
    ----------
    parsed_input : str | Path
        Path to a Hail Table (.ht) or TSV with chrom/pos/ref/alt columns.
    ctx : hvantk.core.models.BuildContext
        Platform-provided context; supplies provenance.
    **params
        Required: config_path (str) — path to AlphaGenome YAML config.
        Optional: no_resume (bool, default False).

    Notes
    -----
    Phase K delegation stub — the legacy builder writes to a temporary
    output directory; predictions.json and checkpoint files are produced
    there. The resulting Hail Table is then re-checkpointed to a persistent
    Hail temp file before the tempdir is cleaned up.
    Phase L cleanup should inline the AlphaGenomeStreamer call directly.
    """
    from hvantk.core.builders.table import create_alphagenome_tb
    from hvantk.core.models import AnnotationTable

    # Strip output_path/overwrite — platform owns those
    safe_params = {
        k: v
        for k, v in params.items()
        if k not in ("output_path", "overwrite")
    }

    with tempfile.TemporaryDirectory() as td:
        tmp_out = os.path.join(td, "alphagenome_out")
        os.makedirs(tmp_out, exist_ok=True)
        create_alphagenome_tb(
            input_path=str(parsed_input),
            output_path=tmp_out,
            overwrite=True,
            **safe_params,
        )
        # The legacy builder writes the Hail Table to
        # <output_path>/alphagenome_variants.ht
        ht_path = os.path.join(tmp_out, "alphagenome_variants.ht")
        ht = hl.read_table(ht_path)
        # Re-checkpoint to a persistent Hail temp path before the tempdir
        # is cleaned up, so the lazy Table doesn't reference a deleted path.
        persistent_tmp = hl.utils.new_temp_file(
            prefix="alphagenome_", extension=".ht"
        )
        ht = ht.checkpoint(persistent_tmp, overwrite=True)

    return AnnotationTable.from_hail(
        ht, provenance=ctx.provenance(schema_id="alphagenome-v1")
    )
