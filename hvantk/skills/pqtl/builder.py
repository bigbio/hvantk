"""Hail Table builder for the pQTL (protein quantitative trait loci) resource.

Phase K plugin promotion — delegates to the legacy create_pqtl_tb function
from hvantk.core.builders.table with the Phase B contract. The legacy
function performs complex multi-step processing (gene-symbol → Ensembl ID
mapping, GTEx variant-ID parsing, SE derivation) that is cleanly encapsulated
in the legacy builder. Using the delegation-stub pattern (as documented in the
Phase K spec) rather than inlining, since the pQTL builder requires external
resources (hgnc_ht) and has complex closure state.

The legacy function stays in place for backward compatibility.
"""
from __future__ import annotations

import logging
import os
import tempfile
from typing import Any

import hail as hl

logger = logging.getLogger(__name__)


def build_pqtl_metrics(
    parsed_input,
    ctx,
    **params,
):
    """Phase B builder — returns an AnnotationTable.

    Delegates to create_pqtl_tb (delegation-stub pattern per Phase K spec)
    because the pQTL builder requires external HGNC table resources and
    complex multi-step transforms that share helpers with other builders.

    Parameters
    ----------
    parsed_input : str | Path
        Path to pQTL allpairs file or directory of per-tissue files.
    ctx : hvantk.core.models.BuildContext
        Platform-provided context; supplies provenance.
    **params
        Optional: reference_genome (str, default "GRCh38"),
                  source (str, default "gtex_fang"),
                  tissue (str), hgnc_ht (str), no_gene_map (bool, default False),
                  p_threshold (float), fields (list of str).

    Notes
    -----
    Phase K delegation stub — the delegation to the legacy builder writes to a
    temporary path and reads the checkpointed table back. This adds an
    extra disk write+read at build time; Phase L cleanup can inline the
    transform directly.
    """
    from hvantk.core.builders.table import create_pqtl_tb
    from hvantk.core.models import AnnotationTable

    # Strip output_path/overwrite/export_tsv — platform owns those
    safe_params = {
        k: v
        for k, v in params.items()
        if k not in ("output_path", "overwrite", "export_tsv")
    }

    with tempfile.TemporaryDirectory() as td:
        tmp_out = os.path.join(td, "pqtl.ht")
        create_pqtl_tb(
            input_path=str(parsed_input),
            output_path=tmp_out,
            overwrite=True,
            **safe_params,
        )
        # Read back before tempdir is cleaned up
        ht = hl.read_table(tmp_out)
        # Force materialisation into a new in-memory representation by
        # collecting the schema — the actual data is lazy until artifact.save()
        # calls ht.write(). We need the table to NOT reference the deleted
        # tempdir. Re-checkpoint to a second temp path that persists until
        # Hail's own temp-file cleanup.
        persistent_tmp = hl.utils.new_temp_file(prefix="pqtl_", extension=".ht")
        ht = ht.checkpoint(persistent_tmp, overwrite=True)

    return AnnotationTable.from_hail(
        ht, provenance=ctx.provenance(schema_id="pqtl-v1")
    )
