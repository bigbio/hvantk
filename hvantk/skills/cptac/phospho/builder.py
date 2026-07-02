"""AnnData builder for the CPTAC phosphoproteomics matrix resource.

Owns the Phase B ``build_cptac_phospho`` builder. Turns a wide-format CPTAC
phospho intensity matrix plus its sample metadata into an ``ExpressionMatrix``
(samples x sites).

CPTAC phospho is inherently per-cancer: each cancer type yields its own
matrix+metadata pair and therefore its own AnnData. Like the sibling
``cptac:expression`` dataset, this builder has no ``lifecycle.parse`` stage — the
download already writes the builder's inputs — so ``reprocess`` hands it the raw
download directory and a ``cancer_type`` selects the pair.
"""

from __future__ import annotations

import logging
import os
from collections.abc import Mapping

logger = logging.getLogger(__name__)


def _resolve_phospho_inputs(parsed_input, *, cancer_type):
    """Return ``(expression_csv, metadata_csv)`` paths.

    ``parsed_input`` is either a ``{"expression","metadata"}`` mapping (direct API
    / the per-cancer driver) or a raw download directory (``reprocess``), in which
    case ``cancer_type`` selects the ``cptac-phospho-<ct>-matrix.csv`` /
    ``-metadata.csv`` pair written by :meth:`CPTACPhosphoDataset.download`.
    """
    if isinstance(parsed_input, Mapping):
        return str(parsed_input["expression"]), str(parsed_input["metadata"])

    from hvantk.skills.cptac.shared.constants import CPTAC_CANCER_TYPES

    raw = str(parsed_input)
    if not os.path.isdir(raw):
        raise TypeError(
            "build_cptac_phospho expects a {'expression','metadata'} mapping or a "
            f"raw download directory; got {raw!r}"
        )
    if not cancer_type:
        raise ValueError(
            "cptac:phospho builds one AnnData per cancer type; pass "
            "--plugin-arg cancer_type=<ct> (one of: "
            f"{', '.join(CPTAC_CANCER_TYPES)})"
        )
    expr = os.path.join(raw, f"cptac-phospho-{cancer_type}-matrix.csv")
    meta = os.path.join(raw, f"cptac-phospho-{cancer_type}-metadata.csv")
    missing = [p for p in (expr, meta) if not os.path.exists(p)]
    if missing:
        raise FileNotFoundError(
            f"Missing CPTAC phospho input(s) for cancer_type={cancer_type!r}: "
            f"{missing}. Run `hvantk download cptac-phospho --cancer-type "
            f"{cancer_type}` (or `hvantk reprocess cptac:phospho`) first."
        )
    return expr, meta


def build_cptac_phospho(
    parsed_input,
    ctx,
    *,
    cancer_type: str | None = None,
    site_id_col: str = "Site",
    sample_id_col: str = "SampleID",
    **_ignored,
):
    """Phase B builder — ``ExpressionMatrix`` (samples x sites) for one cancer type.

    ``site_id_col`` defaults to ``"Site"`` to match the matrix CSV header written
    by ``write_matrix_csv`` (``df.index.name = "Site"``). Extra keyword arguments
    that ``reprocess`` forwards from ``--plugin-arg`` (e.g. ``overwrite``, consumed
    only by the download stage) are ignored here.
    """
    import pandas as pd
    from hvantk.core.models import ExpressionMatrix
    from hvantk.core.models.anndata_utils import annotate_column_summary_ad
    from hvantk.skills.cptac.shared.cptac import create_anndata_from_cptac_phospho

    expression_path, metadata_path = _resolve_phospho_inputs(
        parsed_input, cancer_type=cancer_type
    )

    expr_df = pd.read_csv(expression_path, sep=None, engine="python")
    meta_df = pd.read_csv(metadata_path, sep=None, engine="python")
    meta_df = meta_df.set_index(sample_id_col)

    adata = create_anndata_from_cptac_phospho(
        expression_df=expr_df,
        metadata_df=meta_df,
        site_id_col=site_id_col,
        sample_id_col=sample_id_col,
    )

    annotate_column_summary_ad(adata)

    return ExpressionMatrix.from_anndata(
        adata, provenance=ctx.provenance(schema_id="cptac-phospho-v1")
    )
