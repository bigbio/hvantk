"""AnnData builder for the CPTAC phosphoproteomics matrix resource.

Owns the Phase B ``build_cptac_phospho`` builder. Turns a wide-format CPTAC
phospho intensity matrix plus its sample metadata into an ``ExpressionMatrix``
(samples x sites).
"""

from __future__ import annotations

import logging

logger = logging.getLogger(__name__)


def build_cptac_phospho(
    parsed_input,
    ctx,
    *,
    site_id_col: str = "SiteID",
    sample_id_col: str = "SampleID",
):
    """Phase B builder — returns an ExpressionMatrix from CPTAC wide-format phospho + metadata."""
    import pandas as pd
    from hvantk.core.models import ExpressionMatrix
    from hvantk.core.models.anndata_utils import (
        annotate_column_summary_ad,
        build_anndata_metadata,
    )
    from hvantk.skills.cptac.shared.cptac import create_anndata_from_cptac_phospho

    expression_path = str(parsed_input["expression"])
    metadata_path = str(parsed_input["metadata"])

    expr_df = pd.read_csv(expression_path, sep=None, engine="python")
    meta_df = pd.read_csv(metadata_path, sep=None, engine="python")
    meta_df = meta_df.set_index(sample_id_col)

    adata = create_anndata_from_cptac_phospho(
        expression_df=expr_df,
        metadata_df=meta_df,
        site_id_col=site_id_col,
        sample_id_col=sample_id_col,
    )

    adata.uns["hvantk_metadata"] = build_anndata_metadata("CPTAC-phospho", expression_path)
    annotate_column_summary_ad(adata)

    return ExpressionMatrix.from_anndata(
        adata, provenance=ctx.provenance(schema_id="cptac-phospho-v1")
    )
