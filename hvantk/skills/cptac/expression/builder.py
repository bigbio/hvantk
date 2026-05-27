"""AnnData builder for the CPTAC long-format protein expression resource.

Owns the Phase B ``build_cptac_expression`` builder. Turns a long-format
CPTAC expression TSV plus its sample metadata into an ``ExpressionMatrix``
(samples x genes).
"""

from __future__ import annotations

import logging

logger = logging.getLogger(__name__)


def build_cptac_expression(
    parsed_input,
    ctx,
    *,
    gene_id_col: str = "GeneID",
    gene_name_col: str = "Gene Name",
    sample_id_col: str = "SampleID",
    expression_col: str = "Expression",
):
    """Phase B builder — returns an ExpressionMatrix from long-format CPTAC expression + metadata."""
    import pandas as pd
    from hvantk.core.models import ExpressionMatrix
    from hvantk.core.models.anndata_utils import annotate_column_summary_ad
    from hvantk.skills.cptac.shared.cptac import create_anndata_from_cptac_long

    expression_path = str(parsed_input["expression"])
    metadata_path = str(parsed_input["metadata"])

    expr_df = pd.read_csv(expression_path, sep=None, engine="python")
    meta_df = pd.read_csv(metadata_path, sep=None, engine="python")
    meta_df = meta_df.set_index(sample_id_col)

    adata = create_anndata_from_cptac_long(
        expression_df=expr_df,
        metadata_df=meta_df,
        gene_id_col=gene_id_col,
        gene_name_col=gene_name_col,
        sample_id_col=sample_id_col,
        expression_col=expression_col,
    )

    annotate_column_summary_ad(adata)

    return ExpressionMatrix.from_anndata(
        adata, provenance=ctx.provenance(schema_id="cptac-expression-v1")
    )
