"""AnnData builder for the EBI Expression Atlas resource.

Owns the Phase B ``build_expression_atlas`` builder. Turns an Expression
Atlas baseline bulk-RNA-seq expression TSV plus its SDRF metadata file into
an ``ExpressionMatrix`` (samples x genes).

The shared ``annotate_column_summary_ad`` helper lives in
``hvantk/core/models/anndata_utils`` and is reused by every anndata
builder (UCSC, CPTAC, ...).
"""

from __future__ import annotations

import logging

logger = logging.getLogger(__name__)


def build_expression_atlas(
    parsed_input,
    ctx,
    *,
    gene_column: str = "Gene ID",
    gene_name_column: str = "Gene Name",
    delimiter: str = "\t",
):
    """Phase B builder — returns an ExpressionMatrix.

    Parameters
    ----------
    parsed_input : dict
        Must contain keys ``expression_matrix`` (path to gene x sample TSV)
        and ``sdrf`` (path to SDRF metadata file).
    ctx : hvantk.core.models.BuildContext
        Platform-provided context; supplies provenance.

    Returns
    -------
    hvantk.core.models.ExpressionMatrix
        AnnData-backed ExpressionMatrix wrapped with Provenance.
    """
    from hvantk.core.models import ExpressionMatrix
    from hvantk.core.models.anndata_utils import annotate_column_summary_ad
    from hvantk.skills.expression_atlas.shared.expression_atlas import (
        convert_sdrf_to_dataframe,
        create_anndata_from_expression_atlas,
    )

    expression_matrix_path = str(parsed_input["expression_matrix"])
    sdrf_file = str(parsed_input["sdrf"])

    metadata_df = convert_sdrf_to_dataframe(sdrf_file)
    adata = create_anndata_from_expression_atlas(
        expression_matrix_path=expression_matrix_path,
        metadata_df=metadata_df,
        gene_id_column=gene_column,
        gene_name_column=gene_name_column,
        delimiter=delimiter,
    )

    annotate_column_summary_ad(adata)

    return ExpressionMatrix.from_anndata(
        adata, provenance=ctx.provenance(schema_id="expression-atlas-dataset-v1")
    )
