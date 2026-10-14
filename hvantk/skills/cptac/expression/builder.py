"""AnnData builder for the CPTAC long-format protein expression resource.

Owns the ``build_cptac_expression`` builder. Turns a long-format
CPTAC expression TSV plus its sample metadata into an ``ExpressionMatrix``
(samples x genes).

CPTAC expression is inherently per-cancer: each cancer type yields its own
matrix+metadata pair and therefore its own AnnData. Like the sibling
``cptac:phospho`` dataset, this builder has no ``lifecycle.parse`` stage — the
download already writes the builder's inputs — so ``reprocess`` hands it the raw
download directory and a ``cancer_type`` selects the pair. The directory/mapping
resolution is shared with ``cptac:phospho`` via
``hvantk.skills.cptac.shared.datasets.resolve_raw_dir_inputs``.
"""

from __future__ import annotations

import logging

logger = logging.getLogger(__name__)


def build_cptac_expression(
    parsed_input,
    ctx,
    *,
    cancer_type: str | None = None,
    gene_id_col: str = "GeneID",
    gene_name_col: str = "Gene Name",
    sample_id_col: str = "SampleID",
    expression_col: str = "Expression",
    **_ignored,
):
    """Plugin builder — ``ExpressionMatrix`` (samples x genes) for one cancer type.

    ``parsed_input`` is either a ``{"expression","metadata"}`` mapping (direct API)
    or a raw download directory (the ``reprocess`` path), in which case
    ``cancer_type`` selects the ``cptac-expression-<ct>-matrix.csv`` /
    ``-metadata.csv`` pair written by :meth:`CPTACExpressionDataset.download`. Extra
    keyword arguments that ``reprocess`` forwards from ``--plugin-arg`` (e.g.
    ``overwrite``, consumed only by the download stage) are ignored here.
    """
    import pandas as pd
    from hvantk.core.models import ExpressionMatrix
    from hvantk.core.models.anndata_utils import annotate_column_summary_ad
    from hvantk.skills.cptac.shared.cptac import create_anndata_from_cptac_long
    from hvantk.skills.cptac.shared.datasets import resolve_raw_dir_inputs

    expression_path, metadata_path = resolve_raw_dir_inputs(
        parsed_input, cancer_type=cancer_type, prefix="cptac-expression"
    )

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
