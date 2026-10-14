"""AnnData builder for the CPTAC phosphoproteomics matrix resource.

Owns the ``build_cptac_phospho`` builder. Turns a wide-format CPTAC
phospho intensity matrix plus its sample metadata into an ``ExpressionMatrix``
(samples x sites).

CPTAC phospho is inherently per-cancer: each cancer type yields its own
matrix+metadata pair and therefore its own AnnData. Like the sibling
``cptac:expression`` dataset, this builder has no ``lifecycle.parse`` stage — the
download already writes the builder's inputs — so ``reprocess`` hands it the raw
download directory and a ``cancer_type`` selects the pair. The directory/mapping
resolution is shared with ``cptac:expression`` via
``hvantk.skills.cptac.shared.datasets.resolve_raw_dir_inputs``.
"""

from __future__ import annotations

import logging

logger = logging.getLogger(__name__)


def build_cptac_phospho(
    parsed_input,
    ctx,
    *,
    cancer_type: str | None = None,
    site_id_col: str = "Site",
    sample_id_col: str = "SampleID",
    **_ignored,
):
    """Plugin builder — ``ExpressionMatrix`` (samples x sites) for one cancer type.

    ``site_id_col`` defaults to ``"Site"`` to match the matrix CSV header written
    by ``write_matrix_csv`` (``df.index.name = "Site"``). Extra keyword arguments
    that ``reprocess`` forwards from ``--plugin-arg`` (e.g. ``overwrite``, consumed
    only by the download stage) are ignored here.
    """
    import pandas as pd
    from hvantk.core.models import ExpressionMatrix
    from hvantk.core.models.anndata_utils import annotate_column_summary_ad
    from hvantk.skills.cptac.shared.cptac import create_anndata_from_cptac_phospho
    from hvantk.skills.cptac.shared.datasets import resolve_raw_dir_inputs

    expression_path, metadata_path = resolve_raw_dir_inputs(
        parsed_input, cancer_type=cancer_type, prefix="cptac-phospho"
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
