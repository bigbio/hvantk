"""AnnData builder for the CPTAC phosphoproteomics matrix resource.

This module owns ``build_cptac_phospho_ad``, the canonical builder that turns
a wide-format CPTAC phospho intensity matrix plus its sample metadata into an
``anndata.AnnData`` object (samples x sites). It was migrated out of
:mod:`hvantk.tables.matrix_builders` so that everything CPTAC-specific
(builders, helpers, downloader, dataset class, tests, fixtures, SKILLs)
lives under the plugin folder at :mod:`hvantk.skills.cptac`.
"""

from __future__ import annotations

import logging
from typing import Optional

import anndata as ad

logger = logging.getLogger(__name__)


def build_cptac_phospho_ad(
    expression_path: str,
    metadata_path: str,
    output_path: Optional[str] = None,
    site_id_col: str = "SiteID",
    sample_id_col: str = "SampleID",
    overwrite: bool = False,
) -> "ad.AnnData":
    """Build an AnnData from CPTAC phosphoproteomics wide-format matrix + metadata.

    Parameters
    ----------
    expression_path : str
        Path to wide-format phospho expression TSV/CSV (sites x samples).
    metadata_path : str
        Path to sample metadata TSV/CSV.
    output_path : str, optional
        If provided, save the AnnData as ``.h5ad``.
    site_id_col : str
        Column containing site identifiers (e.g. ``TP53_S315``).
    sample_id_col : str
        Column containing sample identifiers in metadata.
    overwrite : bool
        Allow overwriting *output_path* if it exists.

    Returns
    -------
    ad.AnnData
        Expression AnnData (samples x sites) with parsed site annotations
        in ``var`` and metadata in ``obs``.
    """
    import pandas as pd

    from hvantk.skills.cptac.shared.cptac import create_anndata_from_cptac_phospho
    from hvantk.core.models.anndata_utils import (
        build_anndata_metadata,
        annotate_column_summary_ad,
        save_anndata,
    )

    logger.info("Reading CPTAC phospho expression from %s", expression_path)
    expr_df = pd.read_csv(expression_path, sep=None, engine="python")

    logger.info("Reading CPTAC phospho metadata from %s", metadata_path)
    meta_df = pd.read_csv(metadata_path, sep=None, engine="python")
    meta_df = meta_df.set_index(sample_id_col)

    logger.info("Building CPTAC phospho AnnData")
    adata = create_anndata_from_cptac_phospho(
        expression_df=expr_df,
        metadata_df=meta_df,
        site_id_col=site_id_col,
        sample_id_col=sample_id_col,
    )

    adata.uns["hvantk_metadata"] = build_anndata_metadata("CPTAC", expression_path)
    annotate_column_summary_ad(adata)

    if output_path:
        save_anndata(adata, output_path, overwrite=overwrite)

    return adata


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
