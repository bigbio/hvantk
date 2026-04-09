"""
Module for processing CPTAC (Clinical Proteomic Tumor Analysis Consortium) data.

This module provides functions to convert CPTAC expression and metadata data
into AnnData objects for downstream analysis.
"""

import logging
import re
from typing import Optional, List, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

__all__ = [
    "create_anndata_from_cptac_long",
    "create_anndata_from_cptac_phospho",
]

_SITE_ID_RE = re.compile(r"^(.+)_([A-Z])(\d+)$")


def _parse_site_id(site_id: str) -> Tuple[str, str, int]:
    """Parse a phosphorylation site ID like ``TP53_S315``.

    Returns
    -------
    tuple of (gene_symbol, amino_acid, residue_pos)
        If the pattern does not match, returns ``(site_id, "", 0)``.
    """
    m = _SITE_ID_RE.match(site_id)
    if m:
        return m.group(1), m.group(2), int(m.group(3))
    return site_id, "", 0


def create_anndata_from_cptac_long(
    expression_df: pd.DataFrame,
    metadata_df: Optional[pd.DataFrame] = None,
    gene_id_col: str = "GeneID",
    gene_name_col: str = "Gene Name",
    sample_id_col: str = "SampleID",
    expression_col: str = "Expression",
) -> "ad.AnnData":
    """Create an AnnData from long-format CPTAC expression data.

    Parameters
    ----------
    expression_df : pd.DataFrame
        Long-format DataFrame with columns for gene ID, sample ID, and
        expression value.
    metadata_df : pd.DataFrame, optional
        Sample metadata indexed by *sample_id_col*.
    gene_id_col, gene_name_col, sample_id_col, expression_col : str
        Column names in *expression_df*.

    Returns
    -------
    ad.AnnData
        Expression AnnData (samples x genes) with float32 X.
    """
    import anndata as ad

    logger.info("Creating AnnData from CPTAC long-format expression")

    required_cols = [gene_id_col, sample_id_col, expression_col]
    missing = [c for c in required_cols if c not in expression_df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # Pivot to wide: samples x genes
    wide = expression_df.pivot(
        index=sample_id_col, columns=gene_id_col, values=expression_col
    )
    wide.index.name = sample_id_col

    # Build var (gene metadata)
    var = pd.DataFrame(index=wide.columns)
    var.index.name = gene_id_col
    if gene_name_col and gene_name_col in expression_df.columns:
        name_map = (
            expression_df.drop_duplicates(subset=[gene_id_col])
            .set_index(gene_id_col)[gene_name_col]
        )
        var[gene_name_col] = name_map

    # Build obs
    obs = pd.DataFrame(index=wide.index)
    if metadata_df is not None:
        obs = obs.join(metadata_df, how="left")

    adata = ad.AnnData(
        X=wide.values.astype(np.float32),
        obs=obs,
        var=var,
    )

    logger.info(
        "Created AnnData with %d samples and %d genes", adata.n_obs, adata.n_vars
    )
    return adata


def create_anndata_from_cptac_phospho(
    expression_df: pd.DataFrame,
    metadata_df: Optional[pd.DataFrame] = None,
    site_id_col: str = "SiteID",
    sample_id_col: str = "SampleID",
) -> "ad.AnnData":
    """Create an AnnData from CPTAC phosphoproteomics wide-format data.

    Parameters
    ----------
    expression_df : pd.DataFrame
        Wide-format DataFrame with a *site_id_col* column and one column per
        sample containing intensity values.
    metadata_df : pd.DataFrame, optional
        Sample metadata indexed by *sample_id_col*.
    site_id_col : str
        Column containing site identifiers (e.g. ``TP53_S315``).
    sample_id_col : str
        Name to use for the sample index (used for metadata join).

    Returns
    -------
    ad.AnnData
        Expression AnnData (samples x sites) with parsed site annotations
        in ``var``.
    """
    import anndata as ad

    logger.info("Creating AnnData from CPTAC phospho wide-format data")

    # Set site IDs as index, transpose to samples x sites
    df = expression_df.set_index(site_id_col)
    X = df.T.values.astype(np.float32)

    # Parse site IDs into var annotations
    site_ids = df.index.tolist()
    parsed = [_parse_site_id(sid) for sid in site_ids]
    var = pd.DataFrame(
        {
            "gene_symbol": [p[0] for p in parsed],
            "amino_acid": [p[1] for p in parsed],
            "residue_pos": [p[2] for p in parsed],
        },
        index=pd.Index(site_ids, name=site_id_col),
    )

    # Build obs
    sample_ids = df.columns.tolist()
    obs = pd.DataFrame(index=pd.Index(sample_ids, name=sample_id_col))
    if metadata_df is not None:
        obs = obs.join(metadata_df, how="left")

    adata = ad.AnnData(X=X, obs=obs, var=var)

    logger.info(
        "Created AnnData with %d samples and %d sites", adata.n_obs, adata.n_vars
    )
    return adata
