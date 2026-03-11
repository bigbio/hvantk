"""
Hail adapter for Wilcoxon rank-sum marker gene detection.

Thin bridge between Hail MatrixTables and the pure numpy/scipy
implementation in :mod:`hvantk.utils.wilcoxon`.  Handles data extraction
from a MatrixTable (Phase 1 pre-filter + dense matrix collection) and
delegates all statistical computation to the core module.
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional, Set, Tuple, Union

import hail as hl
import numpy as np
import pandas as pd

from hvantk.utils.wilcoxon import (
    WilcoxonParams,
    rank_genes_groups,
    results_to_gene_set_collection,
)
from hvantk.utils.gene_sets import GeneSetCollection

logger = logging.getLogger(__name__)

__all__ = [
    "extract_expression_for_wilcoxon",
    "wilcoxon_markers_from_mt",
]


def extract_expression_for_wilcoxon(
    mt: hl.MatrixTable,
    group_by: List[str],
    expr_field: str = "x",
    gene_id_field: str = "GeneID",
    gene_name_field: Optional[str] = "Gene Name",
    metadata_field: str = "metadata",
    min_cells_per_group: int = 3,
    filter_by: Optional[Dict[str, Union[str, List[str]]]] = None,
    candidate_gene_ids: Optional[Set[str]] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Optional[np.ndarray]]:
    """Extract dense expression matrix from a Hail MatrixTable.

    Parameters
    ----------
    mt : hl.MatrixTable
        Expression MatrixTable (rows=genes, columns=cells).
    group_by : list of str
        Metadata field(s) used to build group labels.
    expr_field : str
        Entry field containing expression values.
    gene_id_field : str
        Row field for gene IDs.
    gene_name_field : str or None
        Row field for gene names.
    metadata_field : str
        Column struct containing cell metadata.
    min_cells_per_group : int
        Skip groups with fewer cells.
    filter_by : dict, optional
        Pre-filter columns by metadata values.
    candidate_gene_ids : set of str, optional
        If provided, restrict to these gene IDs only.

    Returns
    -------
    expression : np.ndarray
        Dense matrix ``(n_cells, n_genes)``.
    group_labels : np.ndarray
        Group label per cell ``(n_cells,)``.
    gene_ids : np.ndarray
        Gene IDs ``(n_genes,)``.
    gene_names : np.ndarray or None
        Gene names ``(n_genes,)`` if available.
    """
    # --- Apply column filters ---
    if filter_by:
        from hvantk.utils.matrix_utils import filter_by_metadata
        mt = filter_by_metadata(mt, filter_by)

    # --- Build group label ---
    if len(group_by) == 1:
        label_expr = hl.str(mt[metadata_field][group_by[0]])
    else:
        label_expr = hl.delimit(
            [hl.str(mt[metadata_field][f]) for f in group_by], "_"
        )
    mt = mt.annotate_cols(_group_label=label_expr)

    # --- Filter groups by min cells ---
    group_counts = mt.aggregate_cols(hl.agg.counter(mt._group_label))
    valid_groups = {
        g for g, n in group_counts.items() if n >= min_cells_per_group
    }
    skipped = {
        g: n for g, n in group_counts.items() if n < min_cells_per_group
    }
    if skipped:
        logger.warning(
            "Skipping %d group(s) with fewer than %d cells: %s",
            len(skipped), min_cells_per_group, skipped,
        )
    if not valid_groups:
        raise ValueError(
            f"No groups have >= {min_cells_per_group} cells. "
            f"Group sizes: {group_counts}"
        )
    mt = mt.filter_cols(hl.literal(valid_groups).contains(mt._group_label))

    # --- Filter rows to candidate genes ---
    if candidate_gene_ids is not None:
        mt = mt.filter_rows(
            hl.literal(candidate_gene_ids).contains(mt[gene_id_field])
        )

    # --- Collect group labels ---
    col_data = mt.cols().select("_group_label")
    labels_list = col_data.aggregate(hl.agg.collect(col_data._group_label))
    group_labels = np.array(labels_list, dtype=str)

    # --- Collect gene IDs and names ---
    row_fields = [gene_id_field]
    has_names = gene_name_field and gene_name_field in mt.row
    if has_names:
        row_fields.append(gene_name_field)

    row_data = mt.rows().select(*row_fields)
    rows_pd = row_data.to_pandas()
    gene_ids = rows_pd[gene_id_field].values.astype(str)
    gene_names = rows_pd[gene_name_field].values.astype(str) if has_names else None

    n_genes = len(gene_ids)
    n_cells = len(group_labels)
    logger.info(
        "Extracting dense matrix: %d genes x %d cells", n_genes, n_cells
    )

    # --- Collect expression values per gene row ---
    # Aggregate per-row: collect all expression values across cells
    expr_per_gene = mt.annotate_rows(
        _expr_values=hl.agg.collect(mt[expr_field])
    )
    expr_lists = expr_per_gene.rows().select("_expr_values")
    expr_pd = expr_lists.to_pandas()

    # Build dense matrix (n_cells x n_genes)
    expression = np.column_stack(
        [np.array(row, dtype=np.float64) for row in expr_pd["_expr_values"]]
    )

    logger.info(
        "Extracted expression matrix: shape %s, %d groups",
        expression.shape, len(valid_groups),
    )

    return expression, group_labels, gene_ids, gene_names


def wilcoxon_markers_from_mt(
    mt: hl.MatrixTable,
    group_by: Union[str, List[str]],
    filter_by: Optional[Dict[str, Union[str, List[str]]]] = None,
    summary: Optional[Union[hl.Table, str]] = None,
    params: Optional[WilcoxonParams] = None,
    expr_field: str = "x",
    gene_id_field: str = "GeneID",
    gene_name_field: Optional[str] = "Gene Name",
    metadata_field: str = "metadata",
    min_cells_per_group: int = 3,
) -> Tuple[pd.DataFrame, GeneSetCollection]:
    """Full Wilcoxon marker detection pipeline from a MatrixTable.

    Two-phase approach:

    1. **Pre-filter** (Hail): If a summary table is provided, use it to
       identify candidate genes by fold-change and fraction expressed.
       Otherwise, pass all genes (limited by ``max_candidates``).
    2. **Wilcoxon** (numpy/scipy): Extract dense matrix for candidates,
       run vectorised rank-sum tests, correct p-values.

    Parameters
    ----------
    mt : hl.MatrixTable
        Expression MatrixTable.
    group_by : str or list of str
        Metadata field(s) for grouping.
    filter_by : dict, optional
        Pre-filter columns by metadata values.
    summary : hl.Table, str, or None
        Pre-computed summary table from ``summarize_expression``.
        If provided, used for candidate pre-filtering.
    params : WilcoxonParams, optional
        Test parameters.
    expr_field, gene_id_field, gene_name_field, metadata_field : str
        Field names matching the MatrixTable schema.
    min_cells_per_group : int
        Minimum cells per group.

    Returns
    -------
    results_df : pd.DataFrame
        Full Wilcoxon results table.
    collection : GeneSetCollection
        Top significant markers per group.
    """
    if params is None:
        params = WilcoxonParams()

    if isinstance(group_by, str):
        group_by = [group_by]

    # --- Phase 1: Identify candidates via summary (optional) ---
    candidate_gene_ids = None
    if summary is not None:
        candidate_gene_ids = _candidates_from_summary(
            summary, params, gene_id_field, gene_name_field,
        )
        logger.info(
            "Pre-filter from summary: %d candidate genes",
            len(candidate_gene_ids),
        )

    # --- Extract dense matrix ---
    expression, group_labels, gene_ids, gene_names = extract_expression_for_wilcoxon(
        mt,
        group_by=group_by,
        expr_field=expr_field,
        gene_id_field=gene_id_field,
        gene_name_field=gene_name_field,
        metadata_field=metadata_field,
        min_cells_per_group=min_cells_per_group,
        filter_by=filter_by,
        candidate_gene_ids=candidate_gene_ids,
    )

    # --- Phase 2: Wilcoxon rank-sum ---
    results_df = rank_genes_groups(
        expression, group_labels, gene_ids, gene_names, params,
    )

    # --- Convert to GeneSetCollection ---
    background_genes = set(gene_ids)
    gene_col = "gene_name" if gene_names is not None else "gene_id"
    collection = results_to_gene_set_collection(
        results_df,
        background_genes=background_genes,
        top_n=params.top_n,
        alpha=params.alpha,
        gene_col=gene_col,
    )

    return results_df, collection


def _candidates_from_summary(
    summary,
    params: WilcoxonParams,
    gene_id_field: str = "gene_id",
    gene_name_field: Optional[str] = "gene_name",
) -> Set[str]:
    """Extract candidate gene IDs from a summary table using pre-filters."""
    from hvantk.utils.gene_sets import _summary_to_dataframe

    df = _summary_to_dataframe(summary, gene_id_field, gene_name_field)

    # Find group columns
    stat_cols = [c for c in df.columns if c.endswith("_mean")]
    groups = [c.rsplit("_mean", 1)[0] for c in stat_cols]

    if not groups:
        logger.warning("No groups found in summary — returning all genes.")
        return set(df[gene_id_field].tolist())

    candidates = set()
    for group in groups:
        mean_col = f"{group}_mean"
        frac_col = f"{group}_fraction_expressed"
        if mean_col not in df.columns or frac_col not in df.columns:
            continue

        other_cols = [f"{g}_mean" for g in groups if g != group]
        if other_cols:
            other_mean = df[other_cols].mean(axis=1)
        else:
            other_mean = df[mean_col]

        denom = other_mean.replace(0, 1e-10)
        fc = df[mean_col] / denom

        mask = (fc >= params.min_fold_change) & (
            df[frac_col] >= params.min_fraction_expressed
        )
        group_candidates = set(df.loc[mask, gene_id_field].tolist())

        # Limit by max_candidates per group (top by FC)
        if len(group_candidates) > params.max_candidates:
            top_idx = fc[mask].nlargest(params.max_candidates).index
            group_candidates = set(df.loc[top_idx, gene_id_field].tolist())

        candidates.update(group_candidates)

    return candidates
