"""
Utilities for working with Hail MatrixTable objects, specifically for gene expression data.

This module provides functions for:
- Summarizing MatrixTable contents
- Column metadata summarization (annotate_column_summary, describe_expression_mt)
- Expression summarization (summarize_expression)
- Filtering MatrixTables based on various criteria

Visualization utilities have been moved to hvantk.visualization.hail_expression.

The functions are designed to work with MatrixTables having the following structure:
- Column fields: sample_id (key), metadata (struct with sample attributes)
- Row fields: Gene ID, Gene Name, GeneID (key)
- Entry fields: x (expression values)
"""

from __future__ import annotations

import logging
from typing import Any, List, Dict, Union, Optional

try:
    import hail as hl
except ImportError:  # allow AnnData-only usage when Hail is not installed
    hl = None  # type: ignore[assignment]

import numpy as np
import pandas as pd

try:
    import anndata as ad
except ImportError:  # AnnData is optional
    ad = None  # type: ignore[assignment]

logger = logging.getLogger(__name__)


def annotate_column_summary(
    mt: hl.MatrixTable,
    metadata_field: str = "metadata",
    max_levels: int = 100,
    top_n_levels: int = 10,
) -> hl.MatrixTable:
    """Annotate MatrixTable globals with a column metadata summary.

    Computes a one-time summary of each column metadata field (type,
    cardinality, value range) and stores it in ``mt.globals.column_summary``.
    Reading the summary back is instant (no Spark job) because globals are
    stored in the MatrixTable header.

    Parameters
    ----------
    mt : hl.MatrixTable
        Input MatrixTable with a ``metadata`` column struct.
    metadata_field : str
        Name of the column struct containing sample/cell metadata.
    max_levels : int
        Categorical fields with at most this many unique values have their
        full level list stored.  Fields exceeding this threshold store only
        the top *top_n_levels* by frequency.
    top_n_levels : int
        Number of most-frequent levels to store for high-cardinality
        categorical fields.

    Returns
    -------
    hl.MatrixTable
        The input MatrixTable with an additional ``column_summary`` global.
    """
    if metadata_field not in mt.col:
        logger.warning(
            "No '%s' field in column schema — skipping column summary",
            metadata_field,
        )
        return mt

    metadata_fields = list(mt.col[metadata_field].dtype)
    if not metadata_fields:
        logger.warning("Metadata struct is empty — skipping column summary")
        return mt

    summaries: Dict[str, Dict[str, Any]] = {}

    for field in metadata_fields:
        field_type = mt.col[metadata_field][field].dtype
        is_num = hl.is_numeric(field_type)

        if is_num:
            stats = mt.aggregate_cols(hl.agg.stats(mt[metadata_field][field]))
            summaries[field] = {
                "dtype": "numeric",
                "n_levels": -1,
                "levels": hl.empty_array(hl.tstr),
                "top_levels": hl.empty_array(hl.tstr),
                "truncated": False,
                "min_val": float(stats.min) if stats.min is not None else 0.0,
                "max_val": float(stats.max) if stats.max is not None else 0.0,
                "mean_val": float(stats.mean) if stats.mean is not None else 0.0,
            }
        else:
            counter = mt.aggregate_cols(hl.agg.counter(mt[metadata_field][field]))
            n_levels = len(counter)

            if n_levels <= max_levels:
                levels = sorted(str(k) for k in counter.keys() if k is not None)
                summaries[field] = {
                    "dtype": "categorical",
                    "n_levels": n_levels,
                    "levels": levels,
                    "top_levels": hl.empty_array(hl.tstr),
                    "truncated": False,
                    "min_val": 0.0,
                    "max_val": 0.0,
                    "mean_val": 0.0,
                }
            else:
                sorted_by_freq = sorted(counter.items(), key=lambda x: -x[1])
                top = [
                    str(k) for k, _v in sorted_by_freq[:top_n_levels] if k is not None
                ]
                summaries[field] = {
                    "dtype": "categorical",
                    "n_levels": n_levels,
                    "levels": hl.empty_array(hl.tstr),
                    "top_levels": top,
                    "truncated": True,
                    "min_val": 0.0,
                    "max_val": 0.0,
                    "mean_val": 0.0,
                }

    # Build a Hail struct for each field, then wrap in a dict global.
    # Using hl.struct per entry ensures uniform schema for the dict values.
    entries = []
    for field_name, info in summaries.items():
        entry = hl.struct(
            dtype=info["dtype"],
            n_levels=hl.int32(info["n_levels"]),
            levels=info["levels"],
            top_levels=info["top_levels"],
            truncated=info["truncated"],
            min_val=hl.float64(info["min_val"]),
            max_val=hl.float64(info["max_val"]),
            mean_val=hl.float64(info["mean_val"]),
        )
        entries.append((field_name, entry))

    summary_dict = hl.dict(entries)
    mt = mt.annotate_globals(column_summary=summary_dict)

    logger.info("Annotated column_summary with %d metadata fields", len(summaries))
    return mt


def describe_expression_mt(
    mt: Union[hl.MatrixTable, str],
    metadata_field: str = "metadata",
) -> Dict[str, Any]:
    """Return a human-readable description of an expression MatrixTable.

    Reads ``column_summary`` from globals (instant, no Spark job).  If the
    global is missing (older MT), falls back to computing it on the fly and
    logs a warning suggesting a rebuild.

    Parameters
    ----------
    mt : hl.MatrixTable or str
        MatrixTable object or path to a checkpointed ``.mt`` on disk.
    metadata_field : str
        Name of the column struct containing sample/cell metadata.

    Returns
    -------
    dict
        Dictionary with keys ``n_genes``, ``n_cols``, ``fields`` (list of
        per-field info dicts).
    """
    if isinstance(mt, str):
        mt = hl.read_matrix_table(mt)

    n_rows, n_cols = mt.count()

    # Check if pre-computed summary exists in globals (instant read)
    has_summary = "column_summary" in mt.globals.dtype

    if has_summary:
        raw_summary = hl.eval(mt.column_summary)
    else:
        logger.warning(
            "column_summary not found in globals — computing on the fly. "
            "Rebuild the MT with annotate_column_summary() for instant access."
        )
        mt_tmp = annotate_column_summary(mt, metadata_field=metadata_field)
        raw_summary = hl.eval(mt_tmp.column_summary)

    fields_info: List[Dict[str, Any]] = []
    lines = [f"Expression MatrixTable: {n_rows:,} genes x {n_cols:,} cells"]
    lines.append("Column metadata fields:")

    for field_name, info in sorted(raw_summary.items()):
        entry: Dict[str, Any] = {"name": field_name, "dtype": info.dtype}
        if info.dtype == "categorical":
            entry["n_levels"] = info.n_levels
            if info.truncated:
                entry["top_levels"] = list(info.top_levels)
                preview = ", ".join(info.top_levels[:5])
                lines.append(
                    f"  {field_name:<18s} categorical   "
                    f"{info.n_levels} levels  [{preview}, ...] (truncated)"
                )
            else:
                entry["levels"] = list(info.levels)
                preview = ", ".join(info.levels[:5])
                suffix = ", ..." if len(info.levels) > 5 else ""
                lines.append(
                    f"  {field_name:<18s} categorical   "
                    f"{info.n_levels} levels  [{preview}{suffix}]"
                )
        else:
            entry["min"] = info.min_val
            entry["max"] = info.max_val
            entry["mean"] = info.mean_val
            lines.append(
                f"  {field_name:<18s} numeric       "
                f"range {info.min_val:.4g}\u2013{info.max_val:.4g}  "
                f"mean {info.mean_val:.4g}"
            )
        fields_info.append(entry)

    description = "\n".join(lines)
    logger.info(description)
    print(description)

    return {"n_genes": n_rows, "n_cols": n_cols, "fields": fields_info}


def _validate_group_by(
    mt: hl.MatrixTable,
    group_by: Union[str, List[str]],
    metadata_field: str,
) -> None:
    """Validate group_by fields using column_summary if available.

    Raises clear errors if fields are missing, numeric, or very
    high-cardinality.
    """
    if isinstance(group_by, str):
        group_by = [group_by]

    meta_dtype = mt.col[metadata_field].dtype
    available = list(meta_dtype)

    # Use column_summary for richer validation when present
    has_summary = "column_summary" in mt.globals.dtype
    summary = hl.eval(mt.column_summary) if has_summary else None

    for field in group_by:
        if field not in meta_dtype:
            raise ValueError(
                f"Field '{field}' not found in {metadata_field}. "
                f"Available fields: {available}"
            )

        if summary and field in summary:
            info = summary[field]
            if info.dtype == "numeric":
                cat_fields = [
                    f
                    for f in available
                    if summary.get(f) and summary[f].dtype == "categorical"
                ]
                raise ValueError(
                    f"Field '{field}' is numeric "
                    f"(range {info.min_val:.4g}–{info.max_val:.4g}), "
                    f"not categorical. "
                    f"Categorical fields available for group_by: "
                    f"{cat_fields}"
                )
            if info.dtype == "categorical" and info.n_levels > 100:
                logger.warning(
                    "Field '%s' has %d groups — this will produce %d "
                    "columns in the summary table.",
                    field,
                    info.n_levels,
                    info.n_levels,
                )
        else:
            # Fallback: check type without column_summary
            if hl.is_numeric(meta_dtype[field]):
                raise ValueError(
                    f"Field '{field}' is numeric, not categorical. "
                    f"Available fields: {available}"
                )


def summarize_expression(
    mt: hl.MatrixTable,
    group_by: Union[str, List[str]],
    filter_by: Optional[Dict[str, Union[str, List[str]]]] = None,
    expr_field: str = "x",
    gene_id_field: str = "GeneID",
    gene_name_field: Optional[str] = "Gene Name",
    min_cells_per_group: int = 50,
    metadata_field: str = "metadata",
    output_path: Optional[str] = None,
    overwrite: bool = False,
) -> hl.Table:
    """Collapse an expression MatrixTable into a gene-level summary Table.

    Groups cells/samples by one or more metadata fields and computes
    per-gene expression statistics (mean, fraction expressed, cell count)
    for each group.

    Parameters
    ----------
    mt : hl.MatrixTable
        Expression MatrixTable (rows = genes, columns = cells/samples).
    group_by : str or list of str
        One or more column metadata fields to group by.  Multiple fields
        are concatenated (e.g., ``["cell_type", "region"]`` → ``"CM_LV"``).
    filter_by : dict, optional
        Pre-filter on column metadata before grouping.  Keys are metadata
        field names; values are a single value or list of values to keep.
    expr_field : str
        Entry field containing expression values (default ``"x"``).
    gene_id_field : str
        Row field for gene IDs (becomes the Table key).
    gene_name_field : str or None
        Row field for gene names.  Set to None to omit.
    min_cells_per_group : int
        Skip groups with fewer cells than this threshold.
    metadata_field : str
        Column struct containing sample/cell metadata.
    output_path : str, optional
        If provided, checkpoint the Table to this path.
    overwrite : bool
        Overwrite existing output if ``output_path`` is given.

    Returns
    -------
    hl.Table
        Table keyed by ``gene_id`` with a ``stats`` dict mapping group
        labels to ``struct{mean, fraction_expressed, n_cells}``.
    """
    if isinstance(group_by, str):
        group_by = [group_by]

    # --- Validate ---
    _validate_group_by(mt, group_by, metadata_field)

    # --- Filter columns ---
    if filter_by:
        mt = filter_by_metadata(mt, filter_by)

    # --- Build group label ---
    if len(group_by) == 1:
        label_expr = hl.str(mt[metadata_field][group_by[0]])
    else:
        label_expr = hl.delimit([hl.str(mt[metadata_field][f]) for f in group_by], "_")
    mt = mt.annotate_cols(_group_label=label_expr)

    # --- Filter groups by min cells ---
    group_counts = mt.aggregate_cols(hl.agg.counter(mt._group_label))
    valid_groups = {g for g, n in group_counts.items() if n >= min_cells_per_group}
    skipped = {g: n for g, n in group_counts.items() if n < min_cells_per_group}
    if skipped:
        logger.warning(
            "Skipping %d group(s) with fewer than %d cells: %s",
            len(skipped),
            min_cells_per_group,
            skipped,
        )
    if not valid_groups:
        raise ValueError(
            f"No groups have >= {min_cells_per_group} cells. "
            f"Group sizes: {group_counts}"
        )

    mt = mt.filter_cols(hl.literal(valid_groups).contains(mt._group_label))

    # --- Group and aggregate ---
    grouped_mt = mt.group_cols_by(mt._group_label).aggregate(
        mean=hl.agg.mean(mt[expr_field]),
        fraction_expressed=hl.agg.fraction(mt[expr_field] > 0),
        n_cells=hl.agg.count(),
    )

    # --- Collect per-gene stats into a dict ---
    grouped_mt = grouped_mt.annotate_rows(
        stats=hl.dict(
            hl.agg.collect(
                hl.tuple(
                    [
                        grouped_mt._group_label,
                        hl.struct(
                            mean=grouped_mt.mean,
                            fraction_expressed=grouped_mt.fraction_expressed,
                            n_cells=hl.int32(grouped_mt.n_cells),
                        ),
                    ]
                )
            )
        )
    )

    # --- Build output Table ---
    row_fields = {gene_id_field: mt.row[gene_id_field]}
    if gene_name_field and gene_name_field in mt.row:
        row_fields[gene_name_field] = mt.row[gene_name_field]

    tb = grouped_mt.rows()
    select_exprs = {"gene_id": tb[gene_id_field], "stats": tb.stats}
    if gene_name_field and gene_name_field in tb.row:
        select_exprs["gene_name"] = tb[gene_name_field]

    tb = tb.select(**select_exprs)
    tb = tb.key_by("gene_id")

    if output_path:
        logger.info("Checkpointing summary table to %s", output_path)
        tb = tb.checkpoint(output_path, overwrite=overwrite)

    logger.info(
        "Summarized %d genes across %d groups (group_by=%s)",
        tb.count(),
        len(valid_groups),
        group_by,
    )
    return tb


def summarize_matrix(mt: hl.MatrixTable) -> Dict:
    """
    Provide a comprehensive summary of a MatrixTable.

    Args:
        mt: Hail MatrixTable to summarize

    Returns:
        Dict containing summary statistics and information about the MatrixTable
    """
    # Basic counts
    n_samples = mt.count_cols()
    n_genes = mt.count_rows()
    n_entries = n_samples * n_genes

    # Sample metadata summaries
    if "metadata" in mt.col:
        metadata_fields = list(mt.col.metadata.dtype)
        metadata_stats = {}

        for field in metadata_fields:
            # Check if the field is likely numerical
            field_type = mt.col.metadata[field].dtype
            is_numeric = hl.is_numeric(field_type)

            if is_numeric:
                # For numerical fields, compute summary statistics
                field_stats = mt.aggregate_cols(hl.agg.stats(mt.metadata[field]))
                # Convert to a more compact dict format
                field_stats = {
                    "mean": field_stats.mean,
                    "std": field_stats.stdev,
                    "min": field_stats.min,
                    "max": field_stats.max,
                    "n": field_stats.n,
                    "sum": field_stats.sum,
                }
            else:
                # For categorical fields, first collect unique values then count them
                unique_values = mt.aggregate_cols(
                    hl.agg.collect_as_set(mt.metadata[field])
                )
                unique_count = len(unique_values)

                if (
                    unique_count <= 50
                ):  # Only show counter for fields with reasonable cardinality
                    field_stats = mt.aggregate_cols(hl.agg.counter(mt.metadata[field]))
                else:
                    # For high cardinality categorical fields, just show counts
                    field_stats = {
                        "n_unique": unique_count,
                        "most_common": mt.aggregate_cols(
                            hl.agg.take(mt.metadata[field], 5)
                        ),
                    }

            metadata_stats[field] = field_stats
    else:
        metadata_fields = []
        metadata_stats = {}

    # Expression value statistics
    if "x" in mt.entry:
        expr_stats = mt.aggregate_entries(
            hl.struct(
                mean=hl.agg.mean(mt.x),
                std=hl.agg.stats(mt.x).stdev,
                min=hl.agg.min(mt.x),
                max=hl.agg.max(mt.x),
                non_zero=hl.agg.count_where(mt.x > 0),
                zeros=hl.agg.count_where(mt.x == 0),
            )
        )
    else:
        expr_stats = {}

    # Prepare the complete summary
    summary = {
        "dimensions": {
            "n_samples": n_samples,
            "n_genes": n_genes,
            "n_entries": n_entries,
            "sparsity": (
                1 - (expr_stats.get("non_zero", 0) / n_entries)
                if n_entries > 0
                else None
            ),
        },
        "metadata_fields": metadata_fields,
        "metadata_stats": metadata_stats,
        "expression_stats": expr_stats,
    }

    return summary


def filter_by_metadata(
    mt: hl.MatrixTable, filters: Dict[str, Union[str, List[str]]]
) -> hl.MatrixTable:
    """
    Filter the MatrixTable by sample metadata.

    Args:
        mt: Hail MatrixTable to filter
        filters: Dictionary mapping metadata field names to values or lists of values
                 to include

    Returns:
        Filtered MatrixTable
    """
    filtered_mt = mt

    for field, values in filters.items():
        if field not in mt.col.metadata.dtype:
            raise ValueError(f"Field {field} not found in metadata")

        # Convert single values to lists for uniform handling
        if not isinstance(values, list):
            values = [values]

        # Create filter expression
        filter_expr = hl.literal(set(values)).contains(mt.metadata[field])
        filtered_mt = filtered_mt.filter_cols(filter_expr)

    return filtered_mt


def filter_by_gene_list(
    mt: hl.MatrixTable,
    gene_ids: Optional[List[str]] = None,
    gene_names: Optional[List[str]] = None,
    gene_id_field: str = "GeneID",
    gene_name_field: str = "Gene Name",
) -> hl.MatrixTable:
    """
    Filter the MatrixTable to include only specified genes.

    Args:
        mt: Hail MatrixTable to filter
        gene_ids: List of gene IDs to include
        gene_names: List of gene names to include
        gene_id_field: Field name for gene IDs (default: 'GeneID')
        gene_name_field: Field name for gene names (default: 'Gene Name')

    Returns:
        Filtered MatrixTable
    """
    filtered_mt = mt

    if gene_ids is not None:
        if gene_id_field not in mt.row:
            raise ValueError(
                f"Field {gene_id_field} not found in row fields. Available fields: {list(mt.row)}"
            )
        filtered_mt = filtered_mt.filter_rows(
            hl.literal(set(gene_ids)).contains(filtered_mt[gene_id_field])
        )

    if gene_names is not None:
        if gene_name_field not in mt.row:
            raise ValueError(
                f"Field {gene_name_field} not found in row fields. Available fields: {list(mt.row)}"
            )
        filtered_mt = filtered_mt.filter_rows(
            hl.literal(set(gene_names)).contains(filtered_mt[gene_name_field])
        )

    return filtered_mt


def filter_by_expression(
    mt: hl.MatrixTable,
    min_expr: float = None,
    max_expr: float = None,
    min_samples: int = None,
) -> hl.MatrixTable:
    """
    Filter the MatrixTable based on expression values.

    Args:
        mt: Hail MatrixTable to filter
        min_expr: Minimum expression value to include
        max_expr: Maximum expression value to include
        min_samples: Minimum number of samples a gene must be expressed in

    Returns:
        Filtered MatrixTable
    """
    filtered_mt = mt

    if min_expr is not None:
        # For entry-level filtering, where we want to keep the structure
        # but set values below threshold to NA or 0
        filtered_mt = filtered_mt.annotate_entries(
            x_filtered=hl.if_else(filtered_mt.x >= min_expr, filtered_mt.x, 0)
        )
        # Replace the original x with the filtered one
        filtered_mt = filtered_mt.drop("x")
        filtered_mt = filtered_mt.rename({"x_filtered": "x"})

    if max_expr is not None:
        # Apply a maximum threshold
        filtered_mt = filtered_mt.annotate_entries(
            x_filtered=hl.if_else(filtered_mt.x <= max_expr, filtered_mt.x, max_expr)
        )
        filtered_mt = filtered_mt.drop("x")
        filtered_mt = filtered_mt.rename({"x_filtered": "x"})

    if min_samples is not None:
        # Count samples where gene is expressed per row
        # We need to annotate the rows with the count first, then filter
        filtered_mt = filtered_mt.annotate_rows(
            n_expressed_samples=hl.agg.count_where(filtered_mt.x > 0)
        )

        # Filter genes expressed in at least min_samples
        filtered_mt = filtered_mt.filter_rows(
            filtered_mt.n_expressed_samples >= min_samples
        )

        # Clean up the annotation if desired
        filtered_mt = filtered_mt.drop("n_expressed_samples")

    return filtered_mt


def get_top_expressed_genes(
    mt: hl.MatrixTable,
    n: int = 20,
    by_metadata: str = None,
    gene_id_field: str = "GeneID",
    gene_name_field: str = "Gene Name",
) -> pd.DataFrame:
    """
    Get the top expressed genes overall or by a metadata category.

    Args:
        mt: Hail MatrixTable to analyze
        n: Number of top genes to return
        by_metadata: If provided, return top genes for each category in this metadata field
        gene_id_field: Field name for gene IDs (default: 'GeneID')
        gene_name_field: Field name for gene names (default: 'Gene Name')

    Returns:
        Pandas DataFrame with top genes
    """
    # Verify field names exist in the row fields
    row_fields = list(mt.row)
    if gene_id_field not in row_fields:
        raise ValueError(
            f"Field {gene_id_field} not found in row fields. Available fields: {row_fields}"
        )

    # Check if gene name field is available
    has_gene_names = gene_name_field in row_fields

    if by_metadata is None:
        # Calculate mean expression per gene
        mt = mt.annotate_rows(mean_expr=hl.agg.mean(mt.x))

        # Select fields to include in output
        tb = mt.rows().key_by()
        fields_to_select = [gene_id_field, "mean_expr"]
        if has_gene_names:
            fields_to_select.insert(1, gene_name_field)
        df = tb.select(*fields_to_select).to_pandas()

        df = df.sort_values("mean_expr", ascending=False).head(n)

    else:
        # Group by metadata field and get top genes for each group
        if by_metadata not in mt.col.metadata.dtype:
            raise ValueError(f"Field {by_metadata} not found in metadata")

        # First get the categories
        categories = mt.aggregate_cols(hl.agg.collect_as_set(mt.metadata[by_metadata]))

        # Initialize result dataframe
        result_dfs = []

        # For each category, filter and get top genes
        for category in categories:
            mt_filtered = filter_by_metadata(mt, {by_metadata: category})
            mt_filtered = mt_filtered.annotate_rows(
                mean_expr=hl.agg.mean(mt_filtered.x)
            )

            # Get the top n genes
            tb = mt_filtered.rows().key_by()
            fields_to_select = [gene_id_field, "mean_expr"]
            if has_gene_names:
                fields_to_select.insert(1, gene_name_field)
            df_category = tb.select(*fields_to_select).to_pandas()

            df_category = df_category.sort_values("mean_expr", ascending=False).head(n)
            df_category[by_metadata] = category

            result_dfs.append(df_category)

        df = pd.concat(result_dfs)

    return df


# ---------------------------------------------------------------------------
# AnnData-based expression analysis functions
# ---------------------------------------------------------------------------


def describe_expression_ad(adata: "ad.AnnData") -> Dict[str, Any]:
    """Return a summary dict describing an AnnData expression object.

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data matrix.

    Returns
    -------
    dict
        Dictionary with keys ``n_obs``, ``n_vars``, ``fields`` (list of
        per-field info dicts with at least ``name`` and ``dtype``).
    """
    fields_info: List[Dict[str, Any]] = []

    if "column_summary" in adata.uns:
        raw_summary = adata.uns["column_summary"]
        for field_name, info in sorted(raw_summary.items()):
            entry: Dict[str, Any] = {"name": field_name, "dtype": info["dtype"]}
            if info["dtype"] == "categorical":
                entry["n_unique"] = info.get("n_unique")
            else:
                entry["min"] = info.get("min")
                entry["max"] = info.get("max")
            fields_info.append(entry)
    else:
        for col in adata.obs.columns:
            series = adata.obs[col]
            if pd.api.types.is_numeric_dtype(series):
                fields_info.append(
                    {
                        "name": col,
                        "dtype": "numeric",
                        "min": float(np.nanmin(series.values)),
                        "max": float(np.nanmax(series.values)),
                    }
                )
            else:
                fields_info.append(
                    {
                        "name": col,
                        "dtype": "categorical",
                        "n_unique": int(series.nunique()),
                    }
                )

    return {"n_obs": adata.n_obs, "n_vars": adata.n_vars, "fields": fields_info}


def filter_by_metadata_ad(
    adata: "ad.AnnData",
    filters: Dict[str, Union[str, List[str]]],
) -> "ad.AnnData":
    """Filter an AnnData object by obs metadata.

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data matrix.
    filters : dict
        Mapping of obs column name to a value or list of values to keep.

    Returns
    -------
    anndata.AnnData
        Filtered copy of the input.
    """
    mask = np.ones(adata.n_obs, dtype=bool)
    for col, values in filters.items():
        if not isinstance(values, list):
            values = [values]
        mask &= adata.obs[col].isin(values).values
    return adata[mask].copy()


def summarize_expression_ad(
    adata: "ad.AnnData",
    group_by: Union[str, List[str]],
    filter_by: Optional[Dict[str, Union[str, List[str]]]] = None,
    min_cells_per_group: int = 1,
) -> pd.DataFrame:
    """Collapse an AnnData expression matrix into a gene-level summary.

    Parameters
    ----------
    adata : anndata.AnnData
        Expression AnnData (obs = cells, var = genes).
    group_by : str or list of str
        One or more obs columns to group by.  Multiple columns are
        concatenated with ``"_"`` to form a composite label.
    filter_by : dict, optional
        Pre-filter on obs metadata before grouping (passed to
        :func:`filter_by_metadata_ad`).
    min_cells_per_group : int
        Discard groups with fewer cells than this threshold.

    Returns
    -------
    pd.DataFrame
        Long-form DataFrame with columns ``gene_id``, ``group``,
        ``mean``, ``fraction_expressed``, ``n_cells``.
    """
    if filter_by:
        adata = filter_by_metadata_ad(adata, filter_by)

    if isinstance(group_by, str):
        group_by = [group_by]

    # Build composite group labels
    if len(group_by) == 1:
        labels = adata.obs[group_by[0]].astype(str)
    else:
        labels = adata.obs[group_by[0]].astype(str)
        for col in group_by[1:]:
            labels = labels + "_" + adata.obs[col].astype(str)

    adata.obs["_group_label"] = labels.values

    gene_ids = adata.var_names.tolist()
    X = adata.X

    records: List[Dict[str, Any]] = []
    for group_name, idx in adata.obs.groupby("_group_label").groups.items():
        n_cells = len(idx)
        if n_cells < min_cells_per_group:
            continue
        positions = [adata.obs.index.get_loc(i) for i in idx]
        sub = X[positions, :]
        means = np.asarray(sub.mean(axis=0)).ravel()
        frac_expr = np.asarray((sub > 0).mean(axis=0)).ravel()
        for j, gid in enumerate(gene_ids):
            records.append(
                {
                    "gene_id": gid,
                    "group": group_name,
                    "mean": float(means[j]),
                    "fraction_expressed": float(frac_expr[j]),
                    "n_cells": n_cells,
                }
            )

    return pd.DataFrame(records)
