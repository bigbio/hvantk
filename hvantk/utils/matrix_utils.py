"""
Utilities for working with Hail MatrixTable objects, specifically for gene expression data.

This module provides functions for:
- Summarizing MatrixTable contents
- Filtering MatrixTables based on various criteria

Visualization utilities have been moved to hvantk.visualization.hail_expression.

The functions are designed to work with MatrixTables having the following structure:
- Column fields: sample_id (key), metadata (struct with sample attributes)
- Row fields: Gene ID, Gene Name, GeneID (key)
- Entry fields: x (expression values)
"""

from typing import List, Dict, Union, Optional

import hail as hl
import numpy as np
import pandas as pd


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
