import os
import shutil
import tempfile
from pathlib import Path

import hail as hl
import matplotlib.pyplot as plt
import pandas as pd
import pytest

from hvantk.core.constants import UCSC_CELL_ID_COLUMN, UCSC_GENE_COLUMN
from hvantk.utils.matrix_utils import (
    summarize_matrix,
    filter_by_metadata,
    filter_by_gene_list,
    filter_by_expression,
    get_top_expressed_genes,
)
from hvantk.visualization import visualize_expression_distribution

# Mark as Hail-dependent and slow
pytestmark = [pytest.mark.hail, pytest.mark.slow]

# Test data directory
TESTDATA_DIR = Path(__file__).parent / "testdata" / "raw" / "ucsc"
METADATA_FILE_PATH = (TESTDATA_DIR / "meta.test.tsv").resolve()
EXPRESSION_MATRIX_FILE_PATH = (TESTDATA_DIR / "exprMatrix.test.tsv.bgz").resolve()


@pytest.fixture
def temp_dir():
    # Create temporary directory for test outputs
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    # Clean up after test
    shutil.rmtree(temp_dir)


@pytest.fixture
def sample_matrix_table(temp_dir):
    """
    Create a sample MatrixTable for testing using the UCSC test data.

    Builds the MT directly via Hail import functions (the old UCSC wrapper
    functions have been removed in favour of AnnData builders).
    """
    output_path = Path(temp_dir) / "test_matrix_utils.mt"

    # Import metadata
    ht = hl.import_table(str(METADATA_FILE_PATH), delimiter="\t", impute=True)
    fields = list(ht.row)
    first_field = fields[0]
    rename_map = {}
    if first_field != UCSC_CELL_ID_COLUMN:
        rename_map[first_field] = UCSC_CELL_ID_COLUMN
    for f in fields:
        if f == first_field:
            continue
        new_name = f.replace(".", "_")
        if new_name != f:
            rename_map[f] = new_name
    if rename_map:
        ht = ht.rename(rename_map)
    metadata_ht = ht.key_by(UCSC_CELL_ID_COLUMN)

    # Import expression matrix
    mt = hl.import_matrix_table(
        str(EXPRESSION_MATRIX_FILE_PATH),
        delimiter="\t",
        row_fields={UCSC_GENE_COLUMN: hl.tstr},
        row_key=UCSC_GENE_COLUMN,
        min_partitions=5,
        force_bgz=True,
    )
    mt = mt.rename({"col_id": UCSC_CELL_ID_COLUMN})

    # Split gene field (pipe-separated)
    from hvantk.utils.expressions import split_field_expr

    mt = mt.key_rows_by()
    mt = mt.annotate_rows(
        **{UCSC_GENE_COLUMN: split_field_expr(mt, field_name=UCSC_GENE_COLUMN)}
    )
    mt = mt.key_rows_by(mt[UCSC_GENE_COLUMN])

    # Annotate with metadata
    mt = mt.annotate_cols(metadata=metadata_ht[mt.col_key])
    mt = mt.checkpoint(str(output_path), overwrite=True)

    return mt


def test_summarize_matrix(sample_matrix_table):
    """
    Test the summarize_matrix function
    """
    # Get summary of the MatrixTable
    summary = summarize_matrix(sample_matrix_table)

    # Verify summary contains expected keys
    assert "dimensions" in summary
    assert "metadata_fields" in summary
    assert "metadata_stats" in summary
    assert "expression_stats" in summary

    # Verify dimensions are correct
    assert summary["dimensions"]["n_samples"] == sample_matrix_table.count_cols()
    assert summary["dimensions"]["n_genes"] == sample_matrix_table.count_rows()

    # Verify metadata fields exist if present in the MatrixTable
    if "metadata" in sample_matrix_table.col:
        assert len(summary["metadata_fields"]) > 0
        assert len(summary["metadata_stats"]) > 0

    # print summary for debugging
    print("Matrix Summary:" + str(summary))


def test_filter_by_metadata(sample_matrix_table):
    """
    Test the filter_by_metadata function
    """
    if "metadata" not in sample_matrix_table.col:
        pytest.skip("Sample MatrixTable does not contain metadata")

    # Get a metadata field to filter on
    metadata_field = list(sample_matrix_table.col.metadata.dtype)[0]

    # Get a value from that field to filter on
    field_values = sample_matrix_table.aggregate_cols(
        hl.agg.collect_as_set(sample_matrix_table.metadata[metadata_field])
    )

    if not field_values:
        pytest.skip(f"No values found for metadata field {metadata_field}")

    filter_value = list(field_values)[0]

    # Filter the MatrixTable
    filtered_mt = filter_by_metadata(
        sample_matrix_table, {metadata_field: filter_value}
    )

    # Check that all samples in filtered table have the filter value
    all_match = filtered_mt.aggregate_cols(
        hl.agg.all(filtered_mt.metadata[metadata_field] == filter_value)
    )

    assert all_match
    assert filtered_mt.count_cols() <= sample_matrix_table.count_cols()


def test_filter_by_gene_list(sample_matrix_table):
    """
    Test the filter_by_gene_list function
    """

    gene_field_name = "gene"  # This is what we expect to be returned

    # Get a few gene IDs to filter on
    gene_ids = sample_matrix_table.aggregate_rows(
        hl.agg.take(sample_matrix_table[gene_field_name], 5)
    )

    # Filter the MatrixTable by gene IDs
    filtered_mt = filter_by_gene_list(
        sample_matrix_table, gene_ids=gene_ids, gene_id_field=gene_field_name
    )

    # Check that filtered MatrixTable has only the specified genes
    assert filtered_mt.count_rows() == len(gene_ids)

    # Check that all gene IDs in the filtered table are in our list
    all_in_list = filtered_mt.aggregate_rows(
        hl.agg.all(hl.literal(set(gene_ids)).contains(filtered_mt[gene_field_name]))
    )

    assert all_in_list


def test_filter_by_expression(sample_matrix_table):
    """
    Test the filter_by_expression function
    """
    # Get expression statistics to set reasonable thresholds
    expr_stats = sample_matrix_table.aggregate_entries(
        hl.struct(
            mean=hl.agg.mean(sample_matrix_table.x),
            std=hl.agg.stats(sample_matrix_table.x).stdev,
            min=hl.agg.min(sample_matrix_table.x),
            max=hl.agg.max(sample_matrix_table.x),
        )
    )

    # Set min threshold to mean - std (but not less than 0)
    min_threshold = max(0, expr_stats.mean - expr_stats.std)

    # Filter by minimum expression
    min_filtered_mt = filter_by_expression(sample_matrix_table, min_expr=min_threshold)

    # Check that all expression values are either 0 or >= min_threshold
    all_valid = min_filtered_mt.aggregate_entries(
        hl.agg.all((min_filtered_mt.x == 0) | (min_filtered_mt.x >= min_threshold))
    )

    assert all_valid

    # Test filter by min samples
    min_samples = 2
    sample_filtered_mt = filter_by_expression(
        sample_matrix_table, min_samples=min_samples
    )

    # Check that the filtered matrix table has fewer or equal number of rows
    assert sample_filtered_mt.count_rows() <= sample_matrix_table.count_rows()


def test_get_top_expressed_genes(sample_matrix_table):
    """
    Test the get_top_expressed_genes function
    """
    # Get top 10 expressed genes
    top_genes_df = get_top_expressed_genes(
        sample_matrix_table, n=10, gene_id_field="gene"
    )

    # Check that the function returns a pandas DataFrame
    assert isinstance(top_genes_df, pd.DataFrame)

    # Check that we got 10 genes
    assert len(top_genes_df) == 10

    # Check required columns are present - using the field name from the returned DataFrame
    gene_field_name = "gene"  # This is what we expect to be returned
    assert gene_field_name in top_genes_df.columns
    assert "mean_expr" in top_genes_df.columns

    # Check that genes are sorted by expression
    assert top_genes_df["mean_expr"].is_monotonic_decreasing

    # Test by_metadata if metadata is available
    if "metadata" in sample_matrix_table.col:
        metadata_field = list(sample_matrix_table.col.metadata.dtype)[0]

        # Get top genes by a metadata field
        top_by_metadata = get_top_expressed_genes(
            sample_matrix_table, n=5, by_metadata=metadata_field, gene_id_field="gene"
        )

        # Check that the result includes the metadata field
        assert metadata_field in top_by_metadata.columns
        # print for debugging
        print("Top expressed genes by metadata:\n", top_by_metadata)


def test_visualize_expression_distribution(
    sample_matrix_table, temp_dir, show_figure=True
):
    """
    Test the visualize_expression_distribution function

    Args:
        sample_matrix_table: Sample MatrixTable fixture
        temp_dir: Temporary directory fixture
        show_figure: If True, will display the figure (default: False)
    """
    # Generate the plot
    fig = visualize_expression_distribution(
        sample_matrix_table, n_bins=20, log_scale=True
    )

    # Check that the figure is a matplotlib figure
    assert isinstance(fig, plt.Figure)

    # Save the figure to check it was created properly
    fig_path = os.path.join(temp_dir, "expr_dist.png")
    fig.savefig(fig_path)

    # Display the figure if requested
    if show_figure:
        plt.figure(fig.number)
        plt.show()

    # Check that the file was created
    assert os.path.exists(fig_path)

    # Test with log_scale=False
    fig2 = visualize_expression_distribution(sample_matrix_table, log_scale=False)
    assert isinstance(fig2, plt.Figure)

    # Display second figure if requested
    if show_figure:
        plt.figure(fig2.number)
        plt.show()
