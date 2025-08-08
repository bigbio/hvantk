"""
Tests for the CPTAC data processing module.
"""

import pytest
import pandas as pd
import hail as hl
import tempfile
import os
from hvantk.htables.cptac import (
    convert_cptac_expression_to_matrix_table,
    convert_cptac_metadata_to_table,
    create_cptac_matrix_table,
    save_cptac_matrix_table
)

# Initialize Hail
hl.init()

@pytest.fixture
def sample_expression_data():
    """Create sample expression data for testing."""
    return pd.DataFrame({
        'GeneID': ['GENE1', 'GENE2', 'GENE3'],
        'Gene Name': ['Gene One', 'Gene Two', 'Gene Three'],
        'SampleID': ['SAMPLE1', 'SAMPLE1', 'SAMPLE1'],
        'Expression': [1.0, 2.0, 3.0]
    })

@pytest.fixture
def sample_metadata_data():
    """Create sample metadata data for testing."""
    return pd.DataFrame({
        'SampleID': ['SAMPLE1'],
        'TumorType': ['TypeA'],
        'Stage': ['I'],
        'Age': [50],
        'TumorSize': [2.5]
    })

def test_convert_cptac_expression_to_matrix_table(sample_expression_data):
    """Test conversion of expression data to MatrixTable."""
    mt = convert_cptac_expression_to_matrix_table(
        sample_expression_data,
        gene_id_col='GeneID',
        gene_name_col='Gene Name',
        sample_id_col='SampleID',
        expression_col='Expression'
    )
    
    # Check basic properties
    assert mt.count_rows() == 3  # 3 genes
    assert mt.count_cols() == 1  # 1 sample
    
    # Check row keys
    assert mt.row_key.dtype == hl.tstruct(GeneID=hl.tstr)
    
    # Check column keys
    assert mt.col_key.dtype == hl.tstruct(SampleID=hl.tstr)
    
    # Check annotations
    assert 'gene_name' in mt.row
    assert mt.row.gene_name.dtype == hl.tstr
    
    # Check entry field
    assert 'Expression' in mt.entry
    assert mt.entry.Expression.dtype == hl.tfloat64

def test_convert_cptac_expression_to_matrix_table_missing_columns():
    """Test handling of missing required columns."""
    df = pd.DataFrame({
        'GeneID': ['GENE1'],
        'SampleID': ['SAMPLE1']
        # Missing Expression column
    })
    
    with pytest.raises(ValueError, match="Missing required columns"):
        convert_cptac_expression_to_matrix_table(df)

def test_convert_cptac_metadata_to_table(sample_metadata_data):
    """Test conversion of metadata to Table."""
    ht = convert_cptac_metadata_to_table(
        sample_metadata_data,
        sample_id_col='SampleID',
        categorical_cols=['TumorType', 'Stage'],
        numeric_cols=['Age', 'TumorSize']
    )
    
    # Check basic properties
    assert ht.count() == 1  # 1 sample
    
    # Check key
    assert ht.key.dtype == hl.tstruct(SampleID=hl.tstr)
    
    # Check column types
    assert ht.TumorType.dtype == hl.tstr
    assert ht.Stage.dtype == hl.tstr
    assert ht.Age.dtype == hl.tfloat64
    assert ht.TumorSize.dtype == hl.tfloat64

def test_convert_cptac_metadata_to_table_missing_sample_id():
    """Test handling of missing sample ID column."""
    df = pd.DataFrame({
        'TumorType': ['TypeA']
        # Missing SampleID column
    })
    
    with pytest.raises(ValueError, match="Sample ID column"):
        convert_cptac_metadata_to_table(df)

def test_create_cptac_matrix_table(sample_expression_data, sample_metadata_data):
    """Test creation of complete MatrixTable with expression and metadata."""
    mt = create_cptac_matrix_table(
        sample_expression_data,
        sample_metadata_data,
        gene_id_col='GeneID',
        gene_name_col='Gene Name',
        sample_id_col='SampleID',
        expression_col='Expression',
        categorical_cols=['TumorType', 'Stage'],
        numeric_cols=['Age', 'TumorSize']
    )
    
    # Check basic properties
    assert mt.count_rows() == 3  # 3 genes
    assert mt.count_cols() == 1  # 1 sample
    
    # Check metadata annotations
    assert 'TumorType' in mt.col
    assert 'Stage' in mt.col
    assert 'Age' in mt.col
    assert 'TumorSize' in mt.col
    
    # Check metadata types
    assert mt.col.TumorType.dtype == hl.tstr
    assert mt.col.Stage.dtype == hl.tstr
    assert mt.col.Age.dtype == hl.tfloat64
    assert mt.col.TumorSize.dtype == hl.tfloat64

def test_create_cptac_matrix_table_sample_mismatch():
    """Test handling of sample ID mismatches between expression and metadata."""
    expression_df = pd.DataFrame({
        'GeneID': ['GENE1'],
        'SampleID': ['SAMPLE1'],
        'Expression': [1.0]
    })
    
    metadata_df = pd.DataFrame({
        'SampleID': ['SAMPLE2'],  # Different sample ID
        'TumorType': ['TypeA']
    })
    
    # Should raise ValueError due to sample mismatch
    with pytest.raises(ValueError, match="Sample ID mismatches found"):
        create_cptac_matrix_table(expression_df, metadata_df)
