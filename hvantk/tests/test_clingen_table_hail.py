"""
Hail integration tests for ClinGen Gene-Disease Validity table builder.

These tests require Hail and use the test fixture CSV file.
"""

import pytest
import shutil
from pathlib import Path

from hvantk.tables.table_builders import create_clingen_gene_disease_tb

# Mark as Hail-dependent and slow
pytestmark = [pytest.mark.hail, pytest.mark.slow]

# Test data directory
TEST_DIR = Path(__file__).parent / "testdata"
# Temporary directory for testing
TMP_DIR = Path(__file__).parent / "tmp"


@pytest.fixture(autouse=True)
def setup_teardown():
    """Create and clean up the temporary directory for each test."""
    TMP_DIR.mkdir(exist_ok=True, parents=True)
    yield
    if TMP_DIR.exists():
        shutil.rmtree(TMP_DIR)


def test_create_clingen_gene_disease_tb_default_keying():
    """Test building table with default gene_disease keying."""
    input_path = TEST_DIR / "raw/clingen/clingen_test_sample.csv"
    output_path = TMP_DIR / "clingen_gene_disease.ht"

    tb = create_clingen_gene_disease_tb(
        input_path=str(input_path),
        output_path=str(output_path),
        overwrite=True,
    )

    # Check row count (11 data rows in test fixture)
    assert tb.count() == 11

    # Check key fields
    key_fields = list(tb.key.dtype)
    assert "hgnc_id" in key_fields
    assert "mondo_id" in key_fields

    # Check that _SUCCESS file exists
    success_file = Path(output_path) / "_SUCCESS"
    assert success_file.exists()


def test_create_clingen_gene_disease_tb_gene_keying():
    """Test building table with gene-level aggregation."""
    input_path = TEST_DIR / "raw/clingen/clingen_test_sample.csv"
    output_path = TMP_DIR / "clingen_by_gene.ht"

    tb = create_clingen_gene_disease_tb(
        input_path=str(input_path),
        output_path=str(output_path),
        key_by="gene",
        overwrite=True,
    )

    # Check row count (10 unique genes - BRCA1 appears twice with different diseases)
    assert tb.count() == 10

    # Check key fields
    key_fields = list(tb.key.dtype)
    assert "hgnc_id" in key_fields
    assert "mondo_id" not in key_fields

    # Check aggregated fields exist
    row_fields = list(tb.row.dtype)
    assert "disease_labels" in row_fields
    assert "disease_mondo_pairs" in row_fields
    assert "mondo_ids" in row_fields
    assert "classifications" in row_fields
    assert "n_diseases" in row_fields
    assert "max_classification_label" in row_fields


def test_create_clingen_gene_disease_tb_min_classification():
    """Test filtering by minimum classification level."""
    input_path = TEST_DIR / "raw/clingen/clingen_test_sample.csv"
    output_path = TMP_DIR / "clingen_filtered.ht"

    # Filter to Strong or better (Definitive, Strong)
    tb = create_clingen_gene_disease_tb(
        input_path=str(input_path),
        output_path=str(output_path),
        min_classification="Strong",
        overwrite=True,
    )

    # Should have 9 rows (Definitive: 7, Strong: 2; excludes Moderate: 1, Limited: 1)
    assert tb.count() == 9


def test_create_clingen_gene_disease_tb_definitive_only():
    """Test filtering to Definitive only."""
    input_path = TEST_DIR / "raw/clingen/clingen_test_sample.csv"
    output_path = TMP_DIR / "clingen_definitive.ht"

    tb = create_clingen_gene_disease_tb(
        input_path=str(input_path),
        output_path=str(output_path),
        min_classification="Definitive",
        overwrite=True,
    )

    # Should have 6 rows (all Definitive: BRCA1, BRCA2, TP53, MLH1, MSH2, APC)
    assert tb.count() == 6


def test_create_clingen_gene_disease_tb_hgnc_id_cleaned():
    """Test that HGNC: prefix is stripped from hgnc_id."""
    input_path = TEST_DIR / "raw/clingen/clingen_test_sample.csv"
    output_path = TMP_DIR / "clingen_hgnc_check.ht"

    tb = create_clingen_gene_disease_tb(
        input_path=str(input_path),
        output_path=str(output_path),
        overwrite=True,
    )

    # Collect hgnc_ids and verify none have HGNC: prefix
    hgnc_ids = tb.hgnc_id.collect()
    for hgnc_id in hgnc_ids:
        assert not hgnc_id.startswith("HGNC:"), f"HGNC: prefix not stripped: {hgnc_id}"


def test_create_clingen_gene_disease_tb_export_tsv():
    """Test TSV export functionality."""
    input_path = TEST_DIR / "raw/clingen/clingen_test_sample.csv"
    output_path = TMP_DIR / "clingen_export.ht"

    tb = create_clingen_gene_disease_tb(
        input_path=str(input_path),
        output_path=str(output_path),
        overwrite=True,
        export_tsv=True,
    )

    # Check TSV file was created
    tsv_path = Path(str(output_path) + ".tsv.bgz")
    assert tsv_path.exists()


def test_create_clingen_gene_disease_tb_field_selection():
    """Test field selection."""
    input_path = TEST_DIR / "raw/clingen/clingen_test_sample.csv"
    output_path = TMP_DIR / "clingen_fields.ht"

    tb = create_clingen_gene_disease_tb(
        input_path=str(input_path),
        output_path=str(output_path),
        fields=["gene_symbol", "classification"],
        overwrite=True,
    )

    # Check only selected fields are present (plus key fields)
    row_fields = list(tb.row.dtype)
    # Key fields should still be present
    assert "hgnc_id" in row_fields
    assert "mondo_id" in row_fields
    # Selected fields
    assert "gene_symbol" in row_fields
    assert "classification" in row_fields
    # Non-selected fields should not be present
    assert "disease_label" not in row_fields
    assert "report_url" not in row_fields


def test_create_clingen_gene_disease_tb_gene_aggregation_max_classification():
    """Test that gene aggregation correctly identifies max classification."""
    input_path = TEST_DIR / "raw/clingen/clingen_test_sample.csv"
    output_path = TMP_DIR / "clingen_max_class.ht"

    tb = create_clingen_gene_disease_tb(
        input_path=str(input_path),
        output_path=str(output_path),
        key_by="gene",
        overwrite=True,
    )

    # BRCA1 has Definitive and Strong - max should be Definitive
    brca1_row = tb.filter(tb.gene_symbol == "BRCA1").collect()
    assert len(brca1_row) == 1
    assert brca1_row[0].max_classification_label == "Definitive"
    assert brca1_row[0].n_diseases == 2
