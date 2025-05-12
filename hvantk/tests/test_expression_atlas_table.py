import shutil
import tempfile
from pathlib import Path

import pytest


from hvantk.htables.expression_atlas import (
    convert_sdrf_to_hail_table,
    create_mt_from_expression_atlas_matrix,
)

# Test data directory
TESTDATA_DIR = Path(__file__).parent / "testdata" / "raw" / "expression_atlas"
SDRF_FILE_PATH = (TESTDATA_DIR / "E-MTAB-6798.condensed-sdrf.tsv").resolve()
EXPRESSION_MATRIX_FILE_PATH = (TESTDATA_DIR / "E-MTAB-6798-transcripts-tpms.tsv.bgz").resolve()


@pytest.fixture
def temp_dir():
    # Create temporary directory for test downloads
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    # Clean up after test
    shutil.rmtree(temp_dir)


def test_convert_sdrf_to_hail_table(temp_dir):
    """
    Test conversion of SDRF metadata to Hail Table
    """
    output_path = Path(temp_dir) / "expression_atlas_table.ht"
    # Convert the SDRF file to Hail Table
    ht = convert_sdrf_to_hail_table(
        str(SDRF_FILE_PATH),
        output_file=str(output_path)
    )

    ht.describe()

    # Check if the table has the expected number of rows
    assert ht.count() == 317

def test_create_mt_from_expression_atlas_matrix(temp_dir):
    """
    Test creation of Hail MatrixTable from Expression Atlas expression matrix
    """
    # Define output path
    output_path = Path(temp_dir) / "expression_atlas_matrix.mt"

    # Import SDRF file to Hail Table
    metadata_ht = convert_sdrf_to_hail_table(
        str(SDRF_FILE_PATH),
    )

    # Create the MatrixTable from the expression matrix file
    mt = create_mt_from_expression_atlas_matrix(
        expression_matrix_path=str(EXPRESSION_MATRIX_FILE_PATH),
        metadata_ht=metadata_ht
    )

    mt.describe()

    # Check if the MatrixTable has the expected number of rows and columns
    assert mt.count_rows() == 116643
    assert mt.count_cols() == 317

