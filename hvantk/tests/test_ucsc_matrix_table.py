import shutil
import tempfile
from pathlib import Path

import pytest

from hvantk.utils.constants import UCSC_CELL_ID_COLUMN, UCSC_GENE_COLUMN
from hvantk.htables.ucsc import (
    convert_ucsc_metadata_to_hail_table,
    create_mt_from_ucsc_expression_matrix,
)

# Test data directory
TESTDATA_DIR = Path(__file__).parent / "testdata" / "raw" / "ucsc"
METADATA_FILE_PATH = (TESTDATA_DIR / "meta.test.tsv").resolve()
EXPRESSION_MATRIX_FILE_PATH = (TESTDATA_DIR / "exprMatrix.test.tsv.bgz").resolve()


@pytest.fixture
def temp_dir():
    # Create temporary directory for test downloads
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    # Clean up after test
    shutil.rmtree(temp_dir)


def test_convert_ucsc_metadata_to_hail_table(temp_dir):
    """
    Test conversion of UCSC metadata to Hail Table
    """

    # Convert the metadata file to Hail Table
    ht = convert_ucsc_metadata_to_hail_table(
        str(METADATA_FILE_PATH), sep="\t", index_col=0, index_name=UCSC_CELL_ID_COLUMN
    )

    ht.describe()

    # Check if the table has the expected number of rows
    assert ht.count() == 9999


def test_create_mt_from_ucsc_expression_matrix(temp_dir):
    """
    Test creation of Hail MatrixTable from UCSC expression matrix
    """
    # Define output path
    output_path = Path(temp_dir) / "ucsc_expression_matrix.mt"

    # Import metadata file to Hail Table
    metadata_ht = convert_ucsc_metadata_to_hail_table(
        str(METADATA_FILE_PATH), sep="\t", index_col=0, index_name=UCSC_CELL_ID_COLUMN
    )

    # Create the MatrixTable from the expression matrix file
    mt = create_mt_from_ucsc_expression_matrix(
        expression_matrix_path=str(EXPRESSION_MATRIX_FILE_PATH),
        output_path=str(output_path),
        delimiter="\t",
        row_fields=None,
        row_key=UCSC_GENE_COLUMN,
        split_gene_field=True,
        min_partitions=5,
        force_bgz=True,
        overwrite=True,
        metadata_ht=metadata_ht,
    )

    mt.describe()

    # Check if the MatrixTable has the expected number of rows and columns
    assert mt.count_rows() == 199
    assert mt.count_cols() == 9999

    # Check if the _SUCCESS file exists
    success_file = output_path / "_SUCCESS"
    assert success_file.exists()
