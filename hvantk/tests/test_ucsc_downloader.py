import os
import pytest
import tempfile
import shutil
from unittest.mock import patch
from hvantk.data.file_utils import download_file
from hvantk.core.constants import (
    UCSC_CELL_BROWSER_BASE_URL,
    EXPRESSION_MATRIX_FILE_NAME,
    METADATA_FILE_NAME,
)

dataset = "adultPancreas"
url_expression_matrix = (
    f"{UCSC_CELL_BROWSER_BASE_URL}/{dataset}/{EXPRESSION_MATRIX_FILE_NAME}"
)
url_metadata = f"{UCSC_CELL_BROWSER_BASE_URL}/{dataset}/{METADATA_FILE_NAME}"


@pytest.fixture
def temp_dir():
    # Create temporary directory for test downloads
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    # Clean up after test
    shutil.rmtree(temp_dir)


def test_download_expression_matrix(temp_dir):
    """Test download expression matrix file from UCSC Cell Browser"""
    # Mock successful response for unit testing
    with patch("requests.get") as mock_get:
        mock_response = mock_get.return_value
        mock_response.status_code = 200
        mock_response.headers.get.return_value = "1024"  # Content length
        mock_response.iter_content.return_value = [b"test data"]

        download_file(url_expression_matrix, temp_dir, EXPRESSION_MATRIX_FILE_NAME)

        # Verify file exists
        expected_path = os.path.join(temp_dir, EXPRESSION_MATRIX_FILE_NAME)
        assert os.path.exists(expected_path)

        # Verify file contents
        with open(expected_path, "rb") as f:
            assert f.read() == b"test data"


def test_download_metadata(temp_dir):
    """
    Test download metadata file from UCSC Cell Browser

    :return: None
    """
    # Mock successful response for unit testing
    with patch("requests.get") as mock_get:
        mock_response = mock_get.return_value
        mock_response.status_code = 200
        mock_response.headers.get.return_value = "1024"  # Content length
        mock_response.iter_content.return_value = [b"metadata test data"]

        download_file(url_metadata, temp_dir, METADATA_FILE_NAME)
        assert os.path.exists(os.path.join(temp_dir, METADATA_FILE_NAME))
        with open(os.path.join(temp_dir, METADATA_FILE_NAME), "rb") as f:
            assert f.read() == b"metadata test data"
