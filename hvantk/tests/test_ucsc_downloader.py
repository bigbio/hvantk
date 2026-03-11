import os
import pytest
import tempfile
import shutil
from unittest.mock import patch
from click.testing import CliRunner
from hvantk.data.file_utils import download_file
from hvantk.commands.ucsc_downloader import ucsc_downloader
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


# --- Phase 0: Validation and collection warning tests ---


class TestDatasetValidation:
    """Tests for dataset name validation (Phase 0.1)."""

    def test_slashed_child_dataset_passes_validation(self):
        """Forward slashes are allowed for UCSC child dataset paths."""
        runner = CliRunner()
        # Use a non-UCSC base_url to skip URL existence checks;
        # the validation happens before URL construction.
        result = runner.invoke(
            ucsc_downloader,
            ["--dataset", "hoc/all-heart", "--base_url", "http://localhost:9999"],
            catch_exceptions=False,
        )
        # Should NOT see the "Invalid dataset value" error
        assert "Invalid dataset value" not in result.output

    def test_path_traversal_rejected(self):
        """Path traversal attempts with '..' are rejected."""
        runner = CliRunner()
        result = runner.invoke(
            ucsc_downloader,
            ["--dataset", "../etc/passwd"],
        )
        assert "Invalid dataset value" in result.output

    def test_backslash_rejected(self):
        """Backslashes are rejected."""
        runner = CliRunner()
        result = runner.invoke(
            ucsc_downloader,
            ["--dataset", "hoc\\all-heart"],
        )
        assert "Invalid dataset value" in result.output

    def test_whitespace_rejected(self):
        """Dataset names with whitespace are rejected."""
        runner = CliRunner()
        result = runner.invoke(
            ucsc_downloader,
            ["--dataset", "hoc all-heart"],
        )
        assert "Invalid dataset value" in result.output


class TestCollectionWarning:
    """Tests for collection detection and warning (Phase 0.2)."""

    def test_collection_dataset_warns(self):
        """Requesting a known collection prints a warning and exits."""
        runner = CliRunner()
        result = runner.invoke(
            ucsc_downloader,
            ["--dataset", "hoc"],
        )
        assert "collection" in result.output.lower()
        assert "child dataset" in result.output.lower()

    def test_leaf_dataset_no_warning(self):
        """A leaf dataset like 'adultPancreas' does not trigger the collection warning."""
        runner = CliRunner()
        # Use non-UCSC base_url to avoid network calls
        result = runner.invoke(
            ucsc_downloader,
            ["--dataset", "adultPancreas", "--base_url", "http://localhost:9999"],
            catch_exceptions=False,
        )
        assert (
            "collection" not in result.output.lower() or "Warning" not in result.output
        )


class TestListDatasetsSearch:
    """Tests for --list_datasets --search (Phase 2.1)."""

    def test_list_datasets_no_search(self):
        """--list_datasets without search shows all datasets."""
        runner = CliRunner()
        result = runner.invoke(
            ucsc_downloader,
            ["--list_datasets"],
        )
        assert "Available datasets" in result.output
        assert "adultPancreas" in result.output

    def test_list_datasets_with_search(self):
        """--list_datasets --search filters results."""
        runner = CliRunner()
        result = runner.invoke(
            ucsc_downloader,
            ["--list_datasets", "--search", "pancreas"],
        )
        assert "filtered by" in result.output
        assert "pancreas" in result.output.lower()

    def test_list_datasets_search_no_results(self):
        """--search with no matches shows appropriate message."""
        runner = CliRunner()
        result = runner.invoke(
            ucsc_downloader,
            ["--list_datasets", "--search", "zzz_nonexistent_zzz"],
        )
        assert "No datasets matching" in result.output

    def test_list_datasets_shows_collection_info(self):
        """Collections are displayed with dataset count."""
        runner = CliRunner()
        result = runner.invoke(
            ucsc_downloader,
            ["--list_datasets", "--search", "hoc"],
        )
        assert "collection" in result.output.lower()
