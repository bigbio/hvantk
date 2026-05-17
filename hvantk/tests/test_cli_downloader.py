import os
from unittest.mock import patch, MagicMock

import pytest
from click.testing import CliRunner
from hvantk.skills.ucsc_cellbrowser.cli import ucsc_downloader
from hvantk.core.constants import UCSC_CELL_BROWSER_BASE_URL


@pytest.fixture
def mock_download_file():
    with patch("hvantk.data.file_utils.download_file") as mock:
        yield mock


@pytest.fixture
def test_output_dir(tmp_path):
    test_dir = tmp_path / "ucsc_data"
    test_dir.mkdir()
    return test_dir


def test_ucsc_downloader_success(mock_download_file, test_output_dir):
    runner = CliRunner()
    result = runner.invoke(
        ucsc_downloader,
        [
            "--dataset",
            "adultPancreas",
            "--output-dir",
            str(test_output_dir),
            "--base_url",
            "http://example.com",
        ],
    )
    assert result.exit_code == 0
    mock_download_file.assert_any_call(
        "http://example.com/test_dataset/expression_matrix.tsv",
        str(test_output_dir),
        "expression_matrix.tsv",
    )
    mock_download_file.assert_any_call(
        "http://example.com/test_dataset/metadata.tsv",
        str(test_output_dir),
        "metadata.tsv",
    )
    assert "Data downloaded to" in result.output


def test_ucsc_downloader_invalid_base_url(mock_download_file, tmp_path):
    runner = CliRunner()
    result = runner.invoke(
        ucsc_downloader,
        [
            "--dataset",
            "test_dataset",
            "--output-dir",
            str(tmp_path),
            "--base_url",
            "invalid_url",
        ],
    )
    assert result.exit_code != 0
    assert "Invalid URL" in result.output


def test_ucsc_downloader_missing_dataset(mock_download_file, tmp_path):
    runner = CliRunner()
    result = runner.invoke(
        ucsc_downloader,
        ["--output-dir", str(tmp_path), "--base_url", "http://example.com"],
    )
    assert result.exit_code != 0
    assert "Error: Missing option '--dataset'" in result.output


def test_ucsc_downloader_output_dir_creation(mock_download_file, tmp_path):
    test_dir = tmp_path / "new_data_dir"
    runner = CliRunner()
    result = runner.invoke(
        ucsc_downloader,
        [
            "--dataset",
            "test_dataset",
            "--output-dir",
            str(test_dir),
            "--base_url",
            "http://example.com",
        ],
    )
    assert result.exit_code == 0
    assert test_dir.exists()
    mock_download_file.assert_any_call(
        "http://example.com/test_dataset/expression_matrix.tsv",
        str(test_dir),
        "expression_matrix.tsv",
    )
    mock_download_file.assert_any_call(
        "http://example.com/test_dataset/metadata.tsv", str(test_dir), "metadata.tsv"
    )


def test_ucsc_downloader_download_failure(mock_download_file, test_output_dir):
    mock_download_file.side_effect = Exception("Download failed")
    runner = CliRunner()
    result = runner.invoke(
        ucsc_downloader,
        [
            "--dataset",
            "test_dataset",
            "--output-dir",
            str(test_output_dir),
            "--base_url",
            "http://example.com",
        ],
    )
    assert result.exit_code != 0
    assert "Download failed" in result.output


@pytest.mark.parametrize(
    "dataset,expect_invalid",
    [
        ("../etc/passwd", True),
        ("bad\\name", True),
        ("bad name", True),
        ("hoc/all-heart", False),  # forward slashes allowed for child dataset paths
    ],
    ids=["traversal", "backslash", "whitespace", "slash-allowed"],
)
def test_ucsc_downloader_dataset_validation(
    mock_download_file, tmp_path, dataset, expect_invalid
):
    """Test dataset name validation (traversal, backslash, whitespace rejected; slash allowed)."""
    runner = CliRunner()
    result = runner.invoke(
        ucsc_downloader,
        [
            "--dataset",
            dataset,
            "--output-dir",
            str(tmp_path),
            "--base_url",
            "http://localhost:9999",
        ],
    )
    if expect_invalid:
        assert result.exit_code != 0
        assert "Invalid dataset value" in result.output
    else:
        assert "Invalid dataset value" not in result.output
