import os
from pathlib import Path

import pytest
from click.testing import CliRunner
from hvantk.commands import expression_atlas_downloader

# Define the test directory path
TEST_DIR = Path(__file__).parent.parent


@pytest.fixture
def download_path(tmpdir):
    return str(tmpdir.mkdir("downloads"))


@pytest.mark.skipif(
    os.environ.get("CI", "false").lower() == "true",
    reason="Test skipped in CI environment - downloads all files and runs too slow",
)
def test_download_experiments_config(download_path):
    """
    This test download all files from experiments in the config file and
    should run only locally and not in the CI/CD pipeline (too slow).

    :param download_path: Path where downloaded files will be stored
    :return: None
    """
    # Create a CliRunner instance
    runner = CliRunner()

    # Call the download_experiments function with the config file and download path
    config_path = str(TEST_DIR / "resources/expression_atlas.json")
    result = runner.invoke(
        expression_atlas_downloader.download_experiments,
        ["--config_path", config_path, "--download_path", download_path],
    )

    # Assert that the command was successful
    assert result.exit_code == 0

    # Assert that the file was created
    file_path = os.path.join(download_path, "E-GTEX-8.condensed-sdrf.tsv")
    assert os.path.exists(file_path)

    # Assert that the file content is not empty
    assert os.path.getsize(file_path) > 0


def test_download_experiments_accession_with_config(download_path):
    # Create a CliRunner instance
    runner = CliRunner()

    # Call the download_experiments function with the accession and download path
    config_path = str(TEST_DIR / "resources/expression_atlas.json")
    accession = "E-MTAB-6798"
    result = runner.invoke(
        expression_atlas_downloader.download_experiments,
        [
            "--accession",
            accession,
            "--config_path",
            config_path,
            "--download_path",
            download_path,
        ],
    )

    # Assert that the command was successful
    assert result.exit_code == 0

    # Assert that the sdrf file was created
    sdrf_file_path = os.path.join(download_path, f"{accession}.condensed-sdrf.tsv")
    assert os.path.exists(sdrf_file_path)

    # Assert that the tpm-transcript file was created
    tpm_file_path = os.path.join(download_path, f"{accession}-transcripts-tpms.tsv")
    assert os.path.exists(tpm_file_path)

    # Assert that the file content is not empty
    assert os.path.getsize(sdrf_file_path) > 0


def test_download_experiments_without_accession_or_config(download_path):
    # Create a CliRunner instance
    runner = CliRunner()

    # Call the download_experiments function without the accession and download path
    result = runner.invoke(
        expression_atlas_downloader.download_experiments,
        ["--download_path", download_path],
    )

    # Assert that the command was not successful
    assert result.exit_code == 1
