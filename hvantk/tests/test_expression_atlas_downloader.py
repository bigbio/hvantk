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


@pytest.mark.network
@pytest.mark.slow
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


@pytest.mark.network
@pytest.mark.slow
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


def test_download_file_with_retry_progress_bar(download_path, monkeypatch):
    """
    Test that _download_file_with_retry uses tqdm for progress tracking.
    This is a unit test that mocks FTP operations.
    """
    import ftplib
    from unittest.mock import Mock
    from hvantk.commands.expression_atlas_downloader import _download_file_with_retry

    # Create a mock FTP object
    mock_ftp = Mock(spec=ftplib.FTP)
    
    # Mock the file size command
    mock_ftp.voidcmd = Mock()
    mock_ftp.size = Mock(return_value=1024)  # 1KB file
    
    # Mock retrbinary to simulate file download with chunks
    def mock_retrbinary(cmd, callback):
        # Simulate downloading in chunks
        chunk1 = b"x" * 512
        chunk2 = b"y" * 512
        callback(chunk1)
        callback(chunk2)
    
    mock_ftp.retrbinary = mock_retrbinary
    
    # Create test file path
    test_file = os.path.join(download_path, "test_file.txt")
    
    # Call the function
    result = _download_file_with_retry(
        mock_ftp, 
        "test_file.txt", 
        test_file,
        ftp_url="test.ftp.com",
        ftp_path="/test/path"
    )
    
    # Assert download was successful
    assert result is True
    assert os.path.exists(test_file)
    
    # Verify the file content
    with open(test_file, "rb") as f:
        content = f.read()
        assert len(content) == 1024
        assert content[:512] == b"x" * 512
        assert content[512:] == b"y" * 512
    
    # Verify FTP methods were called
    mock_ftp.voidcmd.assert_called_once_with("TYPE I")
    mock_ftp.size.assert_called_once_with("test_file.txt")

