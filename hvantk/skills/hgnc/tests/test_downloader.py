"""
Unit tests for HGNC downloader.

These tests mock the network calls to avoid actual downloads.
"""

import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock
import tempfile
import shutil

from hvantk.skills.hgnc.cli import download_hgnc
from hvantk.core.constants import HGNC_DOWNLOAD_URL


class TestDownloadHgnc:
    """Tests for the download_hgnc function."""

    def test_download_creates_file(self, tmp_path):
        """Test that download creates the output file."""
        output_path = tmp_path / "hgnc_complete_set.txt"

        # Mock urllib.request.urlretrieve
        with patch("urllib.request.urlretrieve") as mock_retrieve:
            # Simulate successful download by creating the file
            def create_file(_url, path):
                Path(path).write_text("test content")
                return (path, None)

            mock_retrieve.side_effect = create_file

            result = download_hgnc(str(output_path))

            mock_retrieve.assert_called_once_with(HGNC_DOWNLOAD_URL, output_path)
            assert result == str(output_path)
            assert output_path.exists()

    def test_download_creates_parent_directories(self, tmp_path):
        """Test that download creates parent directories if needed."""
        output_path = tmp_path / "nested" / "dir" / "hgnc_complete_set.txt"

        with patch("urllib.request.urlretrieve") as mock_retrieve:

            def create_file(_url, path):
                Path(path).write_text("test content")
                return (path, None)

            mock_retrieve.side_effect = create_file

            result = download_hgnc(str(output_path))

            assert output_path.parent.exists()
            assert result == str(output_path)

    def test_download_raises_on_existing_file(self, tmp_path):
        """Test that download raises FileExistsError if file exists and overwrite=False."""
        output_path = tmp_path / "hgnc_complete_set.txt"
        output_path.write_text("existing content")

        with pytest.raises(FileExistsError):
            download_hgnc(str(output_path), overwrite=False)

    def test_download_overwrites_existing_file(self, tmp_path):
        """Test that download overwrites file when overwrite=True."""
        output_path = tmp_path / "hgnc_complete_set.txt"
        output_path.write_text("existing content")

        with patch("urllib.request.urlretrieve") as mock_retrieve:

            def create_file(_url, path):
                Path(path).write_text("new content")
                return (path, None)

            mock_retrieve.side_effect = create_file

            result = download_hgnc(str(output_path), overwrite=True)

            assert result == str(output_path)
            assert output_path.read_text() == "new content"

    def test_download_handles_network_error(self, tmp_path):
        """Test that download wraps network errors in RuntimeError."""
        import urllib.error

        output_path = tmp_path / "hgnc_complete_set.txt"

        with patch("urllib.request.urlretrieve") as mock_retrieve:
            mock_retrieve.side_effect = urllib.error.URLError("Connection refused")

            with pytest.raises(RuntimeError, match="Failed to download HGNC data"):
                download_hgnc(str(output_path))


class TestHgncDownloaderCli:
    """Tests for the CLI command."""

    def test_cli_help(self):
        """Test that CLI help works."""
        from click.testing import CliRunner
        from hvantk.skills.hgnc.cli import download_cmd as hgnc_downloader

        runner = CliRunner()
        result = runner.invoke(hgnc_downloader, ["--help"])

        assert result.exit_code == 0
        assert "Download the HGNC complete gene nomenclature dataset" in result.output

    def test_cli_requires_output(self):
        """Test that CLI requires --output option."""
        from click.testing import CliRunner
        from hvantk.skills.hgnc.cli import download_cmd as hgnc_downloader

        runner = CliRunner()
        result = runner.invoke(hgnc_downloader, [])

        assert result.exit_code != 0
        assert "Missing option" in result.output or "required" in result.output.lower()

    def test_cli_download_success(self, tmp_path):
        """Test successful CLI download."""
        from click.testing import CliRunner
        from hvantk.skills.hgnc.cli import download_cmd as hgnc_downloader

        output_path = tmp_path / "hgnc_complete_set.txt"

        with patch("hvantk.skills.hgnc.cli.download_hgnc") as mock_download:
            mock_download.return_value = str(output_path)

            runner = CliRunner()
            result = runner.invoke(hgnc_downloader, ["--output", str(output_path)])

            assert result.exit_code == 0
            assert "Downloaded to:" in result.output
            mock_download.assert_called_once_with(str(output_path), overwrite=False)

    def test_cli_file_exists_error(self, tmp_path):
        """Test CLI handles file exists error."""
        from click.testing import CliRunner
        from hvantk.skills.hgnc.cli import download_cmd as hgnc_downloader

        output_path = tmp_path / "hgnc_complete_set.txt"
        output_path.write_text("existing")

        with patch("hvantk.skills.hgnc.cli.download_hgnc") as mock_download:
            mock_download.side_effect = FileExistsError(
                f"File already exists: {output_path}"
            )

            runner = CliRunner()
            result = runner.invoke(hgnc_downloader, ["--output", str(output_path)])

            assert result.exit_code == 1
            assert "Use --overwrite" in result.output
