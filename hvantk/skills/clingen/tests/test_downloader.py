"""
Tests for ClinGen Gene-Disease Validity downloader.

Includes unit tests with mocks and optional network tests.
"""

import os
import tempfile
from unittest.mock import patch, MagicMock

import pytest
from click.testing import CliRunner

from hvantk.skills.clingen.cli import clingen_downloader
from hvantk.skills.clingen.shared.datasets import (
    ClinGenGeneDiseaseDataset,
    get_available_versions,
    get_latest_version,
)


class TestClinGenGeneDiseaseDataset:
    """Tests for the ClinGenGeneDiseaseDataset class."""

    def test_from_date_valid(self):
        """Test creating dataset from valid date."""
        dataset = ClinGenGeneDiseaseDataset.from_date("2026-01-15")
        assert dataset.version_date == "2026-01-15"
        assert "2026-01-15" in dataset.file_name
        assert dataset.download_url  # URL should be set

    def test_from_date_invalid_format(self):
        """Test that invalid date format raises ValueError."""
        with pytest.raises(ValueError, match="Invalid version_date format"):
            ClinGenGeneDiseaseDataset.from_date("01-15-2026")

    def test_from_date_invalid_date(self):
        """Test that invalid date raises ValueError."""
        with pytest.raises(ValueError, match="Invalid version_date format"):
            ClinGenGeneDiseaseDataset.from_date("2026-13-45")

    def test_get_metadata(self):
        """Test metadata generation."""
        dataset = ClinGenGeneDiseaseDataset.from_date("2026-01-15")
        metadata = dataset.get_metadata()
        assert metadata["source"] == "ClinGen Gene-Disease Validity"
        assert metadata["version_date"] == "2026-01-15"
        assert "file_name" in metadata
        assert "download_url" in metadata

    def test_str_representation(self):
        """Test string representation."""
        dataset = ClinGenGeneDiseaseDataset.from_date("2026-01-15")
        assert "2026-01-15" in str(dataset)

    def test_download_creates_directory(self):
        """Test that download creates output directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = os.path.join(tmpdir, "new_subdir", "clingen")
            dataset = ClinGenGeneDiseaseDataset.from_date("2026-01-15")

            with patch("hvantk.skills.clingen.shared.datasets.download_file") as mock_dl:
                mock_dl.return_value = os.path.join(output_dir, dataset.file_name)
                dataset.download(output_dir)
                assert os.path.exists(output_dir)

    def test_download_file_exists_error(self):
        """Test that download raises error if file exists without overwrite."""
        with tempfile.TemporaryDirectory() as tmpdir:
            dataset = ClinGenGeneDiseaseDataset.from_date("2026-01-15")
            file_path = os.path.join(tmpdir, dataset.file_name)

            # Create the file
            with open(file_path, "w") as f:
                f.write("test")

            with pytest.raises(FileExistsError):
                dataset.download(tmpdir, overwrite=False)

    def test_download_overwrite(self):
        """Test that download with overwrite succeeds."""
        with tempfile.TemporaryDirectory() as tmpdir:
            dataset = ClinGenGeneDiseaseDataset.from_date("2026-01-15")
            file_path = os.path.join(tmpdir, dataset.file_name)

            # Create the file
            with open(file_path, "w") as f:
                f.write("test")

            with patch("hvantk.skills.clingen.shared.datasets.download_file") as mock_dl:
                mock_dl.return_value = file_path
                result = dataset.download(tmpdir, overwrite=True)
                assert result == file_path
                mock_dl.assert_called_once()


class TestGetAvailableVersions:
    """Tests for get_available_versions function."""

    def test_get_available_versions_reachable(self):
        """Test version check when endpoint is reachable."""
        with patch("urllib.request.urlopen") as mock_urlopen:
            mock_response = MagicMock()
            mock_response.status = 200
            mock_response.__enter__ = MagicMock(return_value=mock_response)
            mock_response.__exit__ = MagicMock(return_value=False)
            mock_urlopen.return_value = mock_response

            versions = get_available_versions()

            assert len(versions) == 1
            # Should return today's date
            assert len(versions[0].split("-")) == 3

    def test_get_available_versions_network_error(self):
        """Test handling of network errors."""
        with patch("urllib.request.urlopen") as mock_urlopen:
            mock_urlopen.side_effect = Exception("Network error")

            versions = get_available_versions()
            assert versions == []

    def test_get_latest_version_mocked(self):
        """Test get_latest_version with mocked versions."""
        with patch(
            "hvantk.skills.clingen.shared.datasets.get_available_versions"
        ) as mock_versions:
            mock_versions.return_value = ["2026-03-09"]

            latest = get_latest_version()
            assert latest == "2026-03-09"

    def test_get_latest_version_empty(self):
        """Test get_latest_version when no versions available."""
        with patch(
            "hvantk.skills.clingen.shared.datasets.get_available_versions"
        ) as mock_versions:
            mock_versions.return_value = []

            latest = get_latest_version()
            assert latest is None


class TestClinGenDownloaderCLI:
    """Tests for the clingen-downloader CLI command."""

    def test_list_versions(self):
        """Test --list-versions flag."""
        runner = CliRunner()
        with patch(
            "hvantk.skills.clingen.shared.datasets.get_available_versions"
        ) as mock_versions:
            mock_versions.return_value = ["2026-03-09"]

            result = runner.invoke(clingen_downloader, ["--list-versions"])

            assert result.exit_code == 0
            assert "real-time" in result.output
            assert "2026-03-09" in result.output

    def test_download_specific_version(self):
        """Test downloading a specific version."""
        runner = CliRunner()
        with tempfile.TemporaryDirectory() as tmpdir:
            with patch(
                "hvantk.skills.clingen.shared.datasets.download_file"
            ) as mock_download:
                mock_download.return_value = os.path.join(
                    tmpdir, "Clingen-Gene-Disease-Summary-2026-01-15.csv"
                )

                result = runner.invoke(
                    clingen_downloader,
                    ["--version", "2026-01-15", "--output-dir", tmpdir],
                )

                assert result.exit_code == 0
                assert "Downloaded to:" in result.output
                mock_download.assert_called_once()

    def test_download_latest_version(self):
        """Test downloading latest version (real-time snapshot)."""
        runner = CliRunner()
        with tempfile.TemporaryDirectory() as tmpdir:
            with patch(
                "hvantk.skills.clingen.shared.datasets.download_file"
            ) as mock_download:
                mock_download.return_value = os.path.join(
                    tmpdir, "Clingen-Gene-Disease-Summary-2026-03-09.csv"
                )

                result = runner.invoke(
                    clingen_downloader,
                    ["--version", "latest", "--output-dir", tmpdir],
                )

                assert result.exit_code == 0
                mock_download.assert_called_once()

    def test_invalid_version_format(self):
        """Test error on invalid version format."""
        runner = CliRunner()
        result = runner.invoke(
            clingen_downloader,
            ["--version", "invalid-date", "--output-dir", "/tmp"],
        )
        assert result.exit_code == 1
        assert "Invalid version_date format" in result.output

    def test_download_overwrite_flag(self):
        """Test --overwrite flag is passed correctly."""
        runner = CliRunner()
        with tempfile.TemporaryDirectory() as tmpdir:
            with patch.object(ClinGenGeneDiseaseDataset, "download") as mock_download:
                mock_download.return_value = os.path.join(tmpdir, "test.csv")

                result = runner.invoke(
                    clingen_downloader,
                    [
                        "--version",
                        "2026-01-15",
                        "--output-dir",
                        tmpdir,
                        "--overwrite",
                    ],
                )

                assert result.exit_code == 0
                mock_download.assert_called_once()
                call_kwargs = mock_download.call_args.kwargs
                assert call_kwargs["overwrite"] is True


@pytest.mark.network
class TestClinGenDownloaderNetwork:
    """Network-dependent tests for ClinGen downloader.

    These tests require network access and are skipped by default.
    Run with: pytest -m network
    """

    def test_get_available_versions_real(self):
        """Test fetching real versions from ClinGen portal."""
        versions = get_available_versions()
        # Should have at least some versions
        assert len(versions) > 0
        # Versions should be in YYYY-MM-DD format
        for version in versions:
            assert len(version.split("-")) == 3

    def test_from_latest_real(self):
        """Test creating dataset from latest real version."""
        dataset = ClinGenGeneDiseaseDataset.from_latest()
        assert dataset.version_date is not None
        assert dataset.download_url is not None
