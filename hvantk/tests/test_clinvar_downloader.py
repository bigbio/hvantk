"""
Tests for ClinVar VCF downloader.

Includes unit tests with mocks and optional network tests.
"""

import hashlib
import os
import tempfile
from unittest.mock import patch, MagicMock

import pytest
from click.testing import CliRunner

from hvantk.commands.clinvar_downloader import clinvar_downloader
from hvantk.datasets.clinvar_datasets import ClinVarDataset


class TestClinVarDataset:
    """Tests for the ClinVarDataset class."""

    def test_latest_grch38(self):
        """Test creating latest GRCh38 dataset."""
        ds = ClinVarDataset.latest()
        assert ds.genome_build == "GRCh38"
        assert ds.version_date is None
        assert ds.file_name == "clinvar.vcf.gz"
        assert "vcf_GRCh38/clinvar.vcf.gz" in ds.download_url

    def test_latest_grch37(self):
        """Test creating latest GRCh37 dataset."""
        ds = ClinVarDataset.latest(genome_build="GRCh37")
        assert ds.genome_build == "GRCh37"
        assert "vcf_GRCh37/clinvar.vcf.gz" in ds.download_url

    def test_latest_invalid_build(self):
        """Test that unsupported genome build raises ValueError."""
        with pytest.raises(ValueError, match="Unsupported genome build"):
            ClinVarDataset.latest(genome_build="hg19")

    def test_from_date_valid(self):
        """Test creating dataset from valid YYYYMMDD date."""
        ds = ClinVarDataset.from_date("20260101")
        assert ds.version_date == "20260101"
        assert ds.file_name == "clinvar_20260101.vcf.gz"
        assert "archive/clinvar_20260101.vcf.gz" in ds.download_url

    def test_from_date_grch37(self):
        """Test from_date with GRCh37 build."""
        ds = ClinVarDataset.from_date("20260101", genome_build="GRCh37")
        assert "vcf_GRCh37/archive/" in ds.download_url

    def test_from_date_invalid_format(self):
        """Test that non-YYYYMMDD strings are rejected."""
        with pytest.raises(ValueError, match="Invalid version_date format"):
            ClinVarDataset.from_date("2026-01-01")

    def test_from_date_invalid_short(self):
        """Test that short date strings are rejected."""
        with pytest.raises(ValueError, match="Invalid version_date format"):
            ClinVarDataset.from_date("202601")

    def test_get_metadata(self):
        """Test metadata generation for latest dataset."""
        ds = ClinVarDataset.latest()
        meta = ds.get_metadata()
        assert meta["source"] == "ClinVar"
        assert meta["genome_build"] == "GRCh38"
        assert meta["version_date"] == "latest"
        assert "file_name" in meta
        assert "download_url" in meta

    def test_get_metadata_versioned(self):
        """Test metadata generation for versioned dataset."""
        ds = ClinVarDataset.from_date("20260101")
        meta = ds.get_metadata()
        assert meta["version_date"] == "20260101"

    def test_str_representation(self):
        """Test string representation."""
        ds = ClinVarDataset.latest()
        assert "latest" in str(ds)
        assert "GRCh38" in str(ds)

    def test_download_creates_directory(self):
        """Test that download creates output directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = os.path.join(tmpdir, "new_subdir", "clinvar")
            ds = ClinVarDataset.latest()

            with patch("hvantk.datasets.clinvar_datasets.download_file") as mock_dl:
                mock_dl.return_value = os.path.join(output_dir, ds.file_name)
                ds.download(output_dir, download_index=False)
                assert os.path.exists(output_dir)

    def test_download_file_exists_error(self):
        """Test that download raises error if file exists without overwrite."""
        with tempfile.TemporaryDirectory() as tmpdir:
            ds = ClinVarDataset.latest()
            file_path = os.path.join(tmpdir, ds.file_name)

            with open(file_path, "w") as f:
                f.write("test")

            with pytest.raises(FileExistsError):
                ds.download(tmpdir, overwrite=False)

    def test_download_overwrite(self):
        """Test that download with overwrite succeeds."""
        with tempfile.TemporaryDirectory() as tmpdir:
            ds = ClinVarDataset.latest()
            file_path = os.path.join(tmpdir, ds.file_name)

            with open(file_path, "w") as f:
                f.write("test")

            with patch("hvantk.datasets.clinvar_datasets.download_file") as mock_dl:
                mock_dl.return_value = file_path
                result = ds.download(tmpdir, overwrite=True, download_index=False)
                assert result == file_path
                mock_dl.assert_called_once()

    def test_download_with_index(self):
        """Test that download_index=True triggers a second download call."""
        with tempfile.TemporaryDirectory() as tmpdir:
            ds = ClinVarDataset.latest()

            with patch("hvantk.datasets.clinvar_datasets.download_file") as mock_dl:
                ds.download(tmpdir, download_index=True)
                # Two calls: one for VCF, one for .tbi
                assert mock_dl.call_count == 2
                tbi_call = mock_dl.call_args_list[1]
                assert tbi_call.kwargs["file_name"].endswith(".tbi")

    def test_download_without_index(self):
        """Test that download_index=False skips the .tbi download."""
        with tempfile.TemporaryDirectory() as tmpdir:
            ds = ClinVarDataset.latest()

            with patch("hvantk.datasets.clinvar_datasets.download_file") as mock_dl:
                ds.download(tmpdir, download_index=False)
                mock_dl.assert_called_once()

    def test_verify_md5_match(self):
        """Test MD5 verification passes when checksums match."""
        with tempfile.TemporaryDirectory() as tmpdir:
            ds = ClinVarDataset.latest()

            # Create a test file and compute its MD5
            test_file = os.path.join(tmpdir, "test.vcf.gz")
            content = b"test content for md5"
            with open(test_file, "wb") as f:
                f.write(content)

            expected_md5 = hashlib.md5(content).hexdigest()
            md5_response = f"{expected_md5}  clinvar.vcf.gz"

            with patch("urllib.request.urlopen") as mock_urlopen:
                mock_resp = MagicMock()
                mock_resp.read.return_value = md5_response.encode("utf-8")
                mock_resp.__enter__ = MagicMock(return_value=mock_resp)
                mock_resp.__exit__ = MagicMock(return_value=False)
                mock_urlopen.return_value = mock_resp

                assert ds.verify_md5(test_file) is True

    def test_verify_md5_mismatch(self):
        """Test MD5 verification fails when checksums don't match."""
        with tempfile.TemporaryDirectory() as tmpdir:
            ds = ClinVarDataset.latest()

            test_file = os.path.join(tmpdir, "test.vcf.gz")
            with open(test_file, "wb") as f:
                f.write(b"test content")

            md5_response = "0000000000000000000000000000dead  clinvar.vcf.gz"

            with patch("urllib.request.urlopen") as mock_urlopen:
                mock_resp = MagicMock()
                mock_resp.read.return_value = md5_response.encode("utf-8")
                mock_resp.__enter__ = MagicMock(return_value=mock_resp)
                mock_resp.__exit__ = MagicMock(return_value=False)
                mock_urlopen.return_value = mock_resp

                assert ds.verify_md5(test_file) is False


class TestClinVarDownloaderCLI:
    """Tests for the clinvar-downloader CLI command."""

    def test_download_latest(self):
        """Test default invocation downloads latest."""
        runner = CliRunner()
        with tempfile.TemporaryDirectory() as tmpdir:
            with patch("hvantk.datasets.clinvar_datasets.download_file") as mock_dl:
                result = runner.invoke(
                    clinvar_downloader,
                    ["--output-dir", tmpdir],
                )
                assert result.exit_code == 0
                assert "Downloaded to:" in result.output
                mock_dl.assert_called()

    def test_download_specific_version(self):
        """Test --version YYYYMMDD works."""
        runner = CliRunner()
        with tempfile.TemporaryDirectory() as tmpdir:
            with patch("hvantk.datasets.clinvar_datasets.download_file") as mock_dl:
                result = runner.invoke(
                    clinvar_downloader,
                    ["--version", "20260101", "--output-dir", tmpdir],
                )
                assert result.exit_code == 0
                assert "clinvar_20260101" in result.output

    def test_download_grch37(self):
        """Test --genome-build GRCh37 works."""
        runner = CliRunner()
        with tempfile.TemporaryDirectory() as tmpdir:
            with patch("hvantk.datasets.clinvar_datasets.download_file") as mock_dl:
                result = runner.invoke(
                    clinvar_downloader,
                    ["--genome-build", "GRCh37", "--output-dir", tmpdir],
                )
                assert result.exit_code == 0
                assert "GRCh37" in result.output

    def test_no_index_flag(self):
        """Test --no-index skips .tbi download."""
        runner = CliRunner()
        with tempfile.TemporaryDirectory() as tmpdir:
            with patch("hvantk.datasets.clinvar_datasets.download_file") as mock_dl:
                result = runner.invoke(
                    clinvar_downloader,
                    ["--no-index", "--output-dir", tmpdir],
                )
                assert result.exit_code == 0
                # Only one call (VCF), no .tbi
                mock_dl.assert_called_once()

    def test_overwrite_flag(self):
        """Test --overwrite is passed correctly."""
        runner = CliRunner()
        with tempfile.TemporaryDirectory() as tmpdir:
            with patch.object(ClinVarDataset, "download") as mock_download:
                mock_download.return_value = os.path.join(tmpdir, "clinvar.vcf.gz")

                result = runner.invoke(
                    clinvar_downloader,
                    ["--overwrite", "--output-dir", tmpdir],
                )

                assert result.exit_code == 0
                call_kwargs = mock_download.call_args.kwargs
                assert call_kwargs["overwrite"] is True

    def test_invalid_version_format(self):
        """Test error on non-YYYYMMDD version string."""
        runner = CliRunner()
        result = runner.invoke(
            clinvar_downloader,
            ["--version", "2026-01-01", "--output-dir", "/tmp"],
        )
        assert result.exit_code == 1
        assert "Invalid version_date format" in result.output

    def test_invalid_genome_build(self):
        """Test error on unsupported genome build (caught by Click Choice)."""
        runner = CliRunner()
        result = runner.invoke(
            clinvar_downloader,
            ["--genome-build", "hg19"],
        )
        assert result.exit_code != 0

    def test_verify_md5_flag(self):
        """Test --verify-md5 triggers checksum verification."""
        runner = CliRunner()
        with tempfile.TemporaryDirectory() as tmpdir:
            with patch("hvantk.datasets.clinvar_datasets.download_file"):
                with patch.object(
                    ClinVarDataset, "verify_md5", return_value=True
                ) as mock_md5:
                    result = runner.invoke(
                        clinvar_downloader,
                        ["--verify-md5", "--output-dir", tmpdir],
                    )
                    assert result.exit_code == 0
                    assert "MD5 checksum verified" in result.output
                    mock_md5.assert_called_once()


@pytest.mark.network
class TestClinVarDownloaderNetwork:
    """Network-dependent tests for ClinVar downloader.

    These tests require network access and are skipped by default.
    Run with: pytest -m network
    """

    def test_latest_url_reachable(self):
        """Test that the latest ClinVar VCF URL is reachable (HEAD request)."""
        import urllib.request

        ds = ClinVarDataset.latest()
        req = urllib.request.Request(
            ds.download_url,
            method="HEAD",
            headers={"User-Agent": "hvantk/1.0"},
        )
        with urllib.request.urlopen(req, timeout=30) as resp:
            assert resp.status == 200
