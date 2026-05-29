"""
ClinVar VCF dataset handling.

This module provides classes for downloading and managing ClinVar VCF datasets.
ClinVar publishes monthly VCF releases on the NCBI FTP server for both GRCh38
and GRCh37 genome builds.

Example usage:
    # Download latest ClinVar VCF
    dataset = ClinVarDataset.latest()
    dataset.download("/data/clinvar")

    # Download a specific archived version
    dataset = ClinVarDataset.from_date("20260101")
    dataset.download("/data/clinvar")

    # Download GRCh37 build
    dataset = ClinVarDataset.latest(genome_build="GRCh37")
    dataset.download("/data/clinvar")
"""

import hashlib
import logging
import os
import re
from dataclasses import dataclass
from typing import Dict, Optional

from hvantk.skills.clinvar.shared.constants import CLINVAR_FTP_BASE, CLINVAR_FTP_BASE_GRCh37
from hvantk.core.utils.file_utils import download_file

logger = logging.getLogger(__name__)

_BASE_URLS = {
    "GRCh38": CLINVAR_FTP_BASE,
    "GRCh37": CLINVAR_FTP_BASE_GRCh37,
}


@dataclass
class ClinVarDataset:
    """
    Represents a ClinVar VCF dataset.

    ClinVar provides a latest VCF (``clinvar.vcf.gz``) and dated archives
    (``clinvar_YYYYMMDD.vcf.gz``) under an ``archive/`` subdirectory.

    Attributes:
        genome_build: Reference genome ("GRCh38" or "GRCh37").
        version_date: Archive date string (YYYYMMDD), or None for latest.
        download_url: Full URL to the VCF file.
        file_name: Name of the VCF file.
    """

    genome_build: str
    version_date: Optional[str]
    download_url: str
    file_name: str

    @classmethod
    def latest(cls, genome_build: str = "GRCh38") -> "ClinVarDataset":
        """
        Reference the current/latest ClinVar VCF.

        Args:
            genome_build: "GRCh38" or "GRCh37".

        Returns:
            ClinVarDataset pointing to the latest VCF.

        Raises:
            ValueError: If genome_build is not supported.
        """
        base = _get_base_url(genome_build)
        file_name = "clinvar.vcf.gz"
        download_url = f"{base}/{file_name}"
        return cls(
            genome_build=genome_build,
            version_date=None,
            download_url=download_url,
            file_name=file_name,
        )

    @classmethod
    def from_date(
        cls, version_date: str, genome_build: str = "GRCh38"
    ) -> "ClinVarDataset":
        """
        Reference an archived ClinVar release by date.

        Args:
            version_date: Date string in YYYYMMDD format.
            genome_build: "GRCh38" or "GRCh37".

        Returns:
            ClinVarDataset pointing to the archived VCF.

        Raises:
            ValueError: If version_date format is invalid or genome_build unsupported.
        """
        if not re.match(r"^\d{8}$", version_date):
            raise ValueError(
                f"Invalid version_date format: {version_date}. Expected YYYYMMDD (e.g. 20260101)"
            )

        base = _get_base_url(genome_build)
        file_name = f"clinvar_{version_date}.vcf.gz"
        download_url = f"{base}/archive/{file_name}"
        return cls(
            genome_build=genome_build,
            version_date=version_date,
            download_url=download_url,
            file_name=file_name,
        )

    def download(
        self,
        output_dir: str,
        overwrite: bool = False,
        download_index: bool = True,
    ) -> str:
        """
        Download the VCF file (and optionally the .tbi index).

        Args:
            output_dir: Directory to save the downloaded file(s).
            overwrite: If True, overwrite existing files.
            download_index: If True, also download the ``.tbi`` tabix index.

        Returns:
            Path to the downloaded VCF file.

        Raises:
            FileExistsError: If file exists and overwrite is False.
            RuntimeError: If the download fails.
        """
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, self.file_name)

        if os.path.exists(output_path) and not overwrite:
            raise FileExistsError(
                f"File already exists: {output_path}. Use overwrite=True to replace."
            )

        logger.info(f"Downloading ClinVar VCF from {self.download_url}")
        try:
            download_file(
                url=self.download_url,
                out_dir=output_dir,
                file_name=self.file_name,
            )
            logger.info(f"Downloaded ClinVar VCF to {output_path}")
        except Exception as e:
            raise RuntimeError(f"Failed to download ClinVar VCF: {e}") from e

        if download_index:
            self._download_index(output_dir, overwrite)

        return output_path

    def verify_md5(self, file_path: str) -> bool:
        """
        Download the ``.md5`` checksum from NCBI and verify against a local file.

        Args:
            file_path: Path to the local VCF file to verify.

        Returns:
            True if the checksum matches, False otherwise.

        Raises:
            RuntimeError: If the MD5 file cannot be downloaded.
        """
        import urllib.request

        md5_url = f"{self.download_url}.md5"
        logger.info(f"Fetching MD5 checksum from {md5_url}")

        try:
            req = urllib.request.Request(md5_url, headers={"User-Agent": "hvantk/1.0"})
            with urllib.request.urlopen(req, timeout=30) as response:
                md5_content = response.read().decode("utf-8").strip()
        except Exception as e:
            raise RuntimeError(f"Failed to fetch MD5 checksum: {e}") from e

        # NCBI .md5 files contain lines like: "d41d8cd9...  clinvar.vcf.gz"
        parts = md5_content.split()
        if not parts:
            raise RuntimeError(f"MD5 file is empty or malformed: {md5_url}")
        expected_hash = parts[0].lower()

        logger.info(f"Computing MD5 of {file_path}")
        # MD5 is the algorithm NCBI publishes for ClinVar VCF integrity (e.g.
        # clinvar.vcf.gz.md5). We are verifying an upstream-published checksum,
        # not generating a cryptographic signature; SHA256 cannot be substituted
        # without the upstream switching first. Suppress the security warning.
        md5 = hashlib.md5(usedforsecurity=False)  # noqa: S324
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(8192), b""):
                md5.update(chunk)
        actual_hash = md5.hexdigest().lower()

        if actual_hash == expected_hash:
            logger.info("MD5 checksum verified successfully")
            return True
        else:
            logger.warning(f"MD5 mismatch: expected {expected_hash}, got {actual_hash}")
            return False

    def get_metadata(self) -> Dict[str, Optional[str]]:
        """Get metadata about this dataset."""
        return {
            "source": "ClinVar",
            "genome_build": self.genome_build,
            "version_date": self.version_date or "latest",
            "download_url": self.download_url,
            "file_name": self.file_name,
            "description": (
                "Clinically relevant variant annotations from ClinVar (VCF format). "
                "Includes pathogenicity classifications, review status, and disease associations."
            ),
        }

    def _download_index(self, output_dir: str, overwrite: bool) -> None:
        """Download the .tbi tabix index file."""
        tbi_name = f"{self.file_name}.tbi"
        tbi_url = f"{self.download_url}.tbi"
        tbi_path = os.path.join(output_dir, tbi_name)

        if os.path.exists(tbi_path) and not overwrite:
            logger.info(f"Index file already exists: {tbi_path}, skipping")
            return

        logger.info(f"Downloading tabix index from {tbi_url}")
        try:
            download_file(url=tbi_url, out_dir=output_dir, file_name=tbi_name)
            logger.info(f"Downloaded tabix index to {tbi_path}")
        except Exception as e:
            logger.warning(f"Failed to download tabix index: {e}")

    def __str__(self) -> str:
        version = self.version_date or "latest"
        return f"ClinVarDataset({self.genome_build}, version={version})"


def _get_base_url(genome_build: str) -> str:
    """Get the FTP base URL for the given genome build."""
    if genome_build not in _BASE_URLS:
        raise ValueError(
            f"Unsupported genome build: {genome_build}. "
            f"Supported: {list(_BASE_URLS.keys())}"
        )
    return _BASE_URLS[genome_build]
