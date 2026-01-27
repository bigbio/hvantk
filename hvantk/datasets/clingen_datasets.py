"""
ClinGen Gene-Disease Validity dataset handling.

This module provides classes for downloading and managing ClinGen Gene-Disease
Validity datasets. ClinGen curates gene-disease associations with evidence-based
classifications (Definitive, Strong, Moderate, Limited, etc.).

Example usage:
    # Get latest dataset
    dataset = ClinGenGeneDiseaseDataset.from_latest()
    dataset.download("/data/clingen")

    # Get specific version by date
    dataset = ClinGenGeneDiseaseDataset.from_date("2026-01-15")
    dataset.download("/data/clingen")

    # List available versions
    versions = get_available_versions()
"""

import logging
import os
import re
from dataclasses import dataclass
from datetime import datetime
from typing import Dict, List, Optional

from hvantk.core.constants import (
    CLINGEN_BASE_URL,
    CLINGEN_DOWNLOADS_URL,
    CLINGEN_FILE_PREFIX,
)
from hvantk.data.file_utils import download_file

logger = logging.getLogger(__name__)


@dataclass
class ClinGenGeneDiseaseDataset:
    """
    Represents a ClinGen Gene-Disease Validity dataset.

    ClinGen provides versioned CSV files containing gene-disease associations
    with evidence-based classifications.

    Attributes:
        version_date: Date of the dataset version (YYYY-MM-DD format)
        download_url: URL to download the CSV file
        file_name: Name of the CSV file
    """

    version_date: str
    download_url: str
    file_name: str

    @classmethod
    def from_date(cls, version_date: str) -> "ClinGenGeneDiseaseDataset":
        """
        Create a dataset reference for a specific version date.

        Args:
            version_date: Date string in YYYY-MM-DD format

        Returns:
            ClinGenGeneDiseaseDataset instance

        Raises:
            ValueError: If version_date is not in valid format
        """
        # Validate date format
        try:
            datetime.strptime(version_date, "%Y-%m-%d")
        except ValueError:
            raise ValueError(
                f"Invalid version_date format: {version_date}. Expected YYYY-MM-DD"
            )

        file_name = f"{CLINGEN_FILE_PREFIX}-{version_date}.csv"
        download_url = f"{CLINGEN_BASE_URL}?file={file_name}"

        return cls(
            version_date=version_date,
            download_url=download_url,
            file_name=file_name,
        )

    @classmethod
    def from_latest(cls) -> "ClinGenGeneDiseaseDataset":
        """
        Create a dataset reference for the latest available version.

        This method fetches available versions from the ClinGen portal
        and returns the most recent one.

        Returns:
            ClinGenGeneDiseaseDataset instance for the latest version

        Raises:
            RuntimeError: If unable to determine the latest version
        """
        versions = get_available_versions()
        if not versions:
            raise RuntimeError(
                "Unable to determine available ClinGen versions. "
                "Network may be unavailable or ClinGen portal structure changed."
            )

        latest_date = versions[0]  # Versions are sorted newest first
        logger.info(f"Latest ClinGen Gene-Disease version: {latest_date}")
        return cls.from_date(latest_date)

    def download(
        self, output_dir: str, overwrite: bool = False
    ) -> str:
        """
        Download the dataset CSV file.

        Args:
            output_dir: Directory to save the downloaded file
            overwrite: If True, overwrite existing file

        Returns:
            Path to the downloaded file

        Raises:
            FileExistsError: If file exists and overwrite=False
            RuntimeError: If download fails
        """
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, self.file_name)

        if os.path.exists(output_path) and not overwrite:
            raise FileExistsError(
                f"File already exists: {output_path}. Use overwrite=True to replace."
            )

        logger.info(f"Downloading ClinGen dataset from {self.download_url}")
        try:
            download_file(
                url=self.download_url,
                out_dir=output_dir,
                file_name=self.file_name,
            )
            logger.info(f"Downloaded ClinGen dataset to {output_path}")
            return output_path
        except Exception as e:
            raise RuntimeError(
                f"Failed to download ClinGen dataset: {str(e)}"
            ) from e

    def get_metadata(self) -> Dict[str, str]:
        """
        Get metadata about this dataset.

        Returns:
            Dictionary with dataset metadata
        """
        return {
            "source": "ClinGen Gene-Disease Validity",
            "version_date": self.version_date,
            "download_url": self.download_url,
            "file_name": self.file_name,
            "description": (
                "Gene-disease validity classifications curated by ClinGen. "
                "Classifications include Definitive, Strong, Moderate, Limited, "
                "Disputed, Refuted, and No Known Disease Relationship."
            ),
        }

    def __str__(self) -> str:
        return f"ClinGenGeneDiseaseDataset(version={self.version_date})"


def get_available_versions() -> List[str]:
    """
    Fetch available ClinGen Gene-Disease dataset versions.

    Scrapes the ClinGen downloads page to find available version dates.
    Returns dates sorted from newest to oldest.

    Returns:
        List of version dates in YYYY-MM-DD format, sorted newest first

    Note:
        This function requires network access. Returns empty list if
        the ClinGen portal is unavailable or format has changed.
    """
    try:
        import urllib.request

        logger.info(f"Fetching available versions from {CLINGEN_DOWNLOADS_URL}")

        req = urllib.request.Request(
            CLINGEN_DOWNLOADS_URL,
            headers={"User-Agent": "hvantk/1.0"},
        )

        with urllib.request.urlopen(req, timeout=30) as response:
            html = response.read().decode("utf-8")

        # Extract dates from file names like "Clingen-Gene-Disease-Summary-YYYY-MM-DD.csv"
        pattern = rf"{CLINGEN_FILE_PREFIX}-(\d{{4}}-\d{{2}}-\d{{2}})\.csv"
        matches = re.findall(pattern, html)

        if not matches:
            logger.warning("No ClinGen versions found on downloads page")
            return []

        # Remove duplicates and sort newest first
        unique_dates = sorted(set(matches), reverse=True)
        logger.info(f"Found {len(unique_dates)} ClinGen versions")
        return unique_dates

    except Exception as e:
        logger.error(f"Failed to fetch ClinGen versions: {e}")
        return []


def get_latest_version() -> Optional[str]:
    """
    Get the latest available ClinGen Gene-Disease dataset version date.

    Returns:
        Version date string (YYYY-MM-DD) or None if unavailable
    """
    versions = get_available_versions()
    return versions[0] if versions else None
