"""
ClinGen Gene-Disease Validity dataset handling.

This module provides classes for downloading and managing ClinGen Gene-Disease
Validity datasets. ClinGen curates gene-disease associations with evidence-based
classifications (Definitive, Strong, Moderate, Limited, etc.).

Example usage:
    # Get latest dataset (today's snapshot)
    dataset = ClinGenGeneDiseaseDataset.from_latest()
    dataset.download("/data/clingen")

    # Check availability
    versions = get_available_versions()
"""

import logging
import os
from dataclasses import dataclass
from datetime import datetime
from typing import Dict, List, Optional

from hvantk.core.constants import (
    CLINGEN_BASE_URL,
    CLINGEN_FILE_PREFIX,
)
from hvantk.core.utils.file_utils import download_file

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
        Create a dataset reference labeled with a specific date.

        Note: ClinGen now provides a single real-time download endpoint.
        The date is used only for labeling the output file. The downloaded
        content will always be the current snapshot regardless of the date
        provided.

        Args:
            version_date: Date string in YYYY-MM-DD format

        Returns:
            ClinGenGeneDiseaseDataset instance

        Raises:
            ValueError: If version_date is not in valid format
        """
        # Validate date format
        try:
            parsed = datetime.strptime(version_date, "%Y-%m-%d")
        except ValueError:
            raise ValueError(
                f"Invalid version_date format: {version_date}. Expected YYYY-MM-DD"
            )

        file_name = f"{CLINGEN_FILE_PREFIX}-{version_date}.csv"

        return cls(
            version_date=version_date,
            download_url=CLINGEN_BASE_URL,
            file_name=file_name,
        )

    @classmethod
    def from_latest(cls) -> "ClinGenGeneDiseaseDataset":
        """
        Create a dataset reference for the latest available snapshot.

        ClinGen provides a real-time download endpoint that generates
        the CSV on-the-fly with current data. The version date is set
        to today's date.

        Returns:
            ClinGenGeneDiseaseDataset instance for today's snapshot
        """
        today = datetime.now().strftime("%Y-%m-%d")
        logger.info(f"Using ClinGen real-time download, labeling as {today}")
        return cls.from_date(today)

    def download(self, output_dir: str, overwrite: bool = False) -> str:
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
            raise RuntimeError(f"Failed to download ClinGen dataset: {str(e)}") from e

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
    Check ClinGen Gene-Disease dataset availability.

    ClinGen no longer provides versioned archives on their downloads page.
    The Gene-Disease Validity CSV is generated in real-time from the endpoint.
    This function verifies the download endpoint is reachable and returns
    today's date as the available version.

    Returns:
        List with today's date if the endpoint is reachable, empty list otherwise

    Note:
        This function requires network access.
    """
    try:
        import urllib.request

        logger.info(f"Checking ClinGen download endpoint: {CLINGEN_BASE_URL}")

        req = urllib.request.Request(
            CLINGEN_BASE_URL,
            method="HEAD",
            headers={"User-Agent": "hvantk/1.0"},
        )

        with urllib.request.urlopen(req, timeout=30) as response:
            if response.status == 200:
                today = datetime.now().strftime("%Y-%m-%d")
                logger.info(
                    "ClinGen endpoint is reachable. "
                    "Dataset is generated in real-time (no versioned archives)."
                )
                return [today]

        return []

    except Exception as e:
        logger.error(f"Failed to reach ClinGen endpoint: {e}")
        return []


def get_latest_version() -> Optional[str]:
    """
    Get the latest available ClinGen Gene-Disease dataset version date.

    Since ClinGen now provides real-time snapshots, this returns today's
    date if the endpoint is reachable.

    Returns:
        Today's date string (YYYY-MM-DD) or None if endpoint is unavailable
    """
    versions = get_available_versions()
    return versions[0] if versions else None
