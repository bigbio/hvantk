"""
GenCC (Gene Curation Coalition) submissions dataset handling.

This module provides classes for downloading and managing GenCC submissions
data. GenCC aggregates gene-disease validity assertions from 12+ submitting
organizations (ClinGen, PanelApp, G2P, Orphanet, etc.).

Example usage:
    dataset = GenCCSubmissionsDataset.from_latest()
    dataset.download("/data/gencc")

    versions = get_available_versions()
"""

import logging
import os
from dataclasses import dataclass
from datetime import datetime
from typing import Dict, List, Optional

from hvantk.core.constants import GENCC_BASE_URL, GENCC_FILE_PREFIX
from hvantk.data.file_utils import download_file

logger = logging.getLogger(__name__)


@dataclass
class GenCCSubmissionsDataset:
    """Represents a GenCC submissions dataset.

    GenCC provides a single TSV download endpoint containing gene-disease
    validity assertions from multiple submitting organizations.

    Attributes:
        version_date: Date label for the snapshot (YYYY-MM-DD format)
        download_url: URL to download the TSV file
        file_name: Name of the TSV file
    """

    version_date: str
    download_url: str
    file_name: str

    @classmethod
    def from_date(cls, version_date: str) -> "GenCCSubmissionsDataset":
        """Create a dataset reference labeled with a specific date.

        Note: GenCC provides a single real-time download endpoint.
        The date is used only for labeling the output file.

        Args:
            version_date: Date string in YYYY-MM-DD format

        Returns:
            GenCCSubmissionsDataset instance

        Raises:
            ValueError: If version_date is not in valid format
        """
        try:
            datetime.strptime(version_date, "%Y-%m-%d")
        except ValueError:
            raise ValueError(
                f"Invalid version_date format: {version_date}. Expected YYYY-MM-DD"
            ) from None

        file_name = f"{GENCC_FILE_PREFIX}-{version_date}.tsv"

        return cls(
            version_date=version_date,
            download_url=GENCC_BASE_URL,
            file_name=file_name,
        )

    @classmethod
    def from_latest(cls) -> "GenCCSubmissionsDataset":
        """Create a dataset reference for the latest available snapshot.

        Returns:
            GenCCSubmissionsDataset instance for today's snapshot
        """
        today = datetime.now().strftime("%Y-%m-%d")
        logger.info(f"Using GenCC real-time download, labeling as {today}")
        return cls.from_date(today)

    def download(self, output_dir: str, overwrite: bool = False) -> str:
        """Download the dataset TSV file.

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

        logger.info(f"Downloading GenCC dataset from {self.download_url}")
        try:
            download_file(
                url=self.download_url,
                out_dir=output_dir,
                file_name=self.file_name,
            )
            logger.info(f"Downloaded GenCC dataset to {output_path}")
            return output_path
        except Exception as e:
            raise RuntimeError(f"Failed to download GenCC dataset: {e!s}") from e

    def get_metadata(self) -> Dict[str, str]:
        """Get metadata about this dataset."""
        return {
            "source": "GenCC (Gene Curation Coalition)",
            "version_date": self.version_date,
            "download_url": self.download_url,
            "file_name": self.file_name,
            "description": (
                "Gene-disease validity assertions aggregated from 12+ submitting "
                "organizations including ClinGen, PanelApp, G2P, and Orphanet."
            ),
        }

    def __str__(self) -> str:
        return f"GenCCSubmissionsDataset(version={self.version_date})"


def get_available_versions() -> List[str]:
    """Check GenCC submissions dataset availability.

    GenCC provides a real-time download endpoint. This function verifies
    the endpoint is reachable and returns today's date as the available
    version.

    Returns:
        List with today's date if the endpoint is reachable, empty list
        otherwise.
    """
    try:
        import urllib.request
        from urllib.parse import urlparse

        logger.info(f"Checking GenCC download endpoint: {GENCC_BASE_URL}")

        parsed = urlparse(GENCC_BASE_URL)
        if parsed.scheme not in ("http", "https"):
            logger.error(f"Invalid URL scheme: {parsed.scheme}")
            return []

        req = urllib.request.Request(
            GENCC_BASE_URL,
            method="HEAD",
            headers={"User-Agent": "hvantk/1.0"},
        )

        with urllib.request.urlopen(req, timeout=30) as response:  # noqa: S310
            if response.status == 200:
                today = datetime.now().strftime("%Y-%m-%d")
                logger.info(
                    "GenCC endpoint is reachable. "
                    "Dataset is generated in real-time (no versioned archives)."
                )
                return [today]

        return []

    except Exception as e:
        logger.error(f"Failed to reach GenCC endpoint: {e}")
        return []


def get_latest_version() -> Optional[str]:
    """Get the latest available GenCC dataset version date.

    Returns:
        Today's date string (YYYY-MM-DD) or None if endpoint is unavailable.
    """
    versions = get_available_versions()
    return versions[0] if versions else None
