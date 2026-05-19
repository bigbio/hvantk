"""
UniProt PTM (Post-Translational Modification) dataset handling.

This module provides classes for downloading and managing human PTM data
from the UniProt REST API. UniProt curates experimentally verified and
predicted post-translational modifications for reviewed (Swiss-Prot) proteins.

Example usage:
    # Get latest dataset (today's snapshot)
    dataset = UniProtPTMDataset.from_latest()
    dataset.download("/data/uniprot_ptm")

    # Check availability
    versions = get_available_versions()
"""

import csv
import logging
import os
from dataclasses import dataclass
from datetime import datetime
from typing import Dict, List, Optional

from hvantk.core.ptm_constants import (
    UNIPROT_API_URL,
    UNIPROT_API_FIELDS,
    UNIPROT_HUMAN_PTM_QUERY,
    UNIPROT_BATCH_SIZE,
)

logger = logging.getLogger(__name__)

# Output TSV column order
_TSV_COLUMNS = [
    "accession",
    "gene_symbol",
    "position",
    "description",
    "amino_acid",
    "ensembl_xrefs",
    "sequence_length",
]


def _build_search_url(api_url: str, query: str, fields: str, size: int) -> str:
    """Build the initial UniProt search URL with query parameters."""
    import urllib.parse

    params = urllib.parse.urlencode(
        {
            "query": query,
            "fields": fields,
            "format": "json",
            "size": size,
        }
    )
    return f"{api_url}?{params}"


def _get_next_link(response) -> Optional[str]:
    """
    Extract the next page URL from the Link header.

    UniProt uses cursor-based pagination via the HTTP Link header:
        Link: <https://rest.uniprot.org/...?cursor=...>; rel="next"

    Args:
        response: An HTTP response object with a ``headers`` attribute.

    Returns:
        The next-page URL, or None if there is no next page.
    """
    import re

    link_header = response.headers.get("Link", "")
    if not link_header:
        return None

    # Use regex to handle URLs that may contain commas in query params
    match = re.search(r'<([^>]+)>\s*;\s*rel="next"', link_header)
    if match:
        return match.group(1)
    return None


def _extract_ptm_records(entry: dict) -> List[dict]:
    """
    Extract PTM records from a single UniProt JSON entry.

    Each MOD_RES feature in the entry produces one row in the output TSV.

    Args:
        entry: A single UniProt result entry (JSON object).

    Returns:
        List of dicts, one per PTM site, with keys matching ``_TSV_COLUMNS``.
    """
    accession = entry.get("primaryAccession", "")

    # Gene symbol: genes[0].geneName.value
    gene_symbol = ""
    genes = entry.get("genes", [])
    if genes:
        gene_name_obj = genes[0].get("geneName")
        if gene_name_obj:
            gene_symbol = gene_name_obj.get("value", "")

    # Sequence and its length
    sequence_obj = entry.get("sequence", {})
    sequence_value = sequence_obj.get("value", "")
    sequence_length = str(sequence_obj.get("length", len(sequence_value)))

    # Ensembl cross-references
    ensembl_ids = []
    for xref in entry.get("uniProtKBCrossReferences", []):
        if xref.get("database") == "Ensembl":
            xref_id = xref.get("id", "")
            if xref_id:
                ensembl_ids.append(xref_id)
    ensembl_xrefs = ";".join(ensembl_ids)

    # MOD_RES features
    records = []
    for feature in entry.get("features", []):
        if feature.get("type") != "Modified residue":
            continue

        location = feature.get("location", {})
        start_pos = location.get("start", {}).get("value")
        if start_pos is None:
            continue

        position = str(start_pos)

        description = feature.get("description", "")

        # Amino acid at the PTM position (1-based index)
        amino_acid = ""
        if sequence_value and start_pos >= 1 and start_pos <= len(sequence_value):
            amino_acid = sequence_value[start_pos - 1]

        records.append(
            {
                "accession": accession,
                "gene_symbol": gene_symbol,
                "position": position,
                "description": description,
                "amino_acid": amino_acid,
                "ensembl_xrefs": ensembl_xrefs,
                "sequence_length": sequence_length,
            }
        )

    return records


def _fetch_json(url: str) -> tuple:
    """
    Fetch a JSON page from UniProt using *requests* (preferred) or *urllib*.

    Args:
        url: The URL to fetch.

    Returns:
        Tuple of (parsed_json_dict, response_object) where response_object
        exposes ``headers`` for pagination link extraction.
    """
    try:
        import requests as _requests

        resp = _requests.get(url, headers={"Accept": "application/json"}, timeout=120)
        resp.raise_for_status()
        return resp.json(), resp
    except ImportError:
        import json
        import urllib.request

        req = urllib.request.Request(
            url, headers={"Accept": "application/json", "User-Agent": "hvantk/1.0"}
        )
        with urllib.request.urlopen(req, timeout=120) as resp:
            data = json.loads(resp.read().decode("utf-8"))
            return data, resp


@dataclass
class UniProtPTMDataset:
    """
    Represents a UniProt post-translational modification dataset.

    The dataset is obtained by querying the UniProt REST API for all
    reviewed human proteins that carry at least one MOD_RES annotation.
    Results are streamed page-by-page and written to a local TSV file.

    Attributes:
        version_date: Date label for the download (YYYY-MM-DD format)
        api_url: UniProt REST API search endpoint
        query: UniProt search query string
        file_name: Name of the output TSV file
    """

    version_date: str
    api_url: str
    query: str
    file_name: str

    @classmethod
    def from_date(cls, version_date: str) -> "UniProtPTMDataset":
        """
        Create a dataset reference labeled with a specific date.

        Note: UniProt provides a live search endpoint. The date is used
        only for labeling the output file. The downloaded content will
        always reflect the current state of UniProt regardless of the
        date provided.

        Args:
            version_date: Date string in YYYY-MM-DD format

        Returns:
            UniProtPTMDataset instance

        Raises:
            ValueError: If version_date is not in valid format
        """
        try:
            datetime.strptime(version_date, "%Y-%m-%d")
        except ValueError:
            raise ValueError(
                f"Invalid version_date format: {version_date}. Expected YYYY-MM-DD"
            ) from None

        file_name = f"uniprot-ptm-human-{version_date}.tsv"

        return cls(
            version_date=version_date,
            api_url=UNIPROT_API_URL,
            query=UNIPROT_HUMAN_PTM_QUERY,
            file_name=file_name,
        )

    @classmethod
    def from_latest(cls) -> "UniProtPTMDataset":
        """
        Create a dataset reference for the latest available snapshot.

        UniProt provides a live search endpoint that returns current data.
        The version date is set to today's date.

        Returns:
            UniProtPTMDataset instance for today's snapshot
        """
        today = datetime.now().strftime("%Y-%m-%d")
        logger.info(f"Using UniProt live search endpoint, labeling as {today}")
        return cls.from_date(today)

    def download(self, output_dir: str, overwrite: bool = False) -> str:
        """
        Download human PTM data from UniProt and write to a TSV file.

        The method queries the UniProt REST API with cursor-based
        pagination, extracts MOD_RES features from each protein, and
        writes one row per PTM site.

        Args:
            output_dir: Directory to save the output TSV file
            overwrite: If True, overwrite an existing file

        Returns:
            Path to the written TSV file

        Raises:
            FileExistsError: If file exists and overwrite is False
            RuntimeError: If the API request fails
        """
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, self.file_name)

        if os.path.exists(output_path) and not overwrite:
            raise FileExistsError(
                f"File already exists: {output_path}. Use overwrite=True to replace."
            )

        url = _build_search_url(
            self.api_url, self.query, UNIPROT_API_FIELDS, UNIPROT_BATCH_SIZE
        )

        logger.info(f"Downloading UniProt PTM data from {self.api_url}")
        total_records = 0
        page_count = 0

        try:
            with open(output_path, "w", newline="") as fh:
                writer = csv.DictWriter(
                    fh, fieldnames=_TSV_COLUMNS, delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()

                while url:
                    page_count += 1
                    logger.debug(f"Fetching page {page_count}: {url}")

                    data, response = _fetch_json(url)

                    for entry in data.get("results", []):
                        rows = _extract_ptm_records(entry)
                        for row in rows:
                            writer.writerow(row)
                            total_records += 1

                    url = _get_next_link(response)

            logger.info(
                f"Downloaded {total_records} PTM records from {page_count} pages "
                f"to {output_path}"
            )
            return output_path

        except Exception as e:
            # Clean up partial file on failure
            if os.path.exists(output_path):
                os.remove(output_path)
            raise RuntimeError(f"Failed to download UniProt PTM data: {str(e)}") from e

    def get_metadata(self) -> Dict[str, str]:
        """
        Get metadata about this dataset.

        Returns:
            Dictionary with dataset metadata
        """
        return {
            "source": "UniProt Post-Translational Modifications",
            "version_date": self.version_date,
            "api_url": self.api_url,
            "query": self.query,
            "file_name": self.file_name,
            "description": (
                "Post-translational modification (MOD_RES) annotations for "
                "reviewed human proteins from UniProt/Swiss-Prot. Each record "
                "represents a single PTM site with its position, description, "
                "amino acid, and Ensembl cross-references."
            ),
        }

    def __str__(self) -> str:
        return f"UniProtPTMDataset(version={self.version_date})"


def get_available_versions() -> List[str]:
    """
    Check UniProt PTM dataset availability.

    UniProt does not provide versioned archives via the REST API. This
    function verifies the search endpoint is reachable and returns
    today's date as the available version.

    Returns:
        List with today's date if the endpoint is reachable, empty list otherwise

    Note:
        This function requires network access.
    """
    try:
        import urllib.request

        test_url = _build_search_url(
            UNIPROT_API_URL, UNIPROT_HUMAN_PTM_QUERY, UNIPROT_API_FIELDS, 1
        )
        logger.info(f"Checking UniProt search endpoint: {test_url}")

        req = urllib.request.Request(
            test_url,
            method="HEAD",
            headers={"User-Agent": "hvantk/1.0"},
        )

        with urllib.request.urlopen(req, timeout=30) as response:
            if response.status == 200:
                today = datetime.now().strftime("%Y-%m-%d")
                logger.info(
                    "UniProt endpoint is reachable. "
                    "Dataset is generated in real-time (no versioned archives)."
                )
                return [today]

        return []

    except Exception as e:
        logger.error(f"Failed to reach UniProt endpoint: {e}")
        return []


def get_latest_version() -> Optional[str]:
    """
    Get the latest available UniProt PTM dataset version date.

    Since UniProt provides real-time results, this returns today's
    date if the endpoint is reachable.

    Returns:
        Today's date string (YYYY-MM-DD) or None if endpoint is unavailable
    """
    versions = get_available_versions()
    return versions[0] if versions else None
