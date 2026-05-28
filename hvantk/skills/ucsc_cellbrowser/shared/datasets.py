from typing import List, Optional
from dataclasses import dataclass, field
import json
import logging

logger = logging.getLogger(__name__)

from hvantk.skills.ucsc_cellbrowser.shared.constants import (
    UCSC_CELL_BROWSER_BASE_URL,
    EXPRESSION_MATRIX_FILE_NAME,
    METADATA_FILE_NAME,
)
from hvantk.core.utils.file_utils import download_file


@dataclass
class DatasetFacets:
    body_parts: List[str]
    organisms: List[str]
    projects: List[str]
    diseases: List[str]
    life_stages: List[str]
    domains: List[str]
    assays: List[str]
    sources: List[str]


@dataclass
class UCSCDataset:
    """
    Represents a single dataset from the UCSC Cell Browser.

    This class encapsulates metadata about a UCSC cell dataset including its
    identifier, facets (like organisms, diseases, etc.), and provides methods
    to download the associated data files.

    Attributes:
        shortLabel: Short descriptive label for the dataset
        name: Unique identifier for the dataset
        md5: MD5 hash of the dataset
        hasFiles: List of available file types
        body_parts: List of associated body parts
        organisms: List of organisms in the dataset
        tags: List of tags categorizing the dataset
        projects: List of projects associated with the dataset
        diseases: List of diseases studied in the dataset
        life_stages: List of life stages represented
        domains: List of scientific domains
        sources: List of data sources
        assays: List of assay types
        facets: Aggregated facets for filtering/categorization
        sampleCount: Number of samples in the dataset
        isCollection: Whether this is a collection of datasets
        collectionCount: Number of collections (if applicable)
        datasetCount: Number of datasets in this collection (if applicable)
    """

    shortLabel: str
    name: str
    md5: str
    hasFiles: Optional[List[str]] = field(default_factory=list)
    body_parts: Optional[List[str]] = field(default_factory=list)
    organisms: Optional[List[str]] = field(default_factory=list)
    tags: Optional[List[str]] = field(default_factory=list)
    projects: Optional[List[str]] = field(default_factory=list)
    diseases: Optional[List[str]] = field(default_factory=list)
    life_stages: Optional[List[str]] = field(default_factory=list)
    domains: Optional[List[str]] = field(default_factory=list)
    sources: Optional[List[str]] = field(default_factory=list)
    assays: Optional[List[str]] = field(default_factory=list)
    facets: Optional[DatasetFacets] = None
    sampleCount: Optional[int] = None
    isCollection: Optional[bool] = False
    collectionCount: Optional[int] = None
    datasetCount: Optional[int] = None
    children: Optional[List["UCSCDataset"]] = field(default=None, repr=False)

    def __post_init__(self):
        if self.facets is None:
            self.facets = DatasetFacets(
                body_parts=self.body_parts,
                organisms=self.organisms,
                projects=self.projects,
                diseases=self.diseases,
                life_stages=self.life_stages,
                domains=self.domains,
                assays=self.assays,
                sources=self.sources,
            )

    def summary(self) -> str:
        """
        Returns a formatted summary of the dataset.

        Returns:
            str: A human-readable summary of dataset information
        """
        summary_lines = [
            f"Dataset: {self.name} ({self.shortLabel})",
            f"Sample count: {self.sampleCount or 'Unknown'}",
            f"Organisms: {', '.join(self.organisms) if self.organisms else 'None'}",
            f"Body parts: {', '.join(self.body_parts) if self.body_parts else 'None'}",
            f"Diseases: {', '.join(self.diseases) if self.diseases else 'None'}",
            f"Assays: {', '.join(self.assays) if self.assays else 'None'}",
        ]

        return "\n".join(summary_lines)

    def download_expression_matrix(self, out_dir: str) -> str:
        """
        Download the expression matrix file for this dataset.

            Args:
                out_dir: Directory where the file will be saved

            Returns:
                str: Local filesystem path to the downloaded expression matrix file

            Raises:
                ValueError: If the download fails
        """
        url_download = (
            f"{UCSC_CELL_BROWSER_BASE_URL}/{self.name}/{EXPRESSION_MATRIX_FILE_NAME}"
        )
        try:
            return download_file(
                url=url_download, out_dir=out_dir, file_name=EXPRESSION_MATRIX_FILE_NAME
            )
        except Exception as e:
            raise ValueError(
                f"Failed to download expression matrix for {self.name}: {str(e)}"
            ) from e

    def download_metadata(self, out_dir: str) -> str:
        """
        Download the metadata file for this dataset.

        Args:
            out_dir: Directory where the file will be saved

        Returns:
            str: Local filesystem path to the downloaded metadata file

        Raises:
            ValueError: If the download fails
        """
        url_download = f"{UCSC_CELL_BROWSER_BASE_URL}/{self.name}/{METADATA_FILE_NAME}"
        try:
            return download_file(
                url=url_download, out_dir=out_dir, file_name=METADATA_FILE_NAME
            )
        except Exception as e:
            raise ValueError(
                f"Failed to download metadata for {self.name}: {str(e)}"
            ) from e

    def fetch_children(self) -> List["UCSCDataset"]:
        """Fetch child datasets from UCSC API for a collection.

        Makes a single HTTP request to ``{base_url}/{name}/dataset.json``
        and parses the ``datasets`` array.  Results are cached on the
        ``children`` attribute so subsequent calls are free.

        Returns empty list if not a collection or if the fetch fails.
        """
        if not self.isCollection:
            return []

        if self.children is not None:
            return self.children

        import requests

        url = f"{UCSC_CELL_BROWSER_BASE_URL}/{self.name}/dataset.json"
        try:
            resp = requests.get(url, timeout=15)
            resp.raise_for_status()
            data = resp.json()
        except Exception as exc:
            logger.warning("Failed to fetch children for %s: %s", self.name, exc)
            return []

        children = []
        for child in data.get("datasets", []):
            children.append(
                UCSCDataset(
                    shortLabel=child.get("shortLabel", ""),
                    name=child.get("name", ""),
                    md5=child.get("md5", ""),
                    sampleCount=child.get("sampleCount"),
                    isCollection=child.get("isCollection", False),
                    datasetCount=child.get("datasetCount"),
                    body_parts=child.get("body_parts", []),
                    organisms=child.get("organisms", []),
                    diseases=child.get("diseases", []),
                )
            )

        self.children = children
        return children


@dataclass
class UCSCDataSetCollection:
    """
    Represents a collection of UCSC Cell Browser datasets.

    This class manages a group of UCSC datasets, providing methods to
    access, filter, and summarize the collection.

    Attributes:
        datasets: List of UCSCDataset objects in this collection
    """

    datasets: List[UCSCDataset]

    @classmethod
    def from_json(cls, json_path: str) -> "UCSCDataSetCollection":
        """
        Loads a UCSC dataset collection from a JSON file.

        Args:
            json_path: Path to the JSON file containing the datasets.

        Returns:
            An instance of UCSCDataSetCollection populated with datasets from the JSON file.

        Raises:
            ValueError: If the JSON file is invalid or cannot be read.
        """
        logger.info(f"Loading UCSC dataset collection from JSON: {json_path}")
        try:
            with open(json_path, "r") as file:
                data = json.load(file)

            if not isinstance(data, dict) or "datasets" not in data:
                raise ValueError(
                    f"JSON file {json_path} must contain a 'datasets' key with a list"
                )

            datasets = []
            for dataset_dict in data["datasets"]:
                datasets.append(UCSCDataset(**dataset_dict))

            logger.info(f"Loaded {len(datasets)} UCSC datasets from {json_path}")
            return cls(datasets=datasets)

        except FileNotFoundError:
            raise ValueError(f"JSON file not found: {json_path}")
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid JSON in file {json_path}: {str(e)}")
        except Exception as e:
            raise ValueError(f"Error loading datasets from {json_path}: {str(e)}")

    def get_by_name(self, name: str) -> Optional[UCSCDataset]:
        """
        Get a dataset by its name.

        Args:
            name: The dataset name to search for

        Returns:
            The dataset with the matching name, or None if not found
        """
        return next((ds for ds in self.datasets if ds.name == name), None)

    def list_dataset_names(self) -> List[str]:
        """
        List all dataset names in the collection.

        Returns:
            List of names for all datasets in the collection
        """
        return [dataset.name for dataset in self.datasets]

    def search(self, query: str) -> "UCSCDataSetCollection":
        """
        Filter datasets by case-insensitive query across name, label, and facets.

        Args:
            query: Search term to match against dataset fields

        Returns:
            New collection containing only matching datasets
        """
        q = query.lower()
        matches = []
        for ds in self.datasets:
            searchable = " ".join(
                [
                    ds.name,
                    ds.shortLabel,
                    " ".join(ds.body_parts or []),
                    " ".join(ds.organisms or []),
                    " ".join(ds.diseases or []),
                ]
            ).lower()
            if q in searchable:
                matches.append(ds)
        return UCSCDataSetCollection(datasets=matches)

    def filter_by_organism(self, organism: str) -> "UCSCDataSetCollection":
        """
        Filter datasets by organism.

        Args:
            organism: Organism to filter by

        Returns:
            New collection containing only datasets with the specified organism
        """
        filtered_datasets = [
            ds for ds in self.datasets if organism in (ds.organisms or [])
        ]
        return UCSCDataSetCollection(datasets=filtered_datasets)

    def summary(self) -> str:
        """
        Returns a formatted summary of the collection.

        Returns:
            str: A human-readable summary of the collection
        """
        if not self.datasets:
            return "Empty UCSC dataset collection"

        organism_counts = {}
        for dataset in self.datasets:
            for organism in dataset.organisms or []:
                organism_counts[organism] = organism_counts.get(organism, 0) + 1

        summary_lines = [
            f"UCSC Dataset Collection",
            f"Total datasets: {len(self.datasets)}",
            "Organisms:",
        ]

        for organism, count in organism_counts.items():
            summary_lines.append(f"  {organism}: {count}")

        return "\n".join(summary_lines)


def load_ucsc_datasets(json_path: Optional[str] = None) -> List[UCSCDataset]:
    """
    Load UCSC datasets from JSON file.

    Args:
        json_path: Path to JSON file. If None, uses default resource file.

    Returns:
        List of UCSCDataset objects
    """
    if json_path is None:
        import os

        json_path = os.path.join(
            os.path.dirname(__file__), "..", "resources", "cells_ucsc_datasets.json"
        )

    collection = UCSCDataSetCollection.from_json(json_path)
    return collection.datasets
