from typing import List, Optional
from dataclasses import dataclass, field
import json
import logging

logger = logging.getLogger(__name__)

from hvantk.utils.constants import (
    UCSC_CELL_BROWSER_BASE_URL,
    EXPRESSION_MATRIX_FILE_NAME,
    METADATA_FILE_NAME,
)
from hvantk.utils.file_utils import download_file


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
                Path to the downloaded file

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
            Path to the downloaded file

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


@dataclass
class UCSCDataSetCollection:
    """
    Represents a collection of UCSC cell datasets.

     This class manages a group of related UCSC datasets, providing methods to
     access, filter, and summarize the collection. It can be instantiated from
     a JSON file containing dataset metadata.

     Attributes:
         shortLabel: Short descriptive label for the collection
         abstract: Description of the dataset collection
         inDir: Input directory information
         name: Unique identifier for the collection
         datasets: List of UCSCDataset objects in this collection
    """

    shortLabel: str
    abstract: str
    inDir: str
    name: str
    datasets: List[UCSCDataset]

    @classmethod
    def from_json(cls, json_path: str) -> "UCSCDataSetCollection":
        """
        Loads a UCSC dataset collection and its datasets from a JSON file.
        
        Reads the specified JSON file, validates required fields, and constructs a UCSCDataSetCollection instance with its datasets. Raises a ValueError if the file is missing required keys, contains invalid data, or cannot be read.
         
        Args:
            json_path: Path to the JSON file containing the collection metadata and datasets.
        
        Returns:
            An instance of UCSCDataSetCollection populated with datasets from the JSON file.
        
        Raises:
            ValueError: If the JSON file is invalid, missing required fields, or cannot be read.
        """
        logger.info(f"Loading UCSC dataset collection from JSON: {json_path}")
        try:
            with open(json_path, "r") as file:
                data = json.load(file)

            required_keys = ["shortLabel", "abstract", "inDir", "name", "datasets"]
            missing_keys = [key for key in required_keys if key not in data]
            if missing_keys:
                raise ValueError(
                    f"Missing required keys in JSON: {', '.join(missing_keys)}"
                )

            if not isinstance(data["datasets"], list):
                raise ValueError("The 'datasets' field must be a list")

            if not data["datasets"]:
                raise ValueError("The 'datasets' list cannot be empty")

            datasets = [UCSCDataset(**ds) for ds in data["datasets"]]
            return cls(
                shortLabel=data["shortLabel"],
                abstract=data["abstract"],
                inDir=data["inDir"],
                name=data["name"],
                datasets=datasets,
            )
        except json.JSONDecodeError as e:
            logger.exception(f"Invalid JSON format: {str(e)}")
            raise ValueError(f"Invalid JSON format: {str(e)}") from e
        except OSError as e:
            logger.exception(f"Could not read file {json_path}: {str(e)}")
            raise ValueError(f"Could not read file {json_path}: {str(e)}") from e
        except (KeyError, TypeError) as e:
            logger.exception(f"Invalid dataset format in JSON: {str(e)}")
            raise ValueError(f"Invalid dataset format in JSON: {str(e)}") from e

    def get_dataset_by_name(self, dataset_name: str) -> Optional[UCSCDataset]:
        """
        Returns the dataset with the specified name, or None if not found.
        
        Args:
        	dataset_name: The name of the dataset to search for.
        
        Returns:
        	The UCSCDataset instance matching the given name, or None if no match exists.
        """
        logger.debug(f"Getting dataset by name: {dataset_name}")
        """
        Retrieve a dataset by its name.

        Args:
            dataset_name: The unique name identifier of the dataset to find

        Returns:
            The matching UCSCDataset object or None if not found
        """
        for dataset in self.datasets:
            if dataset.name == dataset_name:
                logger.debug(f"Dataset found: {dataset.name}")
                return dataset
        logger.debug(f"Dataset not found: {dataset_name}")
        return None

    def total_samples(self) -> int:
        """
        Returns the total number of samples across all datasets in the collection.
        
        If a dataset's sample count is missing, it is treated as zero.
        """
        return sum(dataset.sampleCount or 0 for dataset in self.datasets)

    def list_dataset_names(self) -> List[str]:
        """
        Returns a list of all dataset names in the collection.
        
        Returns:
            List[str]: The names of all datasets.
        """
        logger.debug("Listing dataset names")
        """
        List the names of all datasets in the collection.

        Returns:
            List[str]: A list of dataset names
        """
        dataset_names = [dataset.name for dataset in self.datasets]
        logger.debug(f"Dataset names: {dataset_names}")
        return dataset_names
