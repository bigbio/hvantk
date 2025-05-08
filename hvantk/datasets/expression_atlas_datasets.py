import json
import logging
from dataclasses import dataclass, field
from typing import List, Optional, Dict

logger = logging.getLogger(__name__)

from hvantk.utils.file_utils import download_file
from hvantk.utils.constants import EXPRESSION_ATLAS_BASE_URL


@dataclass
class ExpressionAtlasDataset:
    """
    Represents a single dataset from the Expression Atlas.

    This class encapsulates metadata about an Expression Atlas dataset including its
    title, accession ID, type, and provides methods to download the associated data files.

    Attributes:
        title: Descriptive title of the dataset
        accession: Unique accession ID for the dataset
        type: Type of the dataset (e.g., "RNA-Seq mRNA baseline")
        pubmedid: PubMed ID for the associated publication
        description: Detailed description of the dataset
        files: List of file information associated with the dataset
    """

    title: str
    accession: str
    type: str
    pubmedid: Optional[str] = None
    description: Optional[str] = ""
    files: Optional[List[Dict[str, str]]] = field(default_factory=list)

    def summary(self) -> str:
        """
        Returns a formatted summary of the dataset.

        Returns:
            str: A human-readable summary of dataset information
        """
        summary_lines = [
            f"Dataset: {self.accession} ({self.title})",
            f"Type: {self.type}",
            f"PubMed ID: {self.pubmedid or 'Not available'}",
            f"Description: {self.description or 'Not available'}",
            f"Files: {', '.join([f['type'] for f in self.files]) if self.files else 'None'}",
        ]

        return "\n".join(summary_lines)

    def download_file(self, file_type: str, out_dir: str) -> str:
        """
        Download a specific file for this dataset by file type.

        Args:
            file_type: Type of the file to download (e.g., "tpm", "fpkm", "sdrf")
            out_dir: Directory where the file will be saved

        Returns:
            Path to the downloaded file

        Raises:
            ValueError: If the file type is not available or download fails
        """
        file_info = next((f for f in self.files if f["type"] == file_type), None)
        if not file_info:
            raise ValueError(f"File type '{file_type}' not available for {self.accession}")

        file_name = file_info["name"]
        url_download = f"{EXPRESSION_ATLAS_BASE_URL}/{self.accession}/download/{file_name}"
        
        try:
            return download_file(
                url=url_download, out_dir=out_dir, file_name=file_name
            )
        except Exception as e:
            raise ValueError(
                f"Failed to download {file_type} file for {self.accession}: {str(e)}"
            ) from e

    def download_expression_data(self, out_dir: str, format: str = "tpm") -> str:
        """
        Download expression data for this dataset.

        Args:
            out_dir: Directory where the file will be saved
            format: Format of the expression data to download ("tpm" or "fpkm")

        Returns:
            Path to the downloaded file

        Raises:
            ValueError: If the specified format is not available or download fails
        """
        if format not in ["tpm", "fpkm"]:
            raise ValueError(f"Invalid format '{format}'. Must be 'tpm' or 'fpkm'")
        
        return self.download_file(format, out_dir)

    def download_metadata(self, out_dir: str) -> str:
        """
        Download the metadata file (sdrf) for this dataset.

        Args:
            out_dir: Directory where the file will be saved

        Returns:
            Path to the downloaded file

        Raises:
            ValueError: If the metadata file is not available or download fails
        """
        return self.download_file("sdrf", out_dir)


@dataclass
class ExpressionAtlasDatasetCollection:
    """
    Represents a collection of Expression Atlas datasets.

    This class manages a group of Expression Atlas datasets, providing methods to
    access, filter, and summarize the collection. It can be instantiated from
    a JSON file containing dataset metadata.

    Attributes:
        datasets: List of ExpressionAtlasDataset objects in this collection
    """

    datasets: List[ExpressionAtlasDataset]

    @classmethod
    def from_json(cls, json_path: str) -> "ExpressionAtlasDatasetCollection":
        """
        Loads an Expression Atlas dataset collection from a JSON file.

        Reads the specified JSON file and constructs an ExpressionAtlasDatasetCollection 
        instance with its datasets.

        Args:
            json_path: Path to the JSON file containing the datasets.

        Returns:
            An instance of ExpressionAtlasDatasetCollection populated with datasets from the JSON file.

        Raises:
            ValueError: If the JSON file is invalid or cannot be read.
        """
        logger.info(f"Loading Expression Atlas dataset collection from JSON: {json_path}")
        try:
            with open(json_path, "r") as file:
                data = json.load(file)

            if not isinstance(data, list):
                raise ValueError("The JSON file should contain a list of datasets")

            if not data:
                raise ValueError("The dataset list cannot be empty")

            datasets = [ExpressionAtlasDataset(**ds) for ds in data]
            return cls(datasets=datasets)
            
        except json.JSONDecodeError as e:
            logger.exception(f"Invalid JSON format: {str(e)}")
            raise ValueError(f"Invalid JSON format: {str(e)}") from e
        except OSError as e:
            logger.exception(f"Could not read file {json_path}: {str(e)}")
            raise ValueError(f"Could not read file {json_path}: {str(e)}") from e
        except (KeyError, TypeError) as e:
            logger.exception(f"Invalid dataset format in JSON: {str(e)}")
            raise ValueError(f"Invalid dataset format in JSON: {str(e)}") from e

    def get_dataset_by_accession(self, accession: str) -> Optional[ExpressionAtlasDataset]:
        """
        Returns the dataset with the specified accession ID, or None if not found.

        Args:
            accession: The accession ID of the dataset to search for.

        Returns:
            The ExpressionAtlasDataset instance matching the given accession, or None if no match exists.
        """
        logger.debug(f"Getting dataset by accession: {accession}")
        for dataset in self.datasets:
            if dataset.accession == accession:
                logger.debug(f"Dataset found: {dataset.accession}")
                return dataset
        logger.debug(f"Dataset not found: {accession}")
        return None

    def list_dataset_accessions(self) -> List[str]:
        """
        Returns a list of all dataset accession IDs in the collection.

        Returns:
            List[str]: The accession IDs of all datasets.
        """
        logger.debug("Listing dataset accessions")
        accessions = [dataset.accession for dataset in self.datasets]
        logger.debug(f"Dataset accessions: {accessions}")
        return accessions

    def filter_by_type(self, dataset_type: str) -> List[ExpressionAtlasDataset]:
        """
        Filters the collection to include only datasets of the specified type.

        Args:
            dataset_type: The type of datasets to include

        Returns:
            List[ExpressionAtlasDataset]: Filtered list of datasets
        """
        return [dataset for dataset in self.datasets if dataset.type == dataset_type]

    def filter_by_organism(self, organism: str) -> List[ExpressionAtlasDataset]:
        """
        Filters the collection to include only datasets related to the specified organism.

        Args:
            organism: Name of the organism to filter by (case-insensitive partial match)

        Returns:
            List[ExpressionAtlasDataset]: Filtered list of datasets
        """
        organism_lower = organism.lower()
        return [
            dataset for dataset in self.datasets 
            if organism_lower in dataset.title.lower() or organism_lower in (dataset.description or "").lower()
        ]