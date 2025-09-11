import json
import logging
from dataclasses import dataclass, field
from typing import List, Optional, Dict

logger = logging.getLogger(__name__)

from hvantk.data.file_utils import download_file
from hvantk.core.constants import EXPRESSION_ATLAS_BASE_URL


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
            file_type: Type of the file to download (e.g., "transcript-tpm", "sdrf")
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

    def download_expression_data(self, out_dir: str) -> str:
        """
        Download expression data (transcript TPM) for this dataset.

        Args:
            out_dir: Directory where the file will be saved

        Returns:
            Path to the downloaded file

        Raises:
            ValueError: If the transcript-tpm file is not available or download fails
        """
        return self.download_file("transcript-tpm", out_dir)

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

    def download_expression_file(self, out_dir: str, file_name: str) -> str:
        """Download a specific expression file by name."""
        url_download = f"{EXPRESSION_ATLAS_BASE_URL}/{self.accession}/download/{file_name}"
        try:
            return download_file(url=url_download, out_dir=out_dir, file_name=file_name)
        except Exception as e:
            raise ValueError(f"Failed to download {file_name} for {self.accession}: {str(e)}") from e

    def download_sdrf_file(self, out_dir: str, file_name: str) -> str:
        """Download a specific SDRF file by name."""
        url_download = f"{EXPRESSION_ATLAS_BASE_URL}/{self.accession}/download/{file_name}"
        try:
            return download_file(url=url_download, out_dir=out_dir, file_name=file_name)
        except Exception as e:
            raise ValueError(f"Failed to download {file_name} for {self.accession}: {str(e)}") from e


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
                raise ValueError(f"JSON file {json_path} must contain a list of datasets")

            datasets = []
            for dataset_dict in data:
                datasets.append(ExpressionAtlasDataset(**dataset_dict))

            logger.info(f"Loaded {len(datasets)} Expression Atlas datasets from {json_path}")
            return cls(datasets=datasets)

        except FileNotFoundError:
            raise ValueError(f"JSON file not found: {json_path}")
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid JSON in file {json_path}: {str(e)}")
        except Exception as e:
            raise ValueError(f"Error loading datasets from {json_path}: {str(e)}")

    def filter_by_type(self, dataset_type: str) -> "ExpressionAtlasDatasetCollection":
        """
        Filter datasets by type.

        Args:
            dataset_type: Type to filter by (e.g., "RNA-Seq mRNA baseline")

        Returns:
            New collection containing only datasets of the specified type
        """
        filtered_datasets = [ds for ds in self.datasets if ds.type == dataset_type]
        return ExpressionAtlasDatasetCollection(datasets=filtered_datasets)

    def get_by_accession(self, accession: str) -> Optional[ExpressionAtlasDataset]:
        """
        Get a dataset by its accession ID.

        Args:
            accession: The accession ID to search for

        Returns:
            The dataset with the matching accession, or None if not found
        """
        return next((ds for ds in self.datasets if ds.accession == accession), None)

    def list_dataset_accessions(self) -> List[str]:
        """
        List all dataset accessions in the collection.

        Returns:
            List of accession IDs for all datasets in the collection
        """
        return [dataset.accession for dataset in self.datasets]

    def summary(self) -> str:
        """
        Returns a formatted summary of the collection.

        Returns:
            str: A human-readable summary of the collection
        """
        if not self.datasets:
            return "Empty Expression Atlas dataset collection"

        type_counts = {}
        for dataset in self.datasets:
            type_counts[dataset.type] = type_counts.get(dataset.type, 0) + 1

        summary_lines = [
            f"Expression Atlas Dataset Collection",
            f"Total datasets: {len(self.datasets)}",
            "Dataset types:",
        ]

        for dtype, count in type_counts.items():
            summary_lines.append(f"  {dtype}: {count}")

        return "\n".join(summary_lines)


def load_expression_atlas_datasets(json_path: Optional[str] = None) -> List[ExpressionAtlasDataset]:
    """
    Load Expression Atlas datasets from JSON file.

    Args:
        json_path: Path to JSON file. If None, uses default resource file.

    Returns:
        List of ExpressionAtlasDataset objects
    """
    if json_path is None:
        import os
        json_path = os.path.join(
            os.path.dirname(__file__),
            "..", "resources", "expression_atlas.json"
        )

    collection = ExpressionAtlasDatasetCollection.from_json(json_path)
    return collection.datasets
