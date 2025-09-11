"""
Dataset Creation Registry for tracking datasets that can successfully create Hail MatrixTables.

This module provides a comprehensive registry of all datasets that have been validated
for Hail MatrixTable creation, including their metadata, creation requirements,
and recommended usage patterns.
"""
from __future__ import annotations

import json
import logging
from dataclasses import dataclass, asdict
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Dict, List, Optional, Any, Union

logger = logging.getLogger(__name__)


class DatasetSource(Enum):
    """Dataset source types."""
    UCSC = "ucsc"
    EXPRESSION_ATLAS = "expression_atlas"
    CPTAC = "cptac"
    CUSTOM = "custom"


class MatrixType(Enum):
    """Types of matrices that can be created."""
    EXPRESSION = "expression"
    METHYLATION = "methylation"
    GENOTYPE = "genotype"
    PROTEOMICS = "proteomics"
    SINGLE_CELL = "single_cell"
    BULK_RNA = "bulk_rna"


class CreationStatus(Enum):
    """Status of dataset matrix creation capability."""
    VALIDATED = "validated"           # Successfully creates matrices
    FAILED = "failed"                # Cannot create matrices
    UNTESTED = "untested"            # Not yet validated
    DEPRECATED = "deprecated"        # No longer supported
    REQUIRES_UPDATE = "requires_update"  # Needs configuration update


@dataclass
class DatasetMetadata:
    """Metadata for a dataset."""
    title: str
    description: str
    organism: str
    tissue_type: Optional[str] = None
    sample_count: Optional[int] = None
    feature_count: Optional[int] = None
    publication: Optional[str] = None
    doi: Optional[str] = None
    release_date: Optional[str] = None
    file_size_mb: Optional[float] = None
    tags: List[str] = None

    def __post_init__(self):
        if self.tags is None:
            self.tags = []


@dataclass
class CreationRequirements:
    """Requirements for creating matrices from this dataset."""
    min_memory_gb: int = 4
    min_partitions: int = 1
    max_partitions: Optional[int] = None
    required_tools: List[str] = None
    preprocessing_steps: List[str] = None
    special_handling: Optional[str] = None

    def __post_init__(self):
        if self.required_tools is None:
            self.required_tools = ["hail"]
        if self.preprocessing_steps is None:
            self.preprocessing_steps = []


@dataclass
class FileInfo:
    """Information about dataset files."""
    expression_matrix: Optional[str] = None
    metadata: Optional[str] = None
    features: Optional[str] = None
    additional_files: Dict[str, str] = None

    def __post_init__(self):
        if self.additional_files is None:
            self.additional_files = {}


@dataclass
class UsageExample:
    """Example of how to use this dataset."""
    title: str
    description: str
    code_snippet: str
    expected_output: Optional[str] = None


@dataclass
class DatasetEntry:
    """Complete dataset registry entry."""
    dataset_id: str
    source: DatasetSource
    matrix_type: MatrixType
    status: CreationStatus
    metadata: DatasetMetadata
    requirements: CreationRequirements
    files: FileInfo
    creation_command: str
    validation_date: str
    last_tested: str
    usage_examples: List[UsageExample] = None
    known_issues: List[str] = None
    alternative_approaches: List[str] = None

    def __post_init__(self):
        if self.usage_examples is None:
            self.usage_examples = []
        if self.known_issues is None:
            self.known_issues = []
        if self.alternative_approaches is None:
            self.alternative_approaches = []


class DatasetCreationRegistry:
    """Registry for tracking dataset creation capabilities."""

    def __init__(self, registry_file: Optional[str] = None):
        """Initialize registry."""
        self.registry_file = registry_file or self._get_default_registry_file()
        self.datasets: Dict[str, DatasetEntry] = {}
        self.load_registry()

    def _get_default_registry_file(self) -> str:
        """Get default registry file path."""
        return str(Path(__file__).parent.parent / "resources" / "dataset_creation_registry.json")

    def load_registry(self) -> None:
        """Load registry from file."""
        try:
            if Path(self.registry_file).exists():
                with open(self.registry_file, 'r') as f:
                    data = json.load(f)
                    self.datasets = {
                        dataset_id: self._dict_to_dataset_entry(entry_data)
                        for dataset_id, entry_data in data.items()
                    }
                logger.info(f"Loaded {len(self.datasets)} datasets from registry")
            else:
                logger.info("Registry file not found, starting with empty registry")
                self.datasets = {}
        except Exception as e:
            logger.error(f"Failed to load registry: {e}")
            self.datasets = {}

    def save_registry(self) -> None:
        """Save registry to file."""
        try:
            # Ensure directory exists
            Path(self.registry_file).parent.mkdir(parents=True, exist_ok=True)

            data = {
                dataset_id: self._dataset_entry_to_dict(entry)
                for dataset_id, entry in self.datasets.items()
            }

            with open(self.registry_file, 'w') as f:
                json.dump(data, f, indent=2, default=str)
            logger.info(f"Saved {len(self.datasets)} datasets to registry")
        except Exception as e:
            logger.error(f"Failed to save registry: {e}")
            raise

    def _dict_to_dataset_entry(self, data: Dict[str, Any]) -> DatasetEntry:
        """Convert dictionary to DatasetEntry."""
        # Convert nested dictionaries to dataclasses
        metadata = DatasetMetadata(**data['metadata'])
        requirements = CreationRequirements(**data['requirements'])
        files = FileInfo(**data['files'])

        usage_examples = [
            UsageExample(**example) for example in data.get('usage_examples', [])
        ]

        return DatasetEntry(
            dataset_id=data['dataset_id'],
            source=DatasetSource(data['source']),
            matrix_type=MatrixType(data['matrix_type']),
            status=CreationStatus(data['status']),
            metadata=metadata,
            requirements=requirements,
            files=files,
            creation_command=data['creation_command'],
            validation_date=data['validation_date'],
            last_tested=data['last_tested'],
            usage_examples=usage_examples,
            known_issues=data.get('known_issues', []),
            alternative_approaches=data.get('alternative_approaches', [])
        )

    def _dataset_entry_to_dict(self, entry: DatasetEntry) -> Dict[str, Any]:
        """Convert DatasetEntry to dictionary."""
        return {
            'dataset_id': entry.dataset_id,
            'source': entry.source.value,
            'matrix_type': entry.matrix_type.value,
            'status': entry.status.value,
            'metadata': asdict(entry.metadata),
            'requirements': asdict(entry.requirements),
            'files': asdict(entry.files),
            'creation_command': entry.creation_command,
            'validation_date': entry.validation_date,
            'last_tested': entry.last_tested,
            'usage_examples': [asdict(example) for example in entry.usage_examples],
            'known_issues': entry.known_issues,
            'alternative_approaches': entry.alternative_approaches
        }

    def add_dataset(self, dataset: DatasetEntry) -> None:
        """Add a dataset to the registry."""
        self.datasets[dataset.dataset_id] = dataset
        logger.info(f"Added dataset {dataset.dataset_id} to registry")

    def get_dataset(self, dataset_id: str) -> Optional[DatasetEntry]:
        """Get a dataset by ID."""
        return self.datasets.get(dataset_id)

    def list_datasets(self,
                     source: Optional[DatasetSource] = None,
                     matrix_type: Optional[MatrixType] = None,
                     status: Optional[CreationStatus] = None) -> List[DatasetEntry]:
        """List datasets with optional filtering."""
        datasets = list(self.datasets.values())

        if source:
            datasets = [d for d in datasets if d.source == source]
        if matrix_type:
            datasets = [d for d in datasets if d.matrix_type == matrix_type]
        if status:
            datasets = [d for d in datasets if d.status == status]

        return datasets

    def get_validated_datasets(self) -> List[DatasetEntry]:
        """Get all validated datasets that can create matrices."""
        return self.list_datasets(status=CreationStatus.VALIDATED)

    def get_failed_datasets(self) -> List[DatasetEntry]:
        """Get all failed datasets."""
        return self.list_datasets(status=CreationStatus.FAILED)

    def get_summary_stats(self) -> Dict[str, Any]:
        """Get summary statistics."""
        total = len(self.datasets)
        by_status = {}
        by_source = {}
        by_matrix_type = {}

        for dataset in self.datasets.values():
            # Count by status
            status_key = dataset.status.value
            by_status[status_key] = by_status.get(status_key, 0) + 1

            # Count by source
            source_key = dataset.source.value
            by_source[source_key] = by_source.get(source_key, 0) + 1

            # Count by matrix type
            matrix_key = dataset.matrix_type.value
            by_matrix_type[matrix_key] = by_matrix_type.get(matrix_key, 0) + 1

        return {
            'total': total,
            'by_status': by_status,
            'by_source': by_source,
            'by_matrix_type': by_matrix_type,
            'validated_count': by_status.get('validated', 0),
            'failed_count': by_status.get('failed', 0)
        }

    def generate_creation_report(self, output_file: Optional[str] = None) -> str:
        """Generate a comprehensive creation report."""
        stats = self.get_summary_stats()
        validated_datasets = self.get_validated_datasets()
        failed_datasets = self.get_failed_datasets()

        report_lines = [
            "Dataset Creation Registry Report",
            "=" * 40,
            f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            "",
            "Summary Statistics:",
            f"  Total Datasets: {stats['total']}",
            f"  Validated (Ready for Use): {stats['validated_count']}",
            f"  Failed: {stats['failed_count']}",
            "",
            "By Source:",
        ]

        for source, count in stats['by_source'].items():
            report_lines.append(f"  {source}: {count}")

        report_lines.extend([
            "",
            "By Matrix Type:",
        ])

        for matrix_type, count in stats['by_matrix_type'].items():
            report_lines.append(f"  {matrix_type}: {count}")

        if validated_datasets:
            report_lines.extend([
                "",
                "✅ VALIDATED DATASETS (Ready for Use):",
                "-" * 50,
            ])

            for dataset in validated_datasets:
                report_lines.extend([
                    f"Dataset ID: {dataset.dataset_id}",
                    f"  Title: {dataset.metadata.title}",
                    f"  Source: {dataset.source.value}",
                    f"  Type: {dataset.matrix_type.value}",
                    f"  Samples: {dataset.metadata.sample_count or 'Unknown'}",
                    f"  Creation Command: {dataset.creation_command}",
                    f"  Last Tested: {dataset.last_tested}",
                    ""
                ])

        if failed_datasets:
            report_lines.extend([
                "",
                "❌ FAILED DATASETS:",
                "-" * 50,
            ])

            for dataset in failed_datasets:
                report_lines.extend([
                    f"Dataset ID: {dataset.dataset_id}",
                    f"  Title: {dataset.metadata.title}",
                    f"  Source: {dataset.source.value}",
                    f"  Known Issues: {'; '.join(dataset.known_issues) if dataset.known_issues else 'None listed'}",
                    ""
                ])

        report = "\n".join(report_lines)

        if output_file:
            with open(output_file, 'w') as f:
                f.write(report)
            logger.info(f"Report saved to {output_file}")

        return report

    def get_creation_examples(self, dataset_id: str) -> List[str]:
        """Get creation command examples for a dataset."""
        dataset = self.get_dataset(dataset_id)
        if not dataset:
            return []

        examples = [
            f"# Create matrix from {dataset.dataset_id}",
            dataset.creation_command,
            ""
        ]

        # Add usage examples
        for example in dataset.usage_examples:
            examples.extend([
                f"# {example.title}",
                f"# {example.description}",
                example.code_snippet,
                ""
            ])

        return examples
