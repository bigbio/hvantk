"""
New unified registry manager for the reorganized hvantk registry system.
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Union
from datetime import datetime

# Import schema validator from same directory
try:
    from schema_validator import SchemaValidator
except ImportError:
    import sys

    sys.path.append(str(Path(__file__).parent))
    from schema_validator import SchemaValidator

logger = logging.getLogger(__name__)


class RegistryManager:
    """Unified manager for accessing and managing omics datasets across all categories."""

    def __init__(self, registry_root: Optional[Path] = None):
        self.registry_root = registry_root or Path(__file__).parent / "registry"
        self.validator = SchemaValidator()
        self.omics_types = ["transcriptomics", "proteomics", "genomics", "epigenomics"]
        self._datasets_cache = {}
        self._load_all_datasets()

    def _load_all_datasets(self):
        """Load all datasets from all omics categories."""
        for omics_type in self.omics_types:
            datasets_file = self.registry_root / omics_type / "datasets.json"
            if datasets_file.exists():
                try:
                    with open(datasets_file, "r") as f:
                        datasets = json.load(f)
                    self._datasets_cache[omics_type] = datasets
                    logger.info(f"Loaded {len(datasets)} {omics_type} datasets")
                except Exception as e:
                    logger.error(f"Error loading {omics_type} datasets: {e}")
                    self._datasets_cache[omics_type] = []
            else:
                self._datasets_cache[omics_type] = []

    def get_all_datasets(self) -> Dict[str, List[Dict]]:
        """Get all datasets organized by omics type."""
        return self._datasets_cache.copy()

    def get_datasets_by_type(self, omics_type: str) -> List[Dict]:
        """Get all datasets for a specific omics type."""
        return self._datasets_cache.get(omics_type, [])

    def search_datasets(
        self,
        query: str = "",
        omics_type: Optional[str] = None,
        organism: Optional[str] = None,
        data_source: Optional[str] = None,
    ) -> List[Dict]:
        """
        Search datasets across all omics types or within a specific type.

        Args:
            query: Text to search in title, description, and accession
            omics_type: Filter by specific omics type
            organism: Filter by organism
            data_source: Filter by data source
        """
        results = []

        search_types = [omics_type] if omics_type else self.omics_types

        for otype in search_types:
            datasets = self._datasets_cache.get(otype, [])

            for dataset in datasets:
                # Text search
                if query:
                    searchable_text = " ".join(
                        [
                            dataset.get("title", ""),
                            dataset.get("description", ""),
                            dataset.get("accession", ""),
                        ]
                    ).lower()

                    if query.lower() not in searchable_text:
                        continue

                # Filter criteria
                if organism and dataset.get("organism", "").lower() != organism.lower():
                    continue

                if (
                    data_source
                    and dataset.get("data_source", "").lower() != data_source.lower()
                ):
                    continue

                # Add omics type annotation
                result = dataset.copy()
                result["_omics_type"] = otype
                results.append(result)

        return results

    def get_dataset_by_accession(self, accession: str) -> Optional[Dict]:
        """Find a dataset by its accession number across all omics types."""
        for omics_type in self.omics_types:
            datasets = self._datasets_cache.get(omics_type, [])
            for dataset in datasets:
                if dataset.get("accession") == accession:
                    result = dataset.copy()
                    result["_omics_type"] = omics_type
                    return result
        return None

    def add_dataset(self, dataset: Dict, omics_type: str) -> bool:
        """
        Add a new dataset to the registry.

        Args:
            dataset: Dataset metadata dictionary
            omics_type: Which omics category to add to

        Returns:
            True if successful, False otherwise
        """
        if omics_type not in self.omics_types:
            logger.error(f"Invalid omics type: {omics_type}")
            return False

        # Validate against schema
        schema_name = f"{omics_type}_schema"
        is_valid, errors = self.validator.validate_dataset(dataset, schema_name)

        if not is_valid:
            logger.error(f"Dataset validation failed: {errors}")
            return False

        # Add timestamp
        dataset["last_updated"] = datetime.now().isoformat()

        # Add to cache
        if omics_type not in self._datasets_cache:
            self._datasets_cache[omics_type] = []

        self._datasets_cache[omics_type].append(dataset)

        # Save to file
        return self._save_datasets(omics_type)

    def update_dataset(self, accession: str, updates: Dict) -> bool:
        """Update an existing dataset by accession."""
        # Find the dataset
        dataset_info = self.get_dataset_by_accession(accession)
        if not dataset_info:
            logger.error(f"Dataset not found: {accession}")
            return False

        omics_type = dataset_info["_omics_type"]
        datasets = self._datasets_cache[omics_type]

        # Find and update
        for i, dataset in enumerate(datasets):
            if dataset.get("accession") == accession:
                # Apply updates
                dataset.update(updates)
                dataset["last_updated"] = datetime.now().isoformat()

                # Validate updated dataset
                schema_name = f"{omics_type}_schema"
                is_valid, errors = self.validator.validate_dataset(dataset, schema_name)

                if not is_valid:
                    logger.error(f"Updated dataset validation failed: {errors}")
                    return False

                # Save to file
                return self._save_datasets(omics_type)

        return False

    def _save_datasets(self, omics_type: str) -> bool:
        """Save datasets for a specific omics type to file."""
        try:
            datasets_file = self.registry_root / omics_type / "datasets.json"
            datasets_file.parent.mkdir(parents=True, exist_ok=True)

            with open(datasets_file, "w") as f:
                json.dump(self._datasets_cache[omics_type], f, indent=2)

            logger.info(
                f"Saved {len(self._datasets_cache[omics_type])} {omics_type} datasets"
            )
            return True
        except Exception as e:
            logger.error(f"Error saving {omics_type} datasets: {e}")
            return False

    def get_registry_stats(self) -> Dict:
        """Get statistics about the registry."""
        stats = {
            "total_datasets": 0,
            "by_omics_type": {},
            "by_organism": {},
            "by_data_source": {},
        }

        for omics_type, datasets in self._datasets_cache.items():
            stats["by_omics_type"][omics_type] = len(datasets)
            stats["total_datasets"] += len(datasets)

            for dataset in datasets:
                # Count by organism
                organism = dataset.get("organism", "Unknown")
                stats["by_organism"][organism] = (
                    stats["by_organism"].get(organism, 0) + 1
                )

                # Count by data source
                data_source = dataset.get("data_source", "Unknown")
                stats["by_data_source"][data_source] = (
                    stats["by_data_source"].get(data_source, 0) + 1
                )

        return stats

    def reload(self):
        """Reload all datasets from files."""
        self._datasets_cache = {}
        self._load_all_datasets()
