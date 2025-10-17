"""
Unified registry interface for hvantk - replaces multiple separate registry files.
"""
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Any
from datetime import datetime

logger = logging.getLogger(__name__)

class HvantkRegistry:
    """Main interface for accessing all omics datasets in hvantk."""

    def __init__(self, registry_root: Optional[Path] = None):
        """Initialize the unified registry."""
        self.registry_root = registry_root or Path(__file__).parent / "registry"
        self.omics_types = ["transcriptomics", "proteomics", "genomics", "epigenomics"]
        self._cache = {}
        self._load_registry()

    def _load_registry(self):
        """Load all registry data into memory cache."""
        for omics_type in self.omics_types:
            datasets_file = self.registry_root / omics_type / "datasets.json"
            if datasets_file.exists():
                try:
                    with open(datasets_file, 'r') as f:
                        self._cache[omics_type] = json.load(f)
                except Exception as e:
                    logger.error(f"Error loading {omics_type}: {e}")
                    self._cache[omics_type] = []
            else:
                self._cache[omics_type] = []

    def list_transcriptomics_datasets(self) -> List[Dict]:
        """Get all transcriptomics datasets."""
        return self._cache.get("transcriptomics", [])

    def list_proteomics_datasets(self) -> List[Dict]:
        """Get all proteomics datasets."""
        return self._cache.get("proteomics", [])

    def list_genomics_datasets(self) -> List[Dict]:
        """Get all genomics datasets."""
        return self._cache.get("genomics", [])

    def list_epigenomics_datasets(self) -> List[Dict]:
        """Get all epigenomics datasets."""
        return self._cache.get("epigenomics", [])

    def get_dataset(self, accession: str) -> Optional[Dict]:
        """Get a specific dataset by accession across all omics types."""
        for omics_type in self.omics_types:
            for dataset in self._cache.get(omics_type, []):
                if dataset.get("accession") == accession:
                    result = dataset.copy()
                    result["_omics_type"] = omics_type
                    return result
        return None

    def search(self, query: str = "", omics_type: Optional[str] = None,
               organism: Optional[str] = None, data_source: Optional[str] = None) -> List[Dict]:
        """Search datasets with various filters."""
        results = []
        search_types = [omics_type] if omics_type else self.omics_types

        for otype in search_types:
            for dataset in self._cache.get(otype, []):
                # Text search in title, description, accession
                if query:
                    searchable = f"{dataset.get('title', '')} {dataset.get('description', '')} {dataset.get('accession', '')}".lower()
                    if query.lower() not in searchable:
                        continue

                # Filter by organism
                if organism and dataset.get("organism", "").lower() != organism.lower():
                    continue

                # Filter by data source
                if data_source and dataset.get("data_source", "").lower() != data_source.lower():
                    continue

                result = dataset.copy()
                result["_omics_type"] = otype
                results.append(result)

        return results

    def get_stats(self) -> Dict:
        """Get registry statistics."""
        stats = {
            "total_datasets": sum(len(datasets) for datasets in self._cache.values()),
            "by_omics_type": {otype: len(datasets) for otype, datasets in self._cache.items()},
            "organisms": {},
            "data_sources": {}
        }

        # Count organisms and data sources
        for omics_type, datasets in self._cache.items():
            for dataset in datasets:
                organism = dataset.get("organism", "Unknown")
                data_source = dataset.get("data_source", "Unknown")

                stats["organisms"][organism] = stats["organisms"].get(organism, 0) + 1
                stats["data_sources"][data_source] = stats["data_sources"].get(data_source, 0) + 1

        return stats

    def reload(self):
        """Reload registry from files."""
        self._cache = {}
        self._load_registry()

# Backward compatibility - maintain existing function signatures
def load_transcriptomics_registry() -> List[Dict]:
    """Backward compatibility function for transcriptomics data."""
    registry = HvantkRegistry()
    return registry.list_transcriptomics_datasets()

def load_expression_atlas_datasets() -> List[Dict]:
    """Backward compatibility function for Expression Atlas data."""
    registry = HvantkRegistry()
    return [d for d in registry.list_transcriptomics_datasets()
            if d.get("data_source") == "Expression_Atlas"]

def load_ucsc_datasets() -> List[Dict]:
    """Backward compatibility function for UCSC data."""
    registry = HvantkRegistry()
    return [d for d in registry.list_transcriptomics_datasets()
            if d.get("data_source") == "UCSC"]
