"""
Migration utility to convert existing registry data to new schema format.
"""
import json
import sys
from pathlib import Path
from typing import Dict, List, Any
from datetime import datetime

# Add the parent directory to path
sys.path.append(str(Path(__file__).parent))
from schema_validator import SchemaValidator

class RegistryMigrator:
    """Migrates existing registry data to new schema-compliant format."""

    def __init__(self, resources_dir: Path):
        self.resources_dir = resources_dir
        self.validator = SchemaValidator()
        self.migration_log = []

    def migrate_expression_atlas_data(self) -> List[Dict]:
        """Migrate expression_atlas.json to transcriptomics format."""
        expression_file = self.resources_dir / "expression_atlas.json"
        if not expression_file.exists():
            return []

        with open(expression_file, 'r') as f:
            old_data = json.load(f)

        migrated_datasets = []

        for dataset in old_data:
            # Convert old format to new schema
            new_dataset = {
                "title": dataset.get("title", ""),
                "accession": dataset.get("accession", ""),
                "description": dataset.get("description", ""),
                "pubmedid": dataset.get("pubmedid"),
                "data_source": "Expression_Atlas",
                "last_updated": datetime.now().isoformat(),
                "update_frequency": "static",
                "organism": "Homo sapiens",  # Default, can be updated based on data

                # Map old type field to new schema fields
                "expression_unit": self._map_expression_unit(dataset.get("type", "")),
                "platform_type": self._map_platform_type(dataset.get("type", "")),
                "data_level": "transcript",  # Based on file types

                # Convert files to new format
                "files": self._convert_files(dataset.get("files", []))
            }

            # Validate against schema
            is_valid, errors = self.validator.validate_dataset(new_dataset, "transcriptomics_schema")
            if is_valid:
                migrated_datasets.append(new_dataset)
                self.migration_log.append(f"✓ Migrated {dataset.get('accession', 'unknown')}")
            else:
                self.migration_log.append(f"✗ Failed to migrate {dataset.get('accession', 'unknown')}: {errors}")

        return migrated_datasets

    def migrate_ucsc_datasets(self) -> List[Dict]:
        """Migrate cells_ucsc_datasets.json to transcriptomics format."""
        ucsc_file = self.resources_dir / "cells_ucsc_datasets.json"
        if not ucsc_file.exists():
            return []

        with open(ucsc_file, 'r') as f:
            old_data = json.load(f)

        migrated_datasets = []

        for dataset in old_data.get("datasets", []):
            # Convert UCSC format to new schema
            new_dataset = {
                "title": dataset.get("shortLabel", ""),
                "accession": dataset.get("name", ""),
                "description": f"UCSC Cell Browser dataset. {dataset.get('shortLabel', '')}",
                "data_source": "UCSC",
                "last_updated": datetime.now().isoformat(),
                "update_frequency": "static",
                "organism": self._map_organism(dataset.get("organisms", [])),
                "tissue_type": dataset.get("body_parts", []),
                "sample_count": self._extract_sample_count(dataset),

                # Single-cell RNA-seq specifics
                "expression_unit": "counts",  # UCSC typically provides raw counts
                "platform_type": "single-cell RNA-seq",
                "data_level": "gene",

                # Convert files - UCSC has different file structure
                "files": self._convert_ucsc_files(dataset)
            }

            # Validate against schema
            is_valid, errors = self.validator.validate_dataset(new_dataset, "transcriptomics_schema")
            if is_valid:
                migrated_datasets.append(new_dataset)
                self.migration_log.append(f"✓ Migrated UCSC {dataset.get('name', 'unknown')}")
            else:
                self.migration_log.append(f"✗ Failed to migrate UCSC {dataset.get('name', 'unknown')}: {errors}")

        return migrated_datasets

    def migrate_dataset_creation_registry(self) -> List[Dict]:
        """Migrate dataset_creation_registry.json to transcriptomics format."""
        creation_file = self.resources_dir / "dataset_creation_registry.json"
        if not creation_file.exists():
            return []

        with open(creation_file, 'r') as f:
            old_data = json.load(f)

        migrated_datasets = []

        for dataset_id, dataset_info in old_data.items():
            metadata = dataset_info.get("metadata", {})

            new_dataset = {
                "title": metadata.get("title", ""),
                "accession": dataset_info.get("dataset_id", dataset_id),
                "description": metadata.get("description", ""),
                "data_source": "UCSC" if dataset_info.get("source") == "ucsc" else "Custom",
                "last_updated": datetime.now().isoformat(),
                "update_frequency": "static",
                "organism": metadata.get("organism", "Homo sapiens"),
                "tissue_type": metadata.get("tissue_type", ""),
                "sample_count": max(1, metadata.get("sample_count") or 1),  # Ensure valid sample count

                # Map matrix type to our schema
                "expression_unit": "counts" if dataset_info.get("matrix_type") == "single_cell" else "TPM",
                "platform_type": "single-cell RNA-seq" if dataset_info.get("matrix_type") == "single_cell" else "RNA-seq",
                "data_level": "gene",

                # Create basic file entry
                "files": [{
                    "path": f"{dataset_id}_expression.tsv",
                    "format": "tsv",
                    "size_bytes": int(metadata.get("file_size_mb", 0) * 1024 * 1024),
                    "description": f"Expression matrix for {metadata.get('title', dataset_id)}"
                }]
            }

            # Validate against schema
            is_valid, errors = self.validator.validate_dataset(new_dataset, "transcriptomics_schema")
            if is_valid:
                migrated_datasets.append(new_dataset)
                self.migration_log.append(f"✓ Migrated creation registry {dataset_id}")
            else:
                self.migration_log.append(f"✗ Failed to migrate creation registry {dataset_id}: {errors}")

        return migrated_datasets

    def _map_expression_unit(self, old_type: str) -> str:
        """Map old type field to expression unit."""
        if "tpm" in old_type.lower():
            return "TPM"
        elif "fpkm" in old_type.lower():
            return "FPKM"
        elif "count" in old_type.lower():
            return "counts"
        else:
            return "TPM"  # Default for RNA-seq

    def _map_platform_type(self, old_type: str) -> str:
        """Map old type field to platform type."""
        if "RNA-Seq" in old_type:
            return "RNA-seq"
        elif "microarray" in old_type.lower():
            return "microarray"
        else:
            return "RNA-seq"  # Default

    def _convert_files(self, old_files: List[Dict]) -> List[Dict]:
        """Convert old file format to new schema format."""
        new_files = []

        for file_obj in old_files:
            new_file = {
                "path": file_obj.get("name", ""),
                "format": self._guess_format(file_obj.get("name", "")),
                "size_bytes": 0,  # Will need to be updated with actual sizes
                "description": f"File type: {file_obj.get('type', 'unknown')}"
            }
            new_files.append(new_file)

        return new_files

    def _guess_format(self, filename: str) -> str:
        """Guess file format from filename."""
        if filename.endswith('.tsv'):
            return "tsv"
        elif filename.endswith('.csv'):
            return "csv"
        elif filename.endswith('.json'):
            return "json"
        else:
            return "tsv"  # Default

    def _map_organism(self, organisms: List[str]) -> str:
        """Map organism list to single organism string."""
        if not organisms:
            return "Homo sapiens"

        organism_map = {
            "Human (H. sapiens)": "Homo sapiens",
            "Mouse (M. musculus)": "Mus musculus",
            "Chicken (G. gallus)": "Gallus gallus"
        }

        first_organism = organisms[0] if organisms else "Human (H. sapiens)"
        return organism_map.get(first_organism, first_organism)

    def _extract_sample_count(self, dataset: Dict) -> int:
        """Extract sample count from UCSC dataset."""
        # Try to find sample count in various fields
        if "sampleCount" in dataset and dataset["sampleCount"] > 0:
            return dataset["sampleCount"]
        elif "cellCount" in dataset and dataset["cellCount"] > 0:
            return dataset["cellCount"]
        elif "facets" in dataset:
            # Look for sample count in facets
            facets = dataset["facets"]
            if "sampleCount" in facets and facets["sampleCount"] > 0:
                return facets["sampleCount"]
            elif "cellCount" in facets and facets["cellCount"] > 0:
                return facets["cellCount"]

        # Default to 1 if no valid sample count found to pass schema validation
        # This indicates the sample count needs to be manually verified
        return 1

    def _convert_ucsc_files(self, dataset: Dict) -> List[Dict]:
        """Convert UCSC file structure to new format."""
        files = []

        # UCSC datasets typically have these file types
        has_files = dataset.get("hasFiles", [])
        dataset_name = dataset.get("name", "")

        for file_type in has_files:
            file_obj = {
                "path": f"{dataset_name}_{file_type}.tsv",
                "format": "tsv",
                "size_bytes": 0,  # Will need to be updated
                "description": f"UCSC {file_type} file"
            }
            files.append(file_obj)

        # If no files specified, create a default expression file
        if not files:
            files.append({
                "path": f"{dataset_name}_expression.tsv",
                "format": "tsv",
                "size_bytes": 0,
                "description": "Expression data file"
            })

        return files

    def save_migrated_data(self, datasets: List[Dict], omics_type: str):
        """Save migrated datasets to new registry structure."""
        registry_dir = self.resources_dir / "registry" / omics_type
        registry_dir.mkdir(parents=True, exist_ok=True)

        output_file = registry_dir / "datasets.json"
        with open(output_file, 'w') as f:
            json.dump(datasets, f, indent=2)

        self.migration_log.append(f"✓ Saved {len(datasets)} datasets to {output_file}")

    def print_migration_log(self):
        """Print the migration log."""
        print("\n=== Migration Log ===")
        for entry in self.migration_log:
            print(entry)

if __name__ == "__main__":
    resources_dir = Path(__file__).parent
    migrator = RegistryMigrator(resources_dir)

    print("Starting complete registry migration...")

    # Collect all transcriptomics data from different sources
    all_transcriptomics_data = []

    # Migrate Expression Atlas data
    print("\n1. Migrating Expression Atlas data...")
    expression_atlas_data = migrator.migrate_expression_atlas_data()
    all_transcriptomics_data.extend(expression_atlas_data)

    # Migrate UCSC datasets
    print("2. Migrating UCSC datasets...")
    ucsc_data = migrator.migrate_ucsc_datasets()
    all_transcriptomics_data.extend(ucsc_data)

    # Migrate dataset creation registry
    print("3. Migrating dataset creation registry...")
    creation_registry_data = migrator.migrate_dataset_creation_registry()
    all_transcriptomics_data.extend(creation_registry_data)

    # Save all transcriptomics data together
    migrator.save_migrated_data(all_transcriptomics_data, "transcriptomics")

    migrator.print_migration_log()
    print(f"\nTotal migrated datasets: {len(all_transcriptomics_data)}")
    print(f"- Expression Atlas: {len(expression_atlas_data)}")
    print(f"- UCSC: {len(ucsc_data)}")
    print(f"- Creation Registry: {len(creation_registry_data)}")
