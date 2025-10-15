"""
Script to correct the registry by moving E-PROT datasets from transcriptomics to proteomics.
"""
import json
import sys
from pathlib import Path
from datetime import datetime

# Add the parent directory to path
sys.path.append(str(Path(__file__).parent))
from schema_validator import SchemaValidator

def move_proteomics_datasets():
    """Move E-PROT datasets from transcriptomics to proteomics registry."""

    resources_dir = Path(__file__).parent
    transcriptomics_file = resources_dir / "registry" / "transcriptomics" / "datasets.json"
    proteomics_file = resources_dir / "registry" / "proteomics" / "datasets.json"

    validator = SchemaValidator()

    # Load current transcriptomics data
    with open(transcriptomics_file, 'r') as f:
        transcriptomics_data = json.load(f)

    # Load current proteomics data (should be empty)
    with open(proteomics_file, 'r') as f:
        proteomics_data = json.load(f)

    # Find and extract E-PROT datasets
    proteomics_datasets = []
    remaining_transcriptomics = []

    for dataset in transcriptomics_data:
        if dataset.get("accession", "").startswith("E-PROT"):
            # Convert to proteomics schema
            proteomics_dataset = convert_to_proteomics_schema(dataset)

            # Validate against proteomics schema
            is_valid, errors = validator.validate_dataset(proteomics_dataset, "proteomics_schema")
            if is_valid:
                proteomics_datasets.append(proteomics_dataset)
                print(f"✓ Moved {dataset['accession']} to proteomics registry")
            else:
                print(f"✗ Failed to validate {dataset['accession']} as proteomics: {errors}")
                remaining_transcriptomics.append(dataset)  # Keep in transcriptomics if validation fails
        else:
            remaining_transcriptomics.append(dataset)

    # Add to existing proteomics data
    proteomics_data.extend(proteomics_datasets)

    # Save updated files
    with open(transcriptomics_file, 'w') as f:
        json.dump(remaining_transcriptomics, f, indent=2)

    with open(proteomics_file, 'w') as f:
        json.dump(proteomics_data, f, indent=2)

    print(f"\nCorrection Summary:")
    print(f"- Moved {len(proteomics_datasets)} datasets to proteomics registry")
    print(f"- Remaining transcriptomics datasets: {len(remaining_transcriptomics)}")
    print(f"- Total proteomics datasets: {len(proteomics_data)}")

def convert_to_proteomics_schema(transcriptomics_dataset):
    """Convert a transcriptomics dataset to proteomics schema format."""

    # Map the dataset to proteomics schema
    proteomics_dataset = {
        "title": transcriptomics_dataset.get("title", ""),
        "accession": transcriptomics_dataset.get("accession", ""),
        "description": transcriptomics_dataset.get("description", ""),
        "pubmedid": transcriptomics_dataset.get("pubmedid"),
        "data_source": transcriptomics_dataset.get("data_source", "Expression_Atlas"),
        "last_updated": datetime.now().isoformat(),
        "update_frequency": transcriptomics_dataset.get("update_frequency", "static"),
        "organism": transcriptomics_dataset.get("organism", "Homo sapiens"),
        "tissue_type": transcriptomics_dataset.get("tissue_type", ""),
        "sample_count": transcriptomics_dataset.get("sample_count", 1),
        "license": transcriptomics_dataset.get("license"),
        "doi": transcriptomics_dataset.get("doi"),
        "contact": transcriptomics_dataset.get("contact"),

        # Proteomics-specific fields
        "quantification_method": "intensity_based",  # Default for Human Protein Atlas
        "protein_inference_method": "parsimony",  # Default assumption
        "abundance_unit": "intensity",  # Default for protein abundance
        "protein_database": "UniProt",  # Standard for Human Protein Atlas
        "database_version": "current",  # Will need to be updated with actual version

        # Convert files
        "files": transcriptomics_dataset.get("files", [])
    }

    return proteomics_dataset

if __name__ == "__main__":
    move_proteomics_datasets()
