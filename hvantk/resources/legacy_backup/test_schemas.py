"""
Test script to validate the schema framework works correctly.
"""
import json
import sys
from pathlib import Path

# Add the parent directory to path to import schema_validator
sys.path.append(str(Path(__file__).parent))
from schema_validator import SchemaValidator

def test_schemas():
    """Test the schema validation framework."""
    validator = SchemaValidator()

    # Test data for each schema type
    test_cases = {
        "common_schema": {
            "title": "Test Dataset",
            "accession": "TEST001",
            "description": "Test dataset for validation",
            "files": [{
                "path": "test.tsv",
                "format": "tsv",
                "size_bytes": 1024
            }],
            "data_source": "Custom"
        },
        "transcriptomics_schema": {
            "title": "Test RNA-seq Dataset",
            "accession": "RNASEQ001",
            "description": "Test RNA-seq dataset",
            "files": [{
                "path": "expression.tsv",
                "format": "tsv",
                "size_bytes": 2048
            }],
            "data_source": "GEO",
            "expression_unit": "TPM",
            "platform_type": "RNA-seq"
        }
    }

    print("Testing schema validation framework...")

    for schema_type, test_data in test_cases.items():
        is_valid, errors = validator.validate_dataset(test_data, schema_type)
        if is_valid:
            print(f"✓ {schema_type}: PASSED")
        else:
            print(f"✗ {schema_type}: FAILED - {errors}")

    print("\nSchema requirements for transcriptomics:")
    requirements = validator.get_schema_requirements("transcriptomics_schema")
    print(f"Required fields: {requirements['required_fields']}")

if __name__ == "__main__":
    test_schemas()
