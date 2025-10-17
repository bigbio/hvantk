#!/usr/bin/env python3
"""
Test script for the dataset validation framework.

This script demonstrates the combined approach for validating datasets
and provides examples of how to use the new validation system.
"""
import logging
import tempfile
from pathlib import Path

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def test_dataset_loading():
    """Test loading datasets from JSON catalogs."""
    print("=" * 60)
    print("Testing Dataset Loading")
    print("=" * 60)

    try:
        # Test UCSC dataset loading
        from hvantk.datasets.ucsc_cell_datasets import load_ucsc_datasets
        ucsc_datasets = load_ucsc_datasets()
        print(f"✓ Loaded {len(ucsc_datasets)} UCSC datasets")

        # Show first few datasets
        for i, dataset in enumerate(ucsc_datasets[:3]):
            print(f"  {i+1}. {dataset.name} - {dataset.shortLabel}")

        # Test Expression Atlas dataset loading
        from hvantk.datasets.expression_atlas_datasets import load_expression_atlas_datasets
        atlas_datasets = load_expression_atlas_datasets()
        print(f"✓ Loaded {len(atlas_datasets)} Expression Atlas datasets")

        # Show first few datasets
        for i, dataset in enumerate(atlas_datasets[:3]):
            print(f"  {i+1}. {dataset.accession} - {dataset.title}")

        return True

    except Exception as e:
        print(f"✗ Dataset loading failed: {e}")
        return False

def test_validation_registry():
    """Test the validation registry functionality."""
    print("\n" + "=" * 60)
    print("Testing Validation Registry")
    print("=" * 60)

    try:
        from hvantk.datasets.validation_registry import (
            DatasetValidationRegistry,
            ValidationResult,
            ValidationStatus,
            FailureType
        )

        # Create a temporary registry for testing
        with tempfile.TemporaryDirectory() as temp_dir:
            registry_file = Path(temp_dir) / "test_registry.json"
            registry = DatasetValidationRegistry(str(registry_file))

            # Create a test validation result
            test_result = ValidationResult(
                dataset_id="test-dataset",
                dataset_type="ucsc",
                status=ValidationStatus.TIER2_PASSED,
                tier1_details={"expression_matrix": {"readable": True, "delimiter_detected": "\t"}},
                tier2_details={"sample_files_created": True, "matrix_creation_success": True}
            )

            # Update registry
            registry.update_result(test_result)
            print("✓ Created and saved test validation result")

            # Test retrieval
            retrieved = registry.get_result("test-dataset")
            if retrieved and retrieved.status == ValidationStatus.TIER2_PASSED:
                print("✓ Successfully retrieved validation result")
            else:
                print("✗ Failed to retrieve validation result")
                return False

            # Test statistics
            stats = registry.get_summary_stats()
            print(f"✓ Registry statistics: {stats}")

            return True

    except Exception as e:
        print(f"✗ Validation registry test failed: {e}")
        return False

def test_file_validation():
    """Test file header validation with sample data."""
    print("\n" + "=" * 60)
    print("Testing File Validation")
    print("=" * 60)

    try:
        from hvantk.datasets.validation_registry import DatasetValidationRegistry

        with tempfile.TemporaryDirectory() as temp_dir:
            registry = DatasetValidationRegistry()

            # Create a sample file for testing
            sample_file = Path(temp_dir) / "test_matrix.txt"
            with open(sample_file, 'w') as f:
                f.write("gene\tsample1\tsample2\tsample3\n")
                f.write("GENE1\t1.5\t2.3\t0.8\n")
                f.write("GENE2\t0.2\t1.1\t3.4\n")

            # Test header validation
            details = registry.validate_file_header(str(sample_file), ["gene"])

            if details["readable"] and details["delimiter_detected"] == "\t":
                print("✓ File header validation successful")
                print(f"  Detected delimiter: '{details['delimiter_detected']}'")
                print(f"  Number of columns: {details['num_columns']}")
                return True
            else:
                print("✗ File header validation failed")
                print(f"  Details: {details}")
                return False

    except Exception as e:
        print(f"✗ File validation test failed: {e}")
        return False

def test_sample_file_creation():
    """Test sample file creation functionality."""
    print("\n" + "=" * 60)
    print("Testing Sample File Creation")
    print("=" * 60)

    try:
        from hvantk.datasets.validation_registry import DatasetValidationRegistry

        with tempfile.TemporaryDirectory() as temp_dir:
            registry = DatasetValidationRegistry()

            # Create a larger test file
            input_file = Path(temp_dir) / "large_test_file.txt"
            with open(input_file, 'w') as f:
                f.write("gene\tsample1\tsample2\n")
                for i in range(200):  # Create 200 lines
                    f.write(f"GENE{i}\t{i*1.1}\t{i*0.9}\n")

            # Create sample file with first 10 lines
            output_file = Path(temp_dir) / "sample_file.txt"
            success = registry.create_sample_file(str(input_file), str(output_file), 10)

            if success and output_file.exists():
                # Count lines in output file
                with open(output_file, 'r') as f:
                    lines = f.readlines()

                if len(lines) == 10:
                    print("✓ Sample file creation successful")
                    print(f"  Input file: {len(open(input_file).readlines())} lines")
                    print(f"  Sample file: {len(lines)} lines")
                    return True
                else:
                    print(f"✗ Wrong number of lines in sample file: {len(lines)}")
                    return False
            else:
                print("✗ Sample file creation failed")
                return False

    except Exception as e:
        print(f"✗ Sample file creation test failed: {e}")
        return False

def demo_basic_usage():
    """Demonstrate basic usage of the validation system."""
    print("\n" + "=" * 60)
    print("Basic Usage Demo")
    print("=" * 60)

    try:
        from hvantk.datasets.dataset_validator import DatasetValidator

        # Create validator with temporary work directory
        with tempfile.TemporaryDirectory() as temp_dir:
            validator = DatasetValidator(
                work_dir=temp_dir,
                sample_lines=50,  # Use small sample for demo
            )

            print(f"✓ Created dataset validator with work directory: {temp_dir}")
            print(f"  Downloads directory: {validator.downloads_dir}")
            print(f"  Samples directory: {validator.samples_dir}")

            # Show registry status
            stats = validator.registry.get_summary_stats()
            print(f"✓ Registry initialized with {stats.get('total', 0)} existing results")

            return True

    except Exception as e:
        print(f"✗ Basic usage demo failed: {e}")
        return False

def main():
    """Run all tests."""
    print("Dataset Validation Framework Test Suite")
    print("=" * 60)

    tests = [
        ("Dataset Loading", test_dataset_loading),
        ("Validation Registry", test_validation_registry),
        ("File Validation", test_file_validation),
        ("Sample File Creation", test_sample_file_creation),
        ("Basic Usage Demo", demo_basic_usage),
    ]

    results = {}
    for test_name, test_func in tests:
        try:
            results[test_name] = test_func()
        except Exception as e:
            print(f"✗ {test_name} crashed: {e}")
            results[test_name] = False

    # Summary
    print("\n" + "=" * 60)
    print("Test Results Summary")
    print("=" * 60)

    passed = sum(results.values())
    total = len(results)

    for test_name, result in results.items():
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{status} {test_name}")

    print(f"\nTests passed: {passed}/{total}")

    if passed == total:
        print("\n🎉 All tests passed! The dataset validation framework is ready to use.")
        print("\nNext steps:")
        print("1. Try: python -m hvantk.commands.dataset_validation_cli list")
        print("2. Try: python -m hvantk.commands.dataset_validation_cli validate --ucsc-datasets cortex-dev")
        print("3. Try: python -m hvantk.commands.dataset_validation_cli status")
    else:
        print(f"\n⚠️  {total - passed} tests failed. Please check the errors above.")

    return passed == total

if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)
