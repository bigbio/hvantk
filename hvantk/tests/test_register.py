"""
Tests for the HVANTK Registry Module

This module contains comprehensive tests for the register functionality.
"""
import json
import tempfile
import shutil
from pathlib import Path
import pytest

from hvantk.register import (
    ValidationStatus, FailureType, ValidationResult,
    DatasetValidationRegistry, WebRegistryGenerator,
    APIEndpointGenerator, RegistryConfig, RegistryManager,
    VALIDATION_TIERS, STATUS_DEFINITIONS
)


class TestRegistryConfig:
    """Test the RegistryConfig class."""

    def test_default_values(self):
        """Test default configuration values."""
        config = RegistryConfig()

        assert config.site_title == "HVANTK Dataset Validation Registry"
        assert config.theme_color == "#667eea"
        assert config.max_datasets_per_run == 100
        assert "tier3_passed" in config.chart_colors


class TestValidationRegistry:
    """Test the DatasetValidationRegistry class."""

    def test_registry_operations(self):
        """Test basic registry operations."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            json.dump({}, f)
            temp_file = f.name

        try:
            registry = DatasetValidationRegistry(temp_file)

            # Test adding a result
            result = ValidationResult(
                dataset_id="test-dataset",
                dataset_type="ucsc",
                status=ValidationStatus.TIER2_PASSED,
                error_message="Test error"
            )

            registry.update_result(result)

            # Test retrieving result
            retrieved = registry.get_result("test-dataset")
            assert retrieved is not None
            assert retrieved.dataset_id == "test-dataset"
            assert retrieved.status == ValidationStatus.TIER2_PASSED

            # Test saving and loading
            registry.save_registry()

            # Create new registry instance and verify data persists
            registry2 = DatasetValidationRegistry(temp_file)
            retrieved2 = registry2.get_result("test-dataset")
            assert retrieved2 is not None
            assert retrieved2.dataset_id == "test-dataset"

        finally:
            Path(temp_file).unlink(missing_ok=True)


class TestWebGenerator:
    """Test the WebRegistryGenerator class."""

    def test_badge_generation(self):
        """Test status badge generation."""
        config = RegistryConfig()
        generator = WebRegistryGenerator(config)

        badge = generator.generate_status_badge("tier3_passed")
        assert "badge" in badge
        assert "#28a745" in badge  # Check for the green color
        assert "Tier3 Passed" in badge  # Check for formatted text

    def test_statistics_generation(self):
        """Test statistics generation."""
        config = RegistryConfig()
        generator = WebRegistryGenerator(config)

        # Create test registry data
        test_registry = {
            "test-dataset-1": {
                "dataset_type": "ucsc",
                "status": "tier3_passed",
                "timestamp": "2025-09-11T10:00:00",
                "error_message": None
            },
            "test-dataset-2": {
                "dataset_type": "expression_atlas",
                "status": "tier1_failed",
                "timestamp": "2025-09-11T11:00:00",
                "error_message": "Test error message"
            }
        }

        stats = generator.generate_summary_stats(test_registry)
        assert stats['total'] == 2
        assert stats['successful'] == 1
        assert stats['failed'] == 1
        assert stats['success_rate'] == 50.0


class TestAPIGenerator:
    """Test the APIEndpointGenerator class."""

    def test_api_generation(self):
        """Test API endpoint generation."""
        config = RegistryConfig()
        generator = APIEndpointGenerator(config)

        # Create test registry data
        test_registry = {
            "test-dataset-1": {
                "dataset_type": "ucsc",
                "status": "tier3_passed",
                "timestamp": "2025-09-11T10:00:00",
                "error_message": None
            },
            "test-dataset-2": {
                "dataset_type": "expression_atlas",
                "status": "tier1_failed",
                "timestamp": "2025-09-11T11:00:00",
                "error_message": "Test error message"
            }
        }

        # Test status API generation
        status_api = generator.generate_status_api(test_registry)
        assert status_api['status'] == 'ok'
        assert status_api['summary']['total_datasets'] == 2
        assert status_api['summary']['successful_validations'] == 1
        assert status_api['summary']['success_rate'] == 50.0

        # Test datasets API generation
        datasets_api = generator.generate_datasets_api(test_registry)
        assert 'datasets' in datasets_api
        assert 'metadata' in datasets_api
        assert datasets_api['metadata']['total_count'] == 2

        # Test stats API generation
        stats_api = generator.generate_stats_api(test_registry)
        assert 'statistics' in stats_api
        assert 'recent_failures' in stats_api
        assert len(stats_api['recent_failures']) == 1  # One failed dataset


class TestRegistryManager:
    """Test the RegistryManager integration."""

    def test_registry_manager_integration(self):
        """Test complete registry manager functionality."""
        # Create temporary directories
        temp_dir = Path(tempfile.mkdtemp())
        registry_file = temp_dir / "test_registry.json"
        output_dir = temp_dir / "output"

        try:
            # Initialize manager
            config = RegistryConfig()
            manager = RegistryManager(str(registry_file), config)

            # Add some test data
            result1 = ValidationResult(
                dataset_id="test-ucsc-1",
                dataset_type="ucsc",
                status=ValidationStatus.TIER3_PASSED
            )

            result2 = ValidationResult(
                dataset_id="test-atlas-1",
                dataset_type="expression_atlas",
                status=ValidationStatus.TIER1_FAILED,
                error_message="Test validation error"
            )

            manager.update_validation_result(result1)
            manager.update_validation_result(result2)
            manager.save_registry()

            # Test statistics
            stats = manager.get_registry_statistics()
            assert stats['total_datasets'] == 2
            assert stats['successful'] == 1
            assert stats['failed'] == 1
            assert stats['success_rate'] == 50.0

            # Test failed datasets
            failed = manager.get_failed_datasets()
            assert len(failed) == 1
            assert failed[0]['dataset_id'] == 'test-atlas-1'

            # Test web interface generation
            web_path = manager.generate_web_interface(str(output_dir))
            assert (output_dir / "index.html").exists()
            assert (output_dir / "styles.css").exists()

            # Test API generation
            api_path = manager.generate_api_endpoints(str(output_dir))
            assert (output_dir / "api" / "status.json").exists()
            assert (output_dir / "api" / "datasets.json").exists()
            assert (output_dir / "api" / "stats.json").exists()

            # Test complete registry generation
            paths = manager.generate_complete_registry(str(output_dir))
            assert 'web_interface' in paths
            assert 'api_endpoints' in paths

        finally:
            shutil.rmtree(temp_dir, ignore_errors=True)
