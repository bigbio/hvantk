"""
Unified Registry Manager for HVANTK Dataset Validation

This module provides a high-level interface for managing the complete
dataset validation registry workflow, including validation, web generation,
and API endpoint creation.
"""

import logging
from pathlib import Path
from typing import Optional, Dict, Any

from .validation_registry import DatasetValidationRegistry, ValidationResult
from .api_generator import APIEndpointGenerator
from .config import RegistryConfig

try:
    from .web_generator import WebRegistryGenerator
except Exception:
    WebRegistryGenerator = None

logger = logging.getLogger(__name__)


class RegistryManager:
    """
    Unified manager for the dataset validation registry system.

    This class orchestrates all aspects of the registry including:
    - Validation result storage and retrieval
    - Web interface generation
    - API endpoint creation
    - Configuration management
    """

    def __init__(
        self,
        registry_file: Optional[str] = None,
        config: Optional[RegistryConfig] = None,
    ):
        """
        Initialize the registry manager.

        Args:
            registry_file: Path to the validation registry JSON file
            config: Registry configuration object
        """
        self.config = config or RegistryConfig()
        self.registry = DatasetValidationRegistry(registry_file)
        self.web_generator = (
            WebRegistryGenerator(self.config)
            if WebRegistryGenerator is not None
            else None
        )
        self.api_generator = APIEndpointGenerator(self.config)

        logger.info("Registry manager initialized")

    def get_validation_result(self, dataset_id: str) -> Optional[ValidationResult]:
        """Get validation result for a specific dataset."""
        return self.registry.get_result(dataset_id)

    def update_validation_result(self, result: ValidationResult) -> None:
        """Update or add a validation result."""
        self.registry.update_result(result)
        logger.info(
            f"Updated validation result for {result.dataset_id}: {result.status.value}"
        )

    def get_all_results(self) -> Dict[str, ValidationResult]:
        """Get all validation results."""
        return self.registry.results

    def save_registry(self) -> None:
        """Save validation results to file."""
        self.registry.save_registry()

    def generate_web_interface(self, output_dir: str) -> Path:
        """
        Generate the complete web interface for the registry.

        Args:
            output_dir: Directory to output web files

        Returns:
            Path to generated web interface

        Raises:
            RuntimeError: If web generator is not configured
        """
        if self.web_generator is None:
            raise RuntimeError(
                "WebRegistryGenerator is not configured; cannot generate web interface"
            )

        logger.info(f"Generating web interface in {output_dir}")
        return self.web_generator.generate_web_registry(
            str(self.registry.registry_file), output_dir
        )

    def generate_api_endpoints(self, output_dir: str) -> Path:
        """
        Generate API endpoints for the registry.

        Args:
            output_dir: Directory to output API files

        Returns:
            Path to generated API directory
        """
        logger.info(f"Generating API endpoints in {output_dir}")
        return self.api_generator.generate_api_endpoints(
            str(self.registry.registry_file), output_dir
        )

    def generate_complete_registry(self, output_dir: str) -> Dict[str, Path]:
        """
        Generate both web interface and API endpoints.

        Args:
            output_dir: Directory to output all files

        Returns:
            Dictionary with paths to generated components
        """
        logger.info(f"Generating complete registry in {output_dir}")

        result = {}

        # Generate web interface (only if available)
        if self.web_generator is not None:
            result["web_interface"] = self.generate_web_interface(output_dir)
        else:
            logger.warning(
                "WebRegistryGenerator not available - skipping web interface generation"
            )

        # Generate API endpoints
        result["api_endpoints"] = self.generate_api_endpoints(output_dir)

        return result

    def get_registry_statistics(self) -> Dict[str, Any]:
        """Get comprehensive statistics about the registry."""
        results = self.registry.results

        if not results:
            return {
                "total_datasets": 0,
                "successful": 0,
                "failed": 0,
                "success_rate": 0.0,
                "by_status": {},
                "by_type": {},
            }

        from collections import Counter

        total = len(results)
        status_counts = Counter(result.status.value for result in results.values())
        type_counts = Counter(result.dataset_type for result in results.values())

        successful = sum(
            count
            for status, count in status_counts.items()
            if status in ["tier3_passed", "tier2_passed"]
        )

        return {
            "total_datasets": total,
            "successful": successful,
            "failed": total - successful,
            "success_rate": round((successful / total * 100) if total > 0 else 0, 1),
            "by_status": dict(status_counts),
            "by_type": dict(type_counts),
        }

    def get_failed_datasets(self, limit: int = 10) -> list:
        """Get list of recently failed datasets."""
        failed_results = [
            (dataset_id, result)
            for dataset_id, result in self.registry.results.items()
            if result.status.value
            in ["tier1_failed", "tier2_failed", "tier3_failed", "download_failed"]
        ]

        # Sort by timestamp (most recent first)
        failed_results.sort(key=lambda x: x[1].timestamp, reverse=True)

        return [
            {
                "dataset_id": dataset_id,
                "status": result.status.value,
                "timestamp": result.timestamp,
                "error_message": result.error_message or "No error message",
            }
            for dataset_id, result in failed_results[:limit]
        ]

    def cleanup_old_results(self, days: int = 30) -> int:
        """
        Remove validation results older than specified days.

        Args:
            days: Number of days to keep results

        Returns:
            Number of results removed
        """
        from datetime import datetime, timedelta

        cutoff_date = datetime.now() - timedelta(days=days)
        initial_count = len(self.registry.results)

        # Filter out old results
        self.registry.results = {
            dataset_id: result
            for dataset_id, result in self.registry.results.items()
            if datetime.fromisoformat(result.timestamp) > cutoff_date
        }

        removed_count = initial_count - len(self.registry.results)

        if removed_count > 0:
            self.save_registry()
            logger.info(f"Cleaned up {removed_count} old validation results")

        return removed_count
