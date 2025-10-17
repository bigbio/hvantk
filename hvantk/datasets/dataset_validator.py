"""
Dataset Validator Orchestrator

This module provides high-level orchestration for validating datasets from UCSC and Expression Atlas.
It integrates with existing dataset classes and the validation registry to provide automated
dataset validation workflows.
"""
from __future__ import annotations

import logging
import os
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from hvantk.datasets.validation_registry import (
    DatasetValidationRegistry,
    ValidationResult,
    ValidationStatus,
    FailureType
)

logger = logging.getLogger(__name__)


class DatasetValidator:
    """
    High-level orchestrator for dataset validation workflows.

    Integrates with existing UCSC and Expression Atlas dataset classes
    to provide automated validation of datasets for matrix creation.
    """

    def __init__(self,
                 work_dir: Optional[str] = None,
                 sample_lines: int = 100,
                 registry_file: Optional[str] = None):
        """
        Initialize the dataset validator.

        Args:
            work_dir: Working directory for downloads and sample files
            sample_lines: Number of lines to use for sample validation
            registry_file: Path to validation registry file
        """
        self.work_dir = Path(work_dir) if work_dir else Path.cwd() / "dataset_validation"
        self.sample_lines = sample_lines
        self.registry = DatasetValidationRegistry(registry_file)

        # Create working directories
        self.downloads_dir = self.work_dir / "downloads"
        self.samples_dir = self.work_dir / "samples"

        self.downloads_dir.mkdir(parents=True, exist_ok=True)
        self.samples_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"Dataset validator initialized with work_dir: {self.work_dir}")

    def validate_ucsc_dataset(self, dataset_name: str, force_revalidate: bool = False) -> ValidationResult:
        """
        Validate a UCSC dataset by name.

        Args:
            dataset_name: Name of the UCSC dataset (e.g., 'cortex-dev')
            force_revalidate: Force revalidation even if already validated

        Returns:
            ValidationResult with validation status and details
        """
        # Check if already validated
        existing_result = self.registry.get_result(dataset_name)
        if existing_result and not force_revalidate:
            logger.info(f"Dataset {dataset_name} already validated: {existing_result.status.value}")
            return existing_result

        logger.info(f"Starting validation of UCSC dataset: {dataset_name}")

        try:
            # Load dataset using existing UCSC dataset class
            from hvantk.datasets.ucsc_cell_datasets import load_ucsc_datasets

            # Load all UCSC datasets to find the one we want
            all_datasets = load_ucsc_datasets()
            dataset = None
            for ds in all_datasets:
                if ds.name == dataset_name:
                    dataset = ds
                    break

            if not dataset:
                result = ValidationResult(
                    dataset_id=dataset_name,
                    dataset_type="ucsc",
                    status=ValidationStatus.DOWNLOAD_FAILED,
                    failure_type=FailureType.DOWNLOAD_ERROR,
                    error_message=f"Dataset {dataset_name} not found in UCSC catalog"
                )
                self.registry.update_result(result)
                return result

            # Download files
            dataset_dir = self.downloads_dir / "ucsc" / dataset_name
            dataset_dir.mkdir(parents=True, exist_ok=True)

            files = {}
            try:
                # Download expression matrix
                expr_path = dataset.download_expression_matrix(str(dataset_dir))
                files["expression_matrix"] = expr_path

                # Download metadata
                meta_path = dataset.download_metadata(str(dataset_dir))
                files["metadata"] = meta_path

            except ValueError as e:
                result = ValidationResult(
                    dataset_id=dataset_name,
                    dataset_type="ucsc",
                    status=ValidationStatus.DOWNLOAD_FAILED,
                    failure_type=FailureType.DOWNLOAD_ERROR,
                    error_message=f"Failed to download files: {e!s}"
                )
                self.registry.update_result(result)
                return result

            # Perform validation
            return self._validate_dataset_files(dataset_name, "ucsc", files)

        except Exception as e:
            logger.error(f"Error validating UCSC dataset {dataset_name}: {e}")
            result = ValidationResult(
                dataset_id=dataset_name,
                dataset_type="ucsc",
                status=ValidationStatus.DOWNLOAD_FAILED,
                failure_type=FailureType.UNKNOWN_ERROR,
                error_message=str(e)
            )
            self.registry.update_result(result)
            return result

    def validate_expression_atlas_dataset(self, accession: str, force_revalidate: bool = False) -> ValidationResult:
        """
        Validate an Expression Atlas dataset by accession.

        Args:
            accession: Expression Atlas accession (e.g., 'E-MTAB-5061')
            force_revalidate: Force revalidation even if already validated

        Returns:
            ValidationResult with validation status and details
        """
        # Check if already validated
        existing_result = self.registry.get_result(accession)
        if existing_result and not force_revalidate:
            logger.info(f"Dataset {accession} already validated: {existing_result.status.value}")
            return existing_result

        logger.info(f"Starting validation of Expression Atlas dataset: {accession}")

        try:
            # Load dataset using existing Expression Atlas dataset class
            from hvantk.datasets.expression_atlas_datasets import load_expression_atlas_datasets

            # Load all Expression Atlas datasets to find the one we want
            all_datasets = load_expression_atlas_datasets()
            dataset = None
            for ds in all_datasets:
                if ds.accession == accession:
                    dataset = ds
                    break

            if not dataset:
                result = ValidationResult(
                    dataset_id=accession,
                    dataset_type="expression_atlas",
                    status=ValidationStatus.DOWNLOAD_FAILED,
                    failure_type=FailureType.DOWNLOAD_ERROR,
                    error_message=f"Dataset {accession} not found in Expression Atlas catalog"
                )
                self.registry.update_result(result)
                return result

            # Download files
            dataset_dir = self.downloads_dir / "expression_atlas" / accession
            dataset_dir.mkdir(parents=True, exist_ok=True)

            files = {}
            try:
                # Download expression data and metadata using the correct methods
                expr_path = dataset.download_expression_data(str(dataset_dir))
                files["expression_matrix"] = expr_path

                meta_path = dataset.download_metadata(str(dataset_dir))
                files["metadata"] = meta_path

                if not files:
                    raise Exception("No valid expression matrix or metadata files found")

            except ValueError as e:
                result = ValidationResult(
                    dataset_id=accession,
                    dataset_type="expression_atlas",
                    status=ValidationStatus.DOWNLOAD_FAILED,
                    failure_type=FailureType.DOWNLOAD_ERROR,
                    error_message=f"Failed to download files: {e!s}"
                )
                self.registry.update_result(result)
                return result

            # Perform validation
            return self._validate_dataset_files(accession, "expression_atlas", files)

        except Exception as e:
            logger.error(f"Error validating Expression Atlas dataset {accession}: {e}")
            result = ValidationResult(
                dataset_id=accession,
                dataset_type="expression_atlas",
                status=ValidationStatus.DOWNLOAD_FAILED,
                failure_type=FailureType.UNKNOWN_ERROR,
                error_message=str(e)
            )
            self.registry.update_result(result)
            return result

    def _validate_dataset_files(self, dataset_id: str, dataset_type: str, files: Dict[str, str]) -> ValidationResult:
        """
        Perform the multi-tier validation on downloaded files.

        Args:
            dataset_id: Unique identifier for the dataset
            dataset_type: Type of dataset ("ucsc" or "expression_atlas")
            files: Dictionary mapping file types to file paths

        Returns:
            ValidationResult with complete validation results
        """
        # Tier 1: Header validation
        logger.info(f"Starting Tier 1 validation for {dataset_id}")
        result = self.registry.validate_dataset_tier1(dataset_id, dataset_type, files)

        if result.status == ValidationStatus.TIER1_PASSED:
            # Tier 2: Sample validation
            logger.info(f"Starting Tier 2 validation for {dataset_id}")
            sample_dir = str(self.samples_dir / dataset_type / dataset_id)
            result = self.registry.validate_dataset_tier2(
                result, files, sample_dir, self.sample_lines
            )

        # Update registry
        self.registry.update_result(result)
        logger.info(f"Validation completed for {dataset_id}: {result.status.value}")

        return result

    def validate_ucsc_datasets_batch(self, dataset_names: Optional[List[str]] = None,
                                   max_datasets: Optional[int] = None) -> Dict[str, ValidationResult]:
        """
        Validate multiple UCSC datasets in batch.

        Args:
            dataset_names: List of specific dataset names to validate.
                          If None, validates all available datasets.
            max_datasets: Maximum number of datasets to validate (for testing)

        Returns:
            Dictionary mapping dataset names to validation results
        """
        from hvantk.datasets.ucsc_cell_datasets import load_ucsc_datasets

        all_datasets = load_ucsc_datasets()

        if dataset_names:
            datasets_to_validate = [ds for ds in all_datasets if ds.name in dataset_names]
        else:
            datasets_to_validate = all_datasets

        if max_datasets:
            datasets_to_validate = datasets_to_validate[:max_datasets]

        logger.info(f"Starting batch validation of {len(datasets_to_validate)} UCSC datasets")

        results = {}
        for i, dataset in enumerate(datasets_to_validate, 1):
            logger.info(f"Validating dataset {i}/{len(datasets_to_validate)}: {dataset.name}")
            try:
                result = self.validate_ucsc_dataset(dataset.name)
                results[dataset.name] = result
            except Exception as e:
                logger.error(f"Failed to validate dataset {dataset.name}: {e}")
                results[dataset.name] = ValidationResult(
                    dataset_id=dataset.name,
                    dataset_type="ucsc",
                    status=ValidationStatus.DOWNLOAD_FAILED,
                    failure_type=FailureType.UNKNOWN_ERROR,
                    error_message=str(e)
                )

        logger.info(f"Batch validation completed. Results: {len(results)} datasets processed")
        return results

    def validate_expression_atlas_datasets_batch(self, accessions: Optional[List[str]] = None,
                                               max_datasets: Optional[int] = None) -> Dict[str, ValidationResult]:
        """
        Validate multiple Expression Atlas datasets in batch.

        Args:
            accessions: List of specific dataset accessions to validate.
                       If None, validates all available datasets.
            max_datasets: Maximum number of datasets to validate (for testing)

        Returns:
            Dictionary mapping dataset accessions to validation results
        """
        from hvantk.datasets.expression_atlas_datasets import load_expression_atlas_datasets

        all_datasets = load_expression_atlas_datasets()

        if accessions:
            datasets_to_validate = [ds for ds in all_datasets if ds.accession in accessions]
        else:
            datasets_to_validate = all_datasets

        if max_datasets:
            datasets_to_validate = datasets_to_validate[:max_datasets]

        logger.info(f"Starting batch validation of {len(datasets_to_validate)} Expression Atlas datasets")

        results = {}
        for i, dataset in enumerate(datasets_to_validate, 1):
            logger.info(f"Validating dataset {i}/{len(datasets_to_validate)}: {dataset.accession}")
            try:
                result = self.validate_expression_atlas_dataset(dataset.accession)
                results[dataset.accession] = result
            except Exception as e:
                logger.error(f"Failed to validate dataset {dataset.accession}: {e}")
                results[dataset.accession] = ValidationResult(
                    dataset_id=dataset.accession,
                    dataset_type="expression_atlas",
                    status=ValidationStatus.DOWNLOAD_FAILED,
                    failure_type=FailureType.UNKNOWN_ERROR,
                    error_message=str(e)
                )

        logger.info(f"Batch validation completed. Results: {len(results)} datasets processed")
        return results

    def generate_validation_report(self, output_file: Optional[str] = None) -> str:
        """
        Generate a comprehensive validation report.

        Args:
            output_file: Path to save the report. If None, returns as string.

        Returns:
            Report content as string
        """
        stats = self.registry.get_summary_stats()
        successful = self.registry.list_successful_datasets()
        failed = self.registry.list_failed_datasets()

        report_lines = [
            "Dataset Validation Report",
            "=" * 50,
            "",
            "Summary Statistics:",
            f"  Total datasets validated: {stats['total']}",
            f"  Successful datasets: {len(successful)}",
            f"  Failed datasets: {len(failed)}",
            "",
            "Status Breakdown:",
        ]

        for status, count in stats.get('by_status', {}).items():
            report_lines.append(f"  {status}: {count}")

        if successful:
            report_lines.extend([
                "",
                f"Successful Datasets ({len(successful)}):",
            ])
            for dataset_id in successful[:10]:  # Show first 10
                result = self.registry.get_result(dataset_id)
                report_lines.append(f"  {dataset_id} ({result.dataset_type})")
            if len(successful) > 10:
                report_lines.append(f"  ... and {len(successful) - 10} more")

        if failed:
            report_lines.extend([
                "",
                f"Failed Datasets ({len(failed)}):",
            ])
            for dataset_id in failed[:10]:  # Show first 10
                result = self.registry.get_result(dataset_id)
                failure_info = f" - {result.failure_type.value}" if result.failure_type else ""
                report_lines.append(f"  {dataset_id} ({result.dataset_type}){failure_info}")
            if len(failed) > 10:
                report_lines.append(f"  ... and {len(failed) - 10} more")

        report_content = "\n".join(report_lines)

        if output_file:
            with open(output_file, 'w') as f:
                f.write(report_content)
            logger.info(f"Validation report saved to {output_file}")

        return report_content
