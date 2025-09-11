#!/usr/bin/env python3
"""
Automated UCSC Dataset Validation Script

This script systematically validates all available UCSC datasets for Hail MatrixTable
creation compatibility. It runs validation in batches with proper error handling,
logging, and progress tracking.
"""

import argparse
import json
import logging
import os
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Any

# Add the project root to Python path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from hvantk.datasets.dataset_validator import DatasetValidator
from hvantk.datasets.ucsc_cell_datasets import load_ucsc_datasets
from hvantk.datasets.validation_registry import ValidationStatus


class UCscDatasetValidationRunner:
    """Automated validation runner for UCSC datasets."""

    def __init__(self, work_dir: str = "./dataset_validation",
                 sample_lines: int = 100,
                 batch_size: int = 5,
                 max_datasets: int = None,
                 continue_on_error: bool = True):
        """
        Initialize the validation runner.

        Args:
            work_dir: Working directory for validation files
            sample_lines: Number of lines to use for sample validation
            batch_size: Number of datasets to process in each batch
            max_datasets: Maximum number of datasets to validate (None = all)
            continue_on_error: Whether to continue validation if a dataset fails
        """
        self.work_dir = work_dir
        self.sample_lines = sample_lines
        self.batch_size = batch_size
        self.max_datasets = max_datasets
        self.continue_on_error = continue_on_error

        # Setup logging
        self.setup_logging()

        # Initialize validator
        self.validator = DatasetValidator(
            work_dir=work_dir,
            sample_lines=sample_lines
        )

        # Load available datasets
        self.ucsc_datasets = load_ucsc_datasets()

        # Validation tracking
        self.validation_results = {}
        self.successful_count = 0
        self.failed_count = 0
        self.skipped_count = 0

    def setup_logging(self):
        """Setup comprehensive logging."""
        # Create logs directory
        log_dir = Path(self.work_dir) / "logs"
        log_dir.mkdir(parents=True, exist_ok=True)

        # Setup file logging
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file = log_dir / f"ucsc_validation_{timestamp}.log"

        # Configure logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )

        self.logger = logging.getLogger(__name__)
        self.logger.info(f"Starting UCSC dataset validation - Log file: {log_file}")

    def get_datasets_to_validate(self) -> List[str]:
        """Get list of dataset names to validate."""
        all_dataset_names = [dataset.name for dataset in self.ucsc_datasets]

        # Check which datasets have already been validated
        already_validated = set()
        try:
            # Get existing validation results using the correct method
            existing_results = self.validator.registry.list_successful_datasets()
            failed_results = self.validator.registry.list_failed_datasets()
            already_validated.update(existing_results)
            already_validated.update(failed_results)

            # Also check for any UCSC datasets in the registry
            registry_stats = self.validator.registry.get_summary_stats()
            if registry_stats.get('total', 0) > 0:
                # Load the registry data to get dataset IDs
                registry_data = getattr(self.validator.registry, '_registry_data', {})
                for dataset_id in registry_data.keys():
                    if dataset_id in all_dataset_names:
                        already_validated.add(dataset_id)

        except Exception as e:
            self.logger.warning(f"Could not load existing validation results: {e}")

        # Filter out already validated datasets
        datasets_to_validate = [name for name in all_dataset_names if name not in already_validated]

        self.logger.info(f"Total UCSC datasets: {len(all_dataset_names)}")
        self.logger.info(f"Already validated: {len(already_validated)}")
        self.logger.info(f"Remaining to validate: {len(datasets_to_validate)}")

        # Apply max_datasets limit if specified
        if self.max_datasets:
            datasets_to_validate = datasets_to_validate[:self.max_datasets]
            self.logger.info(f"Limited to first {self.max_datasets} datasets")

        return datasets_to_validate

    def validate_batch(self, dataset_names: List[str]) -> Dict[str, Any]:
        """Validate a batch of datasets."""
        self.logger.info(f"Validating batch: {dataset_names}")

        try:
            # Use the existing batch validation method
            results = self.validator.validate_ucsc_datasets_batch(
                dataset_names=dataset_names,
                max_datasets=len(dataset_names)
            )

            # Update counters
            for dataset_id, result in results.items():
                if result.status in [ValidationStatus.TIER2_PASSED, ValidationStatus.TIER3_PASSED]:
                    self.successful_count += 1
                    self.logger.info(f"✅ {dataset_id}: PASSED")
                else:
                    self.failed_count += 1
                    self.logger.warning(f"❌ {dataset_id}: FAILED - {result.error_message}")

            return results

        except Exception as e:
            self.logger.error(f"Batch validation failed: {e}")
            if not self.continue_on_error:
                raise

            # Mark all datasets in batch as failed
            for dataset_name in dataset_names:
                self.failed_count += 1
                self.logger.error(f"❌ {dataset_name}: BATCH_FAILED - {str(e)}")

            return {}

    def run_validation(self) -> Dict[str, Any]:
        """Run complete validation process."""
        start_time = time.time()

        # Get datasets to validate
        datasets_to_validate = self.get_datasets_to_validate()

        if not datasets_to_validate:
            self.logger.info("No datasets to validate - all already processed!")
            return self.generate_summary()

        self.logger.info(f"Starting validation of {len(datasets_to_validate)} UCSC datasets")
        self.logger.info(f"Batch size: {self.batch_size}, Sample lines: {self.sample_lines}")

        # Process datasets in batches
        total_batches = (len(datasets_to_validate) + self.batch_size - 1) // self.batch_size

        for batch_idx in range(0, len(datasets_to_validate), self.batch_size):
            batch_num = (batch_idx // self.batch_size) + 1
            batch_datasets = datasets_to_validate[batch_idx:batch_idx + self.batch_size]

            self.logger.info(f"\n{'='*60}")
            self.logger.info(f"BATCH {batch_num}/{total_batches}: Processing {len(batch_datasets)} datasets")
            self.logger.info(f"{'='*60}")

            # Validate batch
            batch_results = self.validate_batch(batch_datasets)
            self.validation_results.update(batch_results)

            # Progress update
            processed = batch_idx + len(batch_datasets)
            progress = (processed / len(datasets_to_validate)) * 100

            self.logger.info(f"\nProgress: {processed}/{len(datasets_to_validate)} ({progress:.1f}%)")
            self.logger.info(f"Success: {self.successful_count}, Failed: {self.failed_count}")

            # Small delay between batches to avoid overwhelming the system
            if batch_num < total_batches:
                self.logger.info("Waiting 5 seconds before next batch...")
                time.sleep(5)

        # Final summary
        elapsed_time = time.time() - start_time
        self.logger.info(f"\n{'='*60}")
        self.logger.info(f"VALIDATION COMPLETE!")
        self.logger.info(f"{'='*60}")
        self.logger.info(f"Total time: {elapsed_time:.1f} seconds")

        return self.generate_summary()

    def generate_summary(self) -> Dict[str, Any]:
        """Generate validation summary."""
        summary = {
            'timestamp': datetime.now().isoformat(),
            'total_processed': self.successful_count + self.failed_count,
            'successful': self.successful_count,
            'failed': self.failed_count,
            'success_rate': (self.successful_count / max(1, self.successful_count + self.failed_count)) * 100,
            'validation_results': self.validation_results
        }

        # Save summary to file
        summary_file = Path(self.work_dir) / "logs" / f"validation_summary_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2, default=str)

        self.logger.info(f"Summary saved to: {summary_file}")
        return summary

    def print_final_report(self):
        """Print a comprehensive final report."""
        print(f"\n{'='*80}")
        print(f"🔬 UCSC DATASET VALIDATION REPORT")
        print(f"{'='*80}")
        print(f"📊 Total Processed: {self.successful_count + self.failed_count}")
        print(f"✅ Successful: {self.successful_count}")
        print(f"❌ Failed: {self.failed_count}")

        if self.successful_count + self.failed_count > 0:
            success_rate = (self.successful_count / (self.successful_count + self.failed_count)) * 100
            print(f"📈 Success Rate: {success_rate:.1f}%")

        print(f"\n📁 Results saved in: {self.work_dir}")
        print(f"📋 View detailed status with:")
        print(f"   python -m hvantk.commands.dataset_validation_cli status")
        print(f"📋 View successful datasets with:")
        print(f"   python -m hvantk.commands.dataset_validation_cli validated")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Automatically validate all UCSC datasets for Hail MatrixTable creation"
    )
    parser.add_argument(
        "--work-dir",
        default="./dataset_validation",
        help="Working directory for validation files (default: ./dataset_validation)"
    )
    parser.add_argument(
        "--sample-lines",
        type=int,
        default=100,
        help="Number of lines for sample validation (default: 100)"
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=5,
        help="Number of datasets to process in each batch (default: 5)"
    )
    parser.add_argument(
        "--max-datasets",
        type=int,
        help="Maximum number of datasets to validate (for testing)"
    )
    parser.add_argument(
        "--continue-on-error",
        action="store_true",
        default=True,
        help="Continue validation even if some datasets fail (default: True)"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be validated without actually running validation"
    )

    args = parser.parse_args()

    # Create validation runner
    runner = UCscDatasetValidationRunner(
        work_dir=args.work_dir,
        sample_lines=args.sample_lines,
        batch_size=args.batch_size,
        max_datasets=args.max_datasets,
        continue_on_error=args.continue_on_error
    )

    if args.dry_run:
        datasets_to_validate = runner.get_datasets_to_validate()
        print(f"Would validate {len(datasets_to_validate)} datasets:")
        for i, dataset in enumerate(datasets_to_validate[:20], 1):  # Show first 20
            print(f"  {i}. {dataset}")
        if len(datasets_to_validate) > 20:
            print(f"  ... and {len(datasets_to_validate) - 20} more")
        return

    try:
        # Run validation
        summary = runner.run_validation()
        runner.print_final_report()

    except KeyboardInterrupt:
        print("\n\n⚠️  Validation interrupted by user")
        runner.print_final_report()
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Validation failed with error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
