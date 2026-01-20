"""
CLI command for dataset validation and creation registry management.

Provides command-line interface for validating UCSC and Expression Atlas datasets
for Hail MatrixTable creation compatibility, and managing the dataset creation registry.
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import List, Optional

from hvantk.datasets.dataset_validator import DatasetValidator
from hvantk.datasets.validation_registry import ValidationStatus
from hvantk.datasets.dataset_creation_registry import (
    DatasetCreationRegistry,
    DatasetSource,
    MatrixType,
    CreationStatus,
)

logger = logging.getLogger(__name__)


def setup_logging(verbose: bool = False):
    """Setup logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )


def validate_datasets_command(args):
    """Handle the validate-datasets command."""
    setup_logging(args.verbose)

    # Initialize validator
    validator = DatasetValidator(
        work_dir=args.work_dir,
        sample_lines=args.sample_lines,
        registry_file=args.registry_file,
    )

    results = {}

    # Handle quick validation mode (metadata-only)
    if getattr(args, "quick_validation", False):
        logger.info(
            "Running in quick validation mode - metadata and availability checks only"
        )

        # Quick validation for UCSC datasets
        if args.ucsc_datasets:
            from hvantk.datasets.ucsc_cell_datasets import load_ucsc_datasets

            try:
                all_datasets = load_ucsc_datasets()
                logger.info(
                    f"✓ Successfully loaded UCSC catalog: {len(all_datasets)} datasets available"
                )

                # Test a few specific datasets for availability
                test_datasets = (
                    args.ucsc_datasets
                    if args.ucsc_datasets != ["all"]
                    else ["cortex-dev", "zeisel2015"]
                )
                for dataset_name in test_datasets[:3]:  # Limit to 3 for quick check
                    dataset = next(
                        (ds for ds in all_datasets if ds.name == dataset_name), None
                    )
                    if dataset:
                        logger.info(f"✓ Dataset {dataset_name} found in catalog")
                        # Create a mock successful result for quick validation
                        from hvantk.datasets.validation_registry import (
                            ValidationResult,
                            ValidationStatus,
                        )

                        results[dataset_name] = ValidationResult(
                            dataset_id=dataset_name,
                            dataset_type="ucsc",
                            status=ValidationStatus.TIER1_PASSED,
                            error_message=None,
                        )
                    else:
                        logger.warning(f"✗ Dataset {dataset_name} not found in catalog")

            except Exception as e:
                logger.error(f"Failed to load UCSC catalog: {e}")

        # Quick validation for Expression Atlas datasets
        if args.expression_atlas_datasets:
            from hvantk.datasets.expression_atlas_datasets import (
                load_expression_atlas_datasets,
            )

            try:
                all_datasets = load_expression_atlas_datasets()
                logger.info(
                    f"✓ Successfully loaded Expression Atlas catalog: {len(all_datasets)} datasets available"
                )

                # Test a few specific datasets
                test_datasets = (
                    args.expression_atlas_datasets
                    if args.expression_atlas_datasets != ["all"]
                    else ["E-MTAB-5061"]
                )
                for accession in test_datasets[:3]:
                    dataset = next(
                        (ds for ds in all_datasets if ds.accession == accession), None
                    )
                    if dataset:
                        logger.info(f"✓ Dataset {accession} found in catalog")
                        results[accession] = ValidationResult(
                            dataset_id=accession,
                            dataset_type="expression_atlas",
                            status=ValidationStatus.TIER1_PASSED,
                            error_message=None,
                        )
                    else:
                        logger.warning(f"✗ Dataset {accession} not found in catalog")

            except Exception as e:
                logger.error(f"Failed to load Expression Atlas catalog: {e}")

        logger.info("Quick validation completed - no data downloads performed")

    else:
        # Standard validation with data downloads
        # Validate UCSC datasets
        if args.ucsc_datasets:
            if args.ucsc_datasets == ["all"]:
                logger.info("Validating all UCSC datasets")
                results.update(
                    validator.validate_ucsc_datasets_batch(
                        max_datasets=args.max_datasets
                    )
                )
            else:
                logger.info(f"Validating specific UCSC datasets: {args.ucsc_datasets}")
                results.update(
                    validator.validate_ucsc_datasets_batch(
                        dataset_names=args.ucsc_datasets, max_datasets=args.max_datasets
                    )
                )

        # Validate Expression Atlas datasets
        if args.expression_atlas_datasets:
            if args.expression_atlas_datasets == ["all"]:
                logger.info("Validating all Expression Atlas datasets")
                results.update(
                    validator.validate_expression_atlas_datasets_batch(
                        max_datasets=args.max_datasets
                    )
                )
            else:
                logger.info(
                    f"Validating specific Expression Atlas datasets: {args.expression_atlas_datasets}"
                )
                results.update(
                    validator.validate_expression_atlas_datasets_batch(
                        accessions=args.expression_atlas_datasets,
                        max_datasets=args.max_datasets,
                    )
                )

    # Print summary
    print(f"\nValidation completed for {len(results)} datasets")

    success_count = sum(
        1
        for r in results.values()
        if r.status
        in [
            ValidationStatus.TIER2_PASSED,
            ValidationStatus.TIER3_PASSED,
            ValidationStatus.TIER1_PASSED,
        ]
    )
    failed_count = len(results) - success_count

    print(f"Successful: {success_count}")
    print(f"Failed: {failed_count}")

    if getattr(args, "quick_validation", False):
        print("Note: Quick validation mode - only metadata and availability checked")

    if args.report_file:
        report = validator.generate_validation_report(args.report_file)
        print(f"Detailed report saved to: {args.report_file}")

    # Return non-zero exit code if there were failures and strict mode is enabled
    if args.strict and failed_count > 0:
        sys.exit(1)


def report_command(args):
    """Handle the report command."""
    setup_logging(args.verbose)

    validator = DatasetValidator(registry_file=args.registry_file)

    if args.output:
        report = validator.generate_validation_report(args.output)
        print(f"Report saved to: {args.output}")
    else:
        report = validator.generate_validation_report()
        print(report)


def list_datasets_command(args):
    """Handle the list-datasets command."""
    setup_logging(args.verbose)

    if args.source == "ucsc":
        from hvantk.datasets.ucsc_cell_datasets import load_ucsc_datasets

        datasets = load_ucsc_datasets()
        print(f"Available UCSC datasets ({len(datasets)}):")
        for dataset in datasets[: args.limit] if args.limit else datasets:
            print(f"  {dataset.name} - {dataset.shortLabel}")
            if args.details:
                print(
                    f"    Organisms: {', '.join(dataset.organisms) if dataset.organisms else 'N/A'}"
                )
                print(f"    Sample count: {dataset.sampleCount or 'N/A'}")

    elif args.source == "expression_atlas":
        from hvantk.datasets.expression_atlas_datasets import (
            load_expression_atlas_datasets,
        )

        datasets = load_expression_atlas_datasets()
        print(f"Available Expression Atlas datasets ({len(datasets)}):")
        for dataset in datasets[: args.limit] if args.limit else datasets:
            print(f"  {dataset.accession} - {dataset.title}")
            if args.details:
                print(f"    Type: {dataset.type}")
                print(f"    Files: {len(dataset.files) if dataset.files else 0}")

    else:  # both
        from hvantk.datasets.ucsc_cell_datasets import load_ucsc_datasets
        from hvantk.datasets.expression_atlas_datasets import (
            load_expression_atlas_datasets,
        )

        ucsc_datasets = load_ucsc_datasets()
        atlas_datasets = load_expression_atlas_datasets()

        print(f"Available datasets:")
        print(f"  UCSC: {len(ucsc_datasets)}")
        print(f"  Expression Atlas: {len(atlas_datasets)}")


def status_command(args):
    """Handle the status command."""
    setup_logging(args.verbose)

    validator = DatasetValidator(registry_file=args.registry_file)
    stats = validator.registry.get_summary_stats()

    if stats["total"] == 0:
        print("No validation results found.")
        return

    print(f"Validation Registry Status")
    print(f"=" * 30)
    print(f"Total datasets: {stats['total']}")
    print()

    print("By Status:")
    for status, count in stats.get("by_status", {}).items():
        print(f"  {status}: {count}")
    print()

    print("By Dataset Type:")
    for dtype, count in stats.get("by_dataset_type", {}).items():
        print(f"  {dtype}: {count}")

    if stats.get("by_failure_type"):
        print()
        print("By Failure Type:")
        for ftype, count in stats["by_failure_type"].items():
            print(f"  {ftype}: {count}")

    if args.show_successful:
        successful = validator.registry.list_successful_datasets()
        if successful:
            print(f"\nSuccessful datasets ({len(successful)}):")
            for dataset_id in successful[:10]:  # Show first 10
                result = validator.registry.get_result(dataset_id)
                print(f"  {dataset_id} ({result.dataset_type})")
            if len(successful) > 10:
                print(f"  ... and {len(successful) - 10} more")


def registry_command(args):
    """Handle the registry command for creation registry management."""
    setup_logging(args.verbose)

    # Initialize creation registry
    creation_registry = DatasetCreationRegistry(
        registry_file=args.creation_registry_file
    )

    if args.action == "list":
        # List datasets in creation registry
        datasets = creation_registry.list_datasets(
            source=DatasetSource(args.filter_source) if args.filter_source else None,
            matrix_type=(
                MatrixType(args.filter_matrix_type) if args.filter_matrix_type else None
            ),
            status=CreationStatus(args.filter_status) if args.filter_status else None,
        )

        print(f"Dataset Creation Registry ({len(datasets)} datasets)")
        print("=" * 50)

        for dataset in datasets:
            status_emoji = {
                "validated": "✅",
                "failed": "❌",
                "untested": "⚪",
                "deprecated": "🚫",
                "requires_update": "⚠️",
            }.get(dataset.status.value, "❓")

            print(f"{status_emoji} {dataset.dataset_id}")
            print(f"    Title: {dataset.metadata.title}")
            print(f"    Source: {dataset.source.value}")
            print(f"    Type: {dataset.matrix_type.value}")
            print(f"    Status: {dataset.status.value}")
            if dataset.metadata.sample_count:
                print(f"    Samples: {dataset.metadata.sample_count:,}")
            if args.show_commands:
                print(f"    Command: {dataset.creation_command}")
            print()

    elif args.action == "show":
        # Show detailed information for specific dataset
        if not args.dataset_id:
            print("Error: --dataset-id required for 'show' action")
            sys.exit(1)

        dataset = creation_registry.get_dataset(args.dataset_id)
        if not dataset:
            print(f"Dataset '{args.dataset_id}' not found in creation registry")
            sys.exit(1)

        # Print detailed dataset information
        print(f"Dataset: {dataset.dataset_id}")
        print("=" * 50)
        print(f"Title: {dataset.metadata.title}")
        print(f"Description: {dataset.metadata.description}")
        print(f"Source: {dataset.source.value}")
        print(f"Matrix Type: {dataset.matrix_type.value}")
        print(f"Status: {dataset.status.value}")
        print(f"Organism: {dataset.metadata.organism}")
        if dataset.metadata.tissue_type:
            print(f"Tissue: {dataset.metadata.tissue_type}")
        if dataset.metadata.sample_count:
            print(f"Sample Count: {dataset.metadata.sample_count:,}")
        print(f"Last Tested: {dataset.last_tested}")
        print()

        print("Creation Command:")
        print(f"  {dataset.creation_command}")
        print()

        print("Requirements:")
        print(f"  Memory: {dataset.requirements.min_memory_gb}+ GB")
        print(
            f"  Partitions: {dataset.requirements.min_partitions}-{dataset.requirements.max_partitions or 'unlimited'}"
        )
        print(f"  Tools: {', '.join(dataset.requirements.required_tools)}")
        if dataset.requirements.special_handling:
            print(f"  Special Handling: {dataset.requirements.special_handling}")
        print()

        if dataset.usage_examples:
            print("Usage Examples:")
            for i, example in enumerate(dataset.usage_examples, 1):
                print(f"  {i}. {example.title}")
                print(f"     {example.description}")
                print(f"     {example.code_snippet}")
                print()

        if dataset.known_issues:
            print("Known Issues:")
            for issue in dataset.known_issues:
                print(f"  • {issue}")
            print()

        if dataset.alternative_approaches:
            print("Alternative Approaches:")
            for approach in dataset.alternative_approaches:
                print(f"  • {approach}")

    elif args.action == "stats":
        # Show registry statistics
        stats = creation_registry.get_summary_stats()
        print("Dataset Creation Registry Statistics")
        print("=" * 40)
        print(f"Total Datasets: {stats['total']}")
        print()

        print("By Status:")
        for status, count in stats["by_status"].items():
            status_emoji = {
                "validated": "✅",
                "failed": "❌",
                "untested": "⚪",
                "deprecated": "🚫",
                "requires_update": "⚠️",
            }.get(status, "❓")
            print(f"  {status_emoji} {status}: {count}")
        print()

        print("By Source:")
        for source, count in stats["by_source"].items():
            print(f"  {source}: {count}")
        print()

        print("By Matrix Type:")
        for matrix_type, count in stats["by_matrix_type"].items():
            print(f"  {matrix_type}: {count}")

    elif args.action == "report":
        # Generate comprehensive report
        output_file = args.output if args.output else None
        report = creation_registry.generate_creation_report(output_file)
        if not output_file:
            print(report)

    elif args.action == "examples":
        # Show creation examples for dataset
        if not args.dataset_id:
            print("Error: --dataset-id required for 'examples' action")
            sys.exit(1)

        examples = creation_registry.get_creation_examples(args.dataset_id)
        if not examples:
            print(f"No examples found for dataset '{args.dataset_id}'")
            sys.exit(1)

        print(f"Creation Examples for {args.dataset_id}")
        print("=" * 50)
        for line in examples:
            print(line)


def create_parser():
    """Create the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="Dataset validation and creation registry management for HVANTK",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Validate specific UCSC datasets
  python -m hvantk.commands.dataset_validation_cli validate-datasets --ucsc-datasets cortex-dev zeisel2015

  # Validate all UCSC datasets (limit to 5 for testing)
  python -m hvantk.commands.dataset_validation_cli validate-datasets --ucsc-datasets all --max-datasets 5

  # Generate validation report
  python -m hvantk.commands.dataset_validation_cli report --output validation_report.txt

  # Check validation status
  python -m hvantk.commands.dataset_validation_cli status

  # List available datasets
  python -m hvantk.commands.dataset_validation_cli list-datasets --source ucsc --details
        """,
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Validate datasets command
    validate_parser = subparsers.add_parser(
        "validate-datasets", help="Validate datasets for matrix creation"
    )
    validate_parser.add_argument(
        "--ucsc-datasets",
        nargs="+",
        metavar="DATASET",
        help='UCSC datasets to validate (use "all" for all datasets)',
    )
    validate_parser.add_argument(
        "--expression-atlas-datasets",
        nargs="+",
        metavar="ACCESSION",
        help='Expression Atlas datasets to validate (use "all" for all datasets)',
    )
    validate_parser.add_argument(
        "--work-dir",
        default="dataset_validation",
        help="Working directory for downloads and samples (default: dataset_validation)",
    )
    validate_parser.add_argument(
        "--sample-lines",
        type=int,
        default=100,
        help="Number of lines for sample validation (default: 100)",
    )
    validate_parser.add_argument(
        "--max-datasets",
        type=int,
        help="Maximum number of datasets to validate (for testing)",
    )
    validate_parser.add_argument(
        "--registry-file", help="Path to validation registry file"
    )
    validate_parser.add_argument(
        "--report-file", help="Generate detailed report to file"
    )
    validate_parser.add_argument(
        "--strict",
        action="store_true",
        help="Exit with error code if any validation fails",
    )
    validate_parser.add_argument(
        "--force-revalidate",
        action="store_true",
        help="Force revalidation of datasets even if they were previously validated",
    )
    validate_parser.add_argument(
        "--quick-validation",
        action="store_true",
        help="Quick validation mode - only check metadata and headers, no full downloads",
    )
    validate_parser.add_argument(
        "--verbose", "-v", action="store_true", help="Enable verbose logging"
    )
    validate_parser.set_defaults(func=validate_datasets_command)

    # Report command
    report_parser = subparsers.add_parser("report", help="Generate validation report")
    report_parser.add_argument(
        "--output", "-o", help="Output file for report (default: print to stdout)"
    )
    report_parser.add_argument(
        "--registry-file", help="Path to validation registry file"
    )
    report_parser.set_defaults(func=report_command)

    # List datasets command
    list_parser = subparsers.add_parser("list-datasets", help="List available datasets")
    list_parser.add_argument(
        "--source",
        choices=["ucsc", "expression_atlas", "both"],
        default="both",
        help="Dataset source to list (default: both)",
    )
    list_parser.add_argument("--limit", type=int, help="Limit number of datasets shown")
    list_parser.add_argument(
        "--details", action="store_true", help="Show detailed information"
    )
    list_parser.set_defaults(func=list_datasets_command)

    # Status command
    status_parser = subparsers.add_parser(
        "status", help="Show validation registry status"
    )
    status_parser.add_argument(
        "--registry-file", help="Path to validation registry file"
    )
    status_parser.add_argument(
        "--show-successful",
        action="store_true",
        help="Show list of successful datasets",
    )
    status_parser.set_defaults(func=status_command)

    # Registry management command
    registry_parser = subparsers.add_parser(
        "registry", help="Manage dataset creation registry"
    )
    registry_parser.add_argument(
        "action", choices=["list", "show", "stats"], help="Registry action to perform"
    )
    registry_parser.add_argument("--dataset-id", help="Dataset ID for show action")
    registry_parser.add_argument(
        "--creation-registry-file", help="Path to creation registry file"
    )
    registry_parser.add_argument(
        "--filter-source",
        choices=["ucsc", "expression_atlas"],
        help="Filter by dataset source",
    )
    registry_parser.add_argument(
        "--filter-matrix-type",
        choices=[m.value for m in MatrixType],
        help="Filter by matrix type",
    )
    registry_parser.add_argument(
        "--filter-status",
        choices=["validated", "failed", "untested", "deprecated", "requires_update"],
        help="Filter by creation status",
    )
    registry_parser.add_argument(
        "--show-commands",
        action="store_true",
        help="Show creation commands in list view",
    )
    registry_parser.set_defaults(func=registry_command)

    return parser


def main():
    """Main CLI entry point."""
    parser = create_parser()
    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    try:
        args.func(args)
    except KeyboardInterrupt:
        print("\nOperation cancelled by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        if args.verbose:
            import traceback

            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
