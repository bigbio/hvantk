# HVANTK Dataset Registry Configuration
# This file contains configuration options for the dataset validation registry

class RegistryConfig:
    """Configuration class for the dataset validation registry."""

    def __init__(self):
        # Registry Settings
        self.validation_schedule = "0 2 * * *"  # Daily at 2 AM UTC

        # Validation limits
        self.max_datasets_per_run = 100
        self.sample_lines = 100
        self.timeout_minutes = 120

        # Notification settings
        self.create_issues_on_failure = True
        self.issue_labels = ["validation", "automation", "bug"]

        # Dataset source URLs (for reference)
        self.data_sources = {
            "ucsc": "https://cells.ucsc.edu/",
            "expression_atlas": "https://www.ebi.ac.uk/gxa/"
        }

# Validation tier descriptions
VALIDATION_TIERS = {
    "tier1": {
        "name": "Header Validation",
        "description": "Validates first 5 lines of files - checks file format, delimiters, column names",
        "sample_size": 5,
        "quick_check": True
    },
    "tier2": {
        "name": "Sample Validation",
        "description": "Tests matrix creation with small sample - identifies schema and format issues",
        "sample_size": 100,
        "creates_matrix": True
    },
    "tier3": {
        "name": "Full Validation",
        "description": "Complete dataset validation - full matrix creation test",
        "sample_size": None,
        "comprehensive": True
    }
}

# Status definitions
STATUS_DEFINITIONS = {
    "tier3_passed": {
        "level": "success",
        "description": "Full dataset validated successfully",
        "confidence": "high"
    },
    "tier2_passed": {
        "level": "warning",
        "description": "Sample validation successful",
        "confidence": "medium"
    },
    "tier1_passed": {
        "level": "info",
        "description": "Header validation successful",
        "confidence": "low"
    },
    "tier1_failed": {
        "level": "error",
        "description": "Header validation failed",
        "confidence": "none"
    },
    "tier2_failed": {
        "level": "error",
        "description": "Sample validation failed",
        "confidence": "none"
    },
    "tier3_failed": {
        "level": "error",
        "description": "Full validation failed",
        "confidence": "none"
    },
    "download_failed": {
        "level": "error",
        "description": "Unable to download dataset",
        "confidence": "none"
    },
    "not_tested": {
        "level": "unknown",
        "description": "Dataset not yet tested",
        "confidence": "none"
    }
}
