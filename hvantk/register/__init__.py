"""
HVANTK Dataset Validation Registry Module

This module provides comprehensive dataset validation registry functionality including:
- Live validation status tracking
- Web interface generation
- API endpoint creation
- GitHub Actions integration
- Multi-tier validation management

The registry system enables users to see which datasets are supported and their
validation status through a self-hosted GitHub Pages interface.
"""

from .validation_registry import (
    ValidationStatus,
    FailureType,
    ValidationResult,
    DatasetValidationRegistry
)

from .web_generator import WebRegistryGenerator
from .api_generator import APIEndpointGenerator
from .config import RegistryConfig, VALIDATION_TIERS, STATUS_DEFINITIONS
from .manager import RegistryManager

__all__ = [
    'ValidationStatus',
    'FailureType',
    'ValidationResult',
    'DatasetValidationRegistry',
    'WebRegistryGenerator',
    'APIEndpointGenerator',
    'RegistryConfig',
    'RegistryManager',
    'VALIDATION_TIERS',
    'STATUS_DEFINITIONS'
]
