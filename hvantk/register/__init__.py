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
    DatasetValidationRegistry,
)

from .api_generator import APIEndpointGenerator
from .config import RegistryConfig, VALIDATION_TIERS, STATUS_DEFINITIONS

# Optional modules (guarded to avoid import-time failures)
try:
    from .web_generator import WebRegistryGenerator
except ImportError:
    WebRegistryGenerator = None

try:
    from .manager import RegistryManager
except ImportError:
    RegistryManager = None

# Build __all__ dynamically based on what was successfully imported
__all__ = [
    "ValidationStatus",
    "FailureType",
    "ValidationResult",
    "DatasetValidationRegistry",
    "APIEndpointGenerator",
    "RegistryConfig",
    "VALIDATION_TIERS",
    "STATUS_DEFINITIONS",
]

# Add optional imports only if they succeeded
if WebRegistryGenerator is not None:
    __all__.append("WebRegistryGenerator")

if RegistryManager is not None:
    __all__.append("RegistryManager")
