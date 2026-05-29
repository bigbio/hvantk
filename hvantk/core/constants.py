# Core constants for the project

from pathlib import Path
import logging

logger = logging.getLogger(__name__)

# Base directory of the package (hvantk)
BASE_DIR = Path(__file__).resolve().parent
logger.debug(f"Base directory: {BASE_DIR}")

# Explicit public API for this module
__all__ = [
    "BASE_DIR",
]
