"""
HGC CLI Utilities Module

Shared utility functions for HGC CLI commands.
"""

import logging
import glob
import os
import tempfile

from hvantk.algorithms.hgc import check_path_exists_and_readable, validate_vds_paths

logger = logging.getLogger(__name__)

DEFAULT_TEMP_DIR = os.environ.get("HGC_TEMP_DIR", tempfile.gettempdir())


def setup_logging_for_hgc(log_level="INFO"):
    """Set up logging for HGC operations.

    .. note::
        Root logging (handlers/formatters) is configured centrally in
        ``hvantk.hvantk.setup_logging``. This helper applies the requested
        log level to the root logger.
    """
    if isinstance(log_level, str):
        level = getattr(logging, log_level.upper(), None)
        if level is None:
            logger.warning(
                "Unknown log level '%s' requested for HGC; defaulting to INFO.",
                log_level,
            )
            level = logging.INFO
    else:
        level = log_level

    logging.getLogger().setLevel(level)


def expand_file_patterns(patterns):
    """Expand file patterns/wildcards into actual file paths."""
    expanded = []
    for pattern in patterns:
        matches = glob.glob(pattern)
        expanded.extend(matches)
    return expanded


def validate_input_files(file_paths, file_type="gvcf"):
    """Validate input files and return (is_valid, errors) tuple."""
    errors = []
    try:
        if file_type == "gvcf":
            for path in file_paths:
                check_path_exists_and_readable(path)
        elif file_type == "vds":
            validate_vds_paths(file_paths)
        else:
            for path in file_paths:
                check_path_exists_and_readable(path)
        return (True, [])
    except Exception as e:
        errors.append(str(e))
        return (False, errors)


def validate_output_path(output_path, create_dirs=False):
    """Validate output path and optionally create parent directories."""
    try:
        parent_dir = os.path.dirname(output_path)
        if parent_dir and not os.path.exists(parent_dir):
            if create_dirs:
                os.makedirs(parent_dir, exist_ok=True)
            else:
                return False
        return True
    except Exception:
        return False


def estimate_resource_requirements(file_paths):
    """Estimate resource requirements for an operation."""
    logger = logging.getLogger(__name__)
    total_size = 0
    for path in file_paths:
        try:
            if os.path.isfile(path):
                total_size += os.path.getsize(path)
            elif os.path.isdir(path):
                for root, dirs, files in os.walk(path):
                    for file in files:
                        file_path = os.path.join(root, file)
                        try:
                            total_size += os.path.getsize(file_path)
                        except OSError as e:
                            logger.warning(
                                f"Failed to get size of file '{file_path}': {e}"
                            )
                            continue
        except OSError as e:
            logger.warning(f"Failed to access path '{path}': {e}")
            continue
        except Exception as e:
            logger.error(f"Unexpected error accessing path '{path}': {e}")
            raise

    total_size_gb = total_size / (1024**3)

    # Simple heuristic estimates
    memory = f"{max(4, int(total_size_gb * 2))}g"
    partitions = max(100, int(total_size_gb * 10))
    estimated_runtime = max(5, int(total_size_gb * 2))

    return {
        "memory": memory,
        "partitions": partitions,
        "estimated_runtime_minutes": estimated_runtime,
        "total_size_gb": round(total_size_gb, 2),
    }
