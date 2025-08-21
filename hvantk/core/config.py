import os
import logging

logger = logging.getLogger(__name__)

# context settings for click
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])
logger.debug(f"Context settings: {CONTEXT_SETTINGS}")

# Global paths (may be configured at runtime)
RAW_DATA_PATH = None
ANNOTATION_DATA_PATH = None


def set_raw_data_path(raw_data_path: str):
    """
    Sets the global RAW_DATA_PATH variable to the specified directory path.

    Args:
        raw_data_path: Directory path to use for raw data.

    Returns:
        The updated RAW_DATA_PATH.

    Raises:
        ValueError: If the provided path is not a valid directory.
    """
    global RAW_DATA_PATH
    if os.path.isdir(raw_data_path):
        RAW_DATA_PATH = raw_data_path
        logger.info(f"Setting RAW_DATA_PATH to {raw_data_path}")
        return RAW_DATA_PATH
    else:
        logger.error(f"Invalid raw_data_path: {raw_data_path}")
        raise ValueError("Invalid raw_data_path: {}".format(raw_data_path))


def set_annotation_data_path(annotation_data_path: str):
    """
    Sets the global annotation data path if the provided directory exists.

    Args:
        annotation_data_path: Path to the annotation data directory.

    Returns:
        The updated annotation data path.

    Raises:
        ValueError: If the provided path is not a valid directory.
    """
    global ANNOTATION_DATA_PATH
    if os.path.isdir(annotation_data_path):
        ANNOTATION_DATA_PATH = annotation_data_path
        logger.info(f"Setting ANNOTATION_DATA_PATH to {annotation_data_path}")
        return ANNOTATION_DATA_PATH
    else:
        logger.error(f"Invalid annotation_data_path: {annotation_data_path}")
        raise ValueError(
            "Invalid annotation_data_path: {}".format(annotation_data_path)
        )


# A dictionary of annotation data paths (reserved for future use)
ANNOTATION_DATA_PATHS = {}
logger.debug(f"Annotation data paths: {ANNOTATION_DATA_PATHS}")
