import glob
import os
import logging

import hail as hl

from typing import Union, List

from hvantk.algorithms.hgc.constants import GVCF_EXTENSION, GVCF_EXTENSION_TBI, VDS_EXTENSION

"""
Utility functions for file handling.
"""

# Get module logger - do not configure logging at import time
logger = logging.getLogger(__name__)


def check_path_exists_and_readable(path: str) -> str:
    """
    Check if a file or directory exists and is readable.

    Parameters:
        path (str): Path to the file or directory (absolute or relative).

    Returns:
        str: The original path if it exists and is readable.

    Raises:
        FileNotFoundError: If the file or directory does not exist or the path does not point to a valid file or directory.
        PermissionError: If the file or directory exists but is not readable.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"File or directory '{path}' not found.")
    if not (os.path.isfile(path) or os.path.isdir(path)):
        raise FileNotFoundError(f"Path '{path}' is not a valid file or directory.")
    if not os.access(path, os.R_OK):
        raise PermissionError(f"File or directory '{path}' is not readable.")

    return path


def validate_vcfs_paths(directory: str, pattern: str = None) -> List[str]:
    """
    Retrieve and validate a list of absolute paths to GVCF files in the given directory.

    This function performs the following:
      - Checks that the provided directory exists.
      - Uses glob to find files in the directory matching the given pattern.
      - For each discovered GVCF file:
          * Converts it to an absolute path.
          * Validates that it ends with the expected GVCF extension.
          * Checks that the GVCF file exists and is readable.
          * Checks that the corresponding TBI file (with the expected TBI extension)
            exists and is readable.

    Parameters:
        directory (str): The directory in which to search for GVCF files.
        pattern (str): The glob pattern to match files. If None, defaults to "*{GVCF_EXTENSION}".

    Returns:
        List[str]: A list of validated absolute paths to GVCF files.

    Raises:
        NotADirectoryError: If the provided directory does not exist.
        ValueError: If a file does not have the expected GVCF extension.
        FileNotFoundError or PermissionError: Propagated from the file existence/readability checks.
    """
    if pattern is None:
        pattern = f"*{GVCF_EXTENSION}"

    if not os.path.isdir(directory):
        raise NotADirectoryError(f"'{directory}' is not a valid directory.")

    search_pattern = os.path.join(directory, pattern)
    matching_files = glob.glob(search_pattern)

    valid_paths = []
    for filepath in matching_files:
        abs_path = os.path.abspath(filepath)

        # Validate the file extension using the GVCF_EXTENSION variable.
        if not abs_path.endswith(GVCF_EXTENSION):
            raise ValueError(f"File '{abs_path}' does not end with '{GVCF_EXTENSION}'.")

        # Validate that the GVCF file exists and is readable.
        check_path_exists_and_readable(abs_path)

        # Compute the corresponding TBI file path by replacing the GVCF extension with the TBI extension.
        base = abs_path[: -len(GVCF_EXTENSION)]
        tbi_path = base + GVCF_EXTENSION_TBI

        # Validate that the TBI file exists and is readable.
        check_path_exists_and_readable(tbi_path)

        valid_paths.append(abs_path)

    return valid_paths


def validate_vds_paths(vdses: Union[str, List[str]]) -> List[str]:
    """
    Validate VDS directories given either as a container directory or a list of VDS paths.

    If a directory is provided, the function checks that it exists and is readable, then iterates
    over all entries in the container. For each entry that is a directory, it verifies the directory
    is accessible and its name ends with the expected VDS_EXTENSION.

    If a list of paths is provided, each path is validated individually in a similar manner.

    Parameters:
        vdses (Union[str, List[str]]): Either a single directory path that contains VDS directories,
                                       or a list of individual VDS directory paths.

    Returns:
        List[str]: List of validated VDS directory paths.

    Raises:
        NotADirectoryError: If a provided path expected to be a directory is not one.
        FileNotFoundError or PermissionError: If a directory does not exist or is not accessible.
        ValueError: If a directory does not end with the expected VDS_EXTENSION.
    """
    validated_paths = []

    if isinstance(vdses, str):
        # vdses is a container directory
        if not os.path.isdir(vdses):
            raise NotADirectoryError(f"'{vdses}' is not a directory.")
        for entry in os.listdir(vdses):
            full_path = os.path.join(vdses, entry)
            if os.path.isdir(full_path):
                check_path_exists_and_readable(full_path)
                if not entry.endswith(VDS_EXTENSION):
                    raise ValueError(
                        f"Directory '{full_path}' does not end with '{VDS_EXTENSION}'."
                    )
                validated_paths.append(full_path)
    elif isinstance(vdses, list):
        # vdses is a list of directory paths
        for path in vdses:
            if not os.path.isdir(path):
                raise NotADirectoryError(f"'{path}' is not a directory.")
            check_path_exists_and_readable(path)
            base_name = os.path.basename(path)
            if not base_name.endswith(VDS_EXTENSION):
                raise ValueError(
                    f"Directory '{path}' does not end with '{VDS_EXTENSION}'."
                )
            validated_paths.append(path)
    else:
        raise TypeError(
            "Input must be either a directory path (str) or a list of directory paths."
        )

    return validated_paths


def sort_mts_cols(
    mts: List[hl.MatrixTable], ref_index: int = 0
) -> List[hl.MatrixTable]:
    """
    Sort the column order of a list of matrix tables to match a reference matrix table.

    The column order of each matrix table is updated to match the order defined by the reference matrix table.
    All matrix tables are assumed to have equal number of columns and the same column keys.

    Parameters:
        mts (List[hl.MatrixTable]): List of matrix table objects.
        ref_index (int): The index of the reference matrix table in the list. The reference table's columns
                         remain unchanged. Defaults to 0.

    Returns:
        List[hl.MatrixTable]: A new list of matrix table objects with columns reordered to match the reference.

    Raises:
        IndexError: If ref_index is out of range for the list of matrix tables.

    Example:
        sorted_mts = sort_mts_cols([mt1, mt2, mt3], ref_index=0)
    """
    if not 0 <= ref_index < len(mts):
        raise IndexError(
            "ref_index is out of range for the provided list of matrix tables."
        )

    # Compute the column order based on the reference matrix table.
    ref_mt = mts[ref_index].add_col_index()

    sorted_mts = []
    for i, mt in enumerate(mts):
        if i == ref_index:
            # Leave the reference matrix table unchanged (original without col_idx).
            sorted_mts.append(mt)
        else:
            mt_indexed = mt.add_col_index()
            new_order = mt_indexed.index_cols(ref_mt.col_key).col_idx.collect()
            # Reorder columns and drop the transient col_idx field before returning
            sorted_mts.append(mt_indexed.choose_cols(new_order).drop("col_idx"))

    return sorted_mts
