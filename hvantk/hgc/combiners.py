import os
import logging
from typing import List, Dict, Any

import hail as hl
from hvantk.hgc.file_utils import sort_mts_cols
from hvantk.hgc.file_utils import validate_vcfs_paths, validate_vds_paths
from hvantk.hgc.constants import HG38_GENOME_REFERENCE


def combine_gvcfs(
    gvcf_dir: str,
    vds_output_path: str,
    tmp_path: str,
    save_path: str,
    vdses: List[str],
    kwargs: Dict[str, Any],
    reference_genome: str = HG38_GENOME_REFERENCE,
) -> None:
    """
    Combines GVCF files into a VDS using Hail's GVCF combiner.

    This function initializes Hail with the specified configuration, validates the input GVCF
    file paths using `get_vcfs_files`, and runs Hail's GVCF combiner to produce a variant
    dataset (VDS). The resulting VDS is written to the specified output path.

    Parameters:
        gvcf_dir (str): Directory containing GVCF files.
        vds_output_path (str): Path where the output VDS will be written.
        tmp_path (str): Temporary directory path used during processing.
        save_path (str): Path to save the combiner plan.
        vdses (List[str]): List of VDS paths to be combined.
        kwargs (dict): Additional keyword arguments to pass to the combiner.
        reference_genome (str): Reference genome to use (default: GRCh38).

    Returns:
        None

    Raises:
        FileNotFoundError, ValueError, PermissionError: If file validations fail.
    """
    try:
        if not (gvcf_dir or vdses):
            raise ValueError("Either GVCF files or VDS files must be provided.")

        validated_gvcfs: List[str] = []
        validated_vdses: List[str] = []

        if gvcf_dir:
            # Validate and retrieve GVCF file paths.
            validated_gvcfs = validate_vcfs_paths(gvcf_dir)
            logging.info(f"Validated {len(validated_gvcfs)} GVCF file paths.")

        if vdses:
            # Validate and retrieve VDS file paths.
            validated_vdses = validate_vds_paths(vdses)
            logging.info(f"Validated {len(validated_vdses)} VDS file paths.")

        # Set up and run the Hail GVCF combiner
        combiner = hl.vds.new_combiner(
            output_path=vds_output_path,
            temp_path=tmp_path,
            gvcf_paths=validated_gvcfs,
            save_path=save_path,
            vds_paths=validated_vdses,
            use_genome_default_intervals=True,  # TODO: Use default intervals for the genome, but must be customized
            reference_genome=reference_genome,
            **kwargs,
        )
        combiner.run()

    except Exception as e:
        logging.exception("An error occurred during GVCF combination.")
        raise

    finally:
        logging.info("GVCF combination completed successfully.")


def combine_vdses(
    vdses_dir: str, output_path: str, validate: bool = True, overwrite: bool = False
) -> None:
    """
    Combine multiple VDS directories into a single VDS.

    This function takes a container directory that holds VDS directories (each VDS is
    represented as a directory with a specific extension), check every path is valid, combines their
    rows using Hail's union_rows, optionally validates the resulting VDS, and writes the result
    to the specified output path.

    Parameters:
        vdses_dir (str): Container directory containing VDS directories.
        output_path (str): Path where the combined VDS will be written.
        validate (bool): Whether to validate the combined VDS using Hail's validator.
        overwrite (bool): Whether to overwrite the output if it already exists.

    Raises:
        ValueError: If no valid VDS directories are found.
        Exception: If an error occurs during VDS combination or writing.
    """
    try:
        # Ensure a container directory is provided and valid.
        if not vdses_dir:
            raise ValueError("No VDS directory provided.")
        if not os.path.isdir(vdses_dir):
            raise NotADirectoryError(f"Provided path '{vdses_dir}' is not a directory.")

        # Validate VDS directories within the container.
        vdses_paths: List[str] = validate_vds_paths(vdses_dir)
        if not vdses_paths:
            raise ValueError(f"No valid VDS directories found in '{vdses_dir}'.")
        logging.info(f"Validated {len(vdses_paths)} VDS directories.")

        # Combine the VDS directories using Hail's union_rows function.
        vdses = [(hl.vds.read_vds(path)) for path in vdses_paths]
        combined_vds = hl.vds.VariantDataset.union_rows(*vdses)
        logging.info("Successfully combined VDS directories.")

        # Optionally validate the combined VDS.
        if validate:
            hl.vds.VariantDataset.validate(combined_vds)
            logging.info("Combined VDS validated successfully.")

        # Write the combined VDS to the output path.
        combined_vds.write(output_path, overwrite=overwrite)

    except Exception as e:
        logging.exception("An error occurred during VDS combination.")
        raise

    finally:
        logging.info(f"Combined VDS written to {output_path}")


def combine_matrix_table_rows(
    mt_paths: List[str],
    output_path: str,
    n_partitions: int,
    force_sort_cols: bool = False,
    overwrite: bool = False,
    kwargs=None,
) -> None:
    """
    Combine multiple MatrixTables by rows into a single MatrixTable.

    This function takes a list of MatrixTable paths, reads each MatrixTable, combines their rows (e.g., variants)
    using Hail's union_rows, optionally validates the resulting MatrixTable, and writes the result
    to the specified output path.

    Parameters:
        mt_paths (List[str]): List of MatrixTable paths to combine.
        output_path (str): Path where the combined MatrixTable will be written.
        n_partitions: Number of partitions to use when writing the combined MatrixTable.
        force_sort_cols (bool): Whether to sort the columns of the MatrixTables before combining.
        overwrite (bool): Whether to overwrite the output if it already exists.
        kwargs: Additional keyword arguments to pass to the Hail's union_rows function.

    Raises:
        ValueError: If no valid MatrixTables are found.
        Exception: If an error occurs during MatrixTable combination or writing.
    """
    if kwargs is None:
        kwargs = {}
    try:
        # Combine the MatrixTables using Hail's union_rows function.
        mts = [hl.read_matrix_table(path) for path in mt_paths]

        if force_sort_cols:
            logging.info("Sorting columns of MatrixTables before combining.")
            mts = sort_mts_cols(mts)

        logging.info("Combining MatrixTables rows...")
        combined_mt = hl.MatrixTable.union_rows(*mts, **kwargs)

        # Repartition the combined MatrixTable and write it to the output path.
        if n_partitions:
            logging.info(
                f"Repartitioning the combined MatrixTable into {n_partitions} partitions."
            )
            combined_mt = combined_mt.repartition(n_partitions)

        # Write the combined MatrixTable to the output path.
        combined_mt.write(output_path, overwrite=overwrite)

    except Exception as e:
        logging.exception("An error occurred during MatrixTable combination.")
        raise

    finally:
        logging.info(f"Combined MatrixTable written to {output_path}")


def combine_matrix_table_cols(
    mt_paths: List[str],
    output_path: str,
    n_partitions: int,
    overwrite: bool = False,
    kwargs=None,
) -> None:
    """
    Combine multiple MatrixTables by columns into a single MatrixTable.

    This function takes a list of MatrixTable paths, reads each MatrixTable, combines their columns (e.g., samples)
    using Hail's union_cols, optionally validates the resulting MatrixTable, and writes the result
    to the specified output path.

    Parameters:
        mt_paths (List[str]): List of MatrixTable paths to combine.
        output_path (str): Path where the combined MatrixTable will be written.
        n_partitions (bool): Number of partitions to use when writing the combined MatrixTable.
        overwrite (bool): Whether to overwrite the output if it already exists.
        kwargs: Additional keyword arguments to pass to the Hail's union_cols function.

    Raises:
        ValueError: If no valid MatrixTables are found.
        Exception: If an error occurs during MatrixTable combination or writing.
    """
    if kwargs is None:
        kwargs = {}
    try:
        # Combine the MatrixTables using Hail's union_cols function.
        mts = [hl.read_matrix_table(path) for path in mt_paths]

        logging.info("Combining MatrixTables columns...")
        combined_mt = hl.MatrixTable.union_cols(*mts, **kwargs)

        # Repartition the combined MatrixTable and write it to the output path.
        if n_partitions:
            logging.info(
                f"Repartitioning the combined MatrixTable into {n_partitions} partitions."
            )
            combined_mt = combined_mt.repartition(n_partitions)

        # Write the combined MatrixTable to the output path.
        combined_mt.write(output_path, overwrite=overwrite)

    except Exception as e:
        logging.exception("An error occurred during MatrixTable combination.")
        raise

    finally:
        logging.info(f"Combined MatrixTable written to {output_path}")
