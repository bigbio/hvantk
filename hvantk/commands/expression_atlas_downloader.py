"""
Downloader and transformer for Expression Atlas datasets.
"""

import click
import logging
import os
import sys
import ftplib
import json
import time
import random

from hvantk.datasets.expression_atlas_datasets import ExpressionAtlasDatasetCollection
from hvantk.core.constants import EXPRESSION_ATLAS_JSON_FILE_PATH

logger = logging.getLogger(__name__)


def _download_file_with_retry(ftp, remote_file, local_file_path, ftp_url=None, ftp_path=None, max_attempts=3, initial_delay=1):
    """
    Download a file with retry logic to handle connection issues.
    
    Args:
        ftp (ftplib.FTP): FTP connection
        remote_file (str): Name of remote file to download
        local_file_path (str): Path where to save the downloaded file
        ftp_url (str, optional): FTP server URL for reconnection
        ftp_path (str, optional): FTP directory path for reconnection
        max_attempts (int): Maximum number of retry attempts
        initial_delay (int): Initial delay between retries in seconds
    
    Returns:
        bool: True if successful, False otherwise
    """
    attempt = 0
    delay = initial_delay

    while attempt < max_attempts:
        try:
            with open(local_file_path, "wb") as f:
                ftp.retrbinary(f"RETR {remote_file}", f.write)
            logger.info(f"Downloaded {remote_file} to {local_file_path}")
            return True  # Success
        except (ConnectionResetError, ftplib.error_temp, ftplib.error_proto, EOFError) as e:
            attempt += 1
            if attempt >= max_attempts:
                logger.error(f"Failed to download {remote_file} after {max_attempts} attempts: {str(e)}")
                return False

            # Log the retry attempt
            logger.warning(
                f"FTP connection error while downloading {remote_file}: {str(e)}. "
                f"Retrying ({attempt}/{max_attempts}) in {delay} seconds..."
            )

            # Wait before retrying with exponential backoff and jitter
            time.sleep(delay + random.uniform(0, 1))
            delay *= 2  # Exponential backoff

            # Re-establish FTP connection if URL is provided
            if ftp_url:
                import contextlib
                with contextlib.suppress(Exception):
                    ftp.quit()

                try:
                    ftp.connect(ftp_url)
                    ftp.login()
                    ftp.set_pasv(True)
                except (ftplib.error_temp, ConnectionRefusedError) as e:
                    logger.error(f"Failed to reconnect to FTP server: {str(e)}")
                    return False

                # Navigate back to the correct directory if path is provided
                if ftp_path:
                    try:
                        ftp.cwd(ftp_path)
                    except ftplib.error_perm as e:
                        logger.error(f"Cannot access {ftp_path}: {e}")
                        return False

    return False  # If we get here, all attempts failed


def _print_dataset_accessions():
    """
    Prints the accessions of all available Expression Atlas datasets in a formatted list.

    :raises FileNotFoundError: If the JSON file containing Expression Atlas dataset
                                information is not found at the specified path.
    :raises JSONDecodeError: If the JSON file contains invalid JSON format.
    """
    try:
        collection = ExpressionAtlasDatasetCollection.from_json(EXPRESSION_ATLAS_JSON_FILE_PATH)
        dataset_accessions = collection.list_dataset_accessions()
        if not dataset_accessions:
            click.echo("No datasets available.")
        else:
            click.echo(f"Available datasets ({len(dataset_accessions)}):")
            for accession in dataset_accessions:
                click.echo(f"- {accession}")
    except ValueError as e:
        click.echo(f"Error: {e}")


@click.command("expression-atlas-downloader", short_help="Download Expression Atlas dataset")
@click.option(
    "--config_path",
    required=False,
    help="The path to the registry config file (optional). Either --config_path or --accession must be supplied.",
)
@click.option(
    "--accession", required=False, help="The accession of the experiment to download."
)
@click.option(
    "--download_path",
    type=click.Path(exists=False, file_okay=False, dir_okay=True),
    required=True,
    help="The path to download the data to."
)
@click.option(
    "--list_datasets",
    is_flag=True,
    required=False,
    help="List all available datasets in the Expression Atlas collection.",
)
def download_experiments(config_path, accession, download_path, list_datasets):
    """
    This function downloads experiments from the Expression Atlas database using FTP.

    :param config_path: The optional path to an expression_atlas.json file. This JSON
        configuration file contains details on specific files to download for each
        experiment.
    :type config_path: str, optional
    :param accession: The optional accession of the experiment to download. If provided,
        only files from this specific experiment are downloaded unless filtered further
        via config_path.
    :type accession: str, optional
    :param download_path: The required path to save downloaded data. This directory will
        be created if it does not already exist.
    :type download_path: str
    :return: None
    :param list_datasets: If set, lists all available datasets in the Expression Atlas
        collection instead of downloading.
    :type list_datasets: bool
    :return: None
    :raises ftplib.error_perm: Raised if there are FTP permission issues while accessing
        the FTP server or its directories.
    :raises Exception: Raised for any other errors encountered during the operation.
    """

    if list_datasets:
        logger.warning("--accession and --config_path options are ignored when --list_datasets is specified.")
        _print_dataset_accessions()
        return

    ftp_url = "ftp.ebi.ac.uk"

    # Create download directory if it doesn't exist
    os.makedirs(download_path, exist_ok=True)

    try:
        with ftplib.FTP(ftp_url, timeout=60) as ftp:  # Increased timeout
            ftp.login()
            ftp.set_pasv(True)  # Use passive mode to help with firewalls

            if accession:
                experiment_id = accession
                ftp_path = (
                    f"/pub/databases/microarray/data/atlas/experiments/{experiment_id}/"
                )
                logger.info(
                    f"Downloading experiment {experiment_id} from {ftp_url}{ftp_path}"
                )
                ftp.cwd(ftp_path)
                files = ftp.nlst()

                if config_path:
                    with open(config_path, "r") as f:
                        config = json.load(f)
                    experiment_config = next(
                        (item for item in config if item["accession"] == accession),
                        None,
                    )
                    if experiment_config:
                        files_to_download = [
                            f["name"] for f in experiment_config["files"]
                        ]
                        logger.info(
                            f"Downloading only the files specified in the config file: {files_to_download}"
                        )
                        files = [f for f in files if f in files_to_download]
                    else:
                        logger.warning(
                            f"No configuration found for accession {accession} in {config_path}. Downloading all files."
                        )

                else:
                    logger.info("Downloading all files from the given accession.")

                # Download each file with retry logic
                for file in files:
                    file_path = os.path.join(download_path, file)
                    _download_file_with_retry(ftp, file, file_path, ftp_url, ftp_path)

            elif config_path:
                with open(config_path, "r") as f:
                    config = json.load(f)

                for experiment in config:
                    experiment_id = experiment["accession"]
                    ftp_path = f"/pub/databases/microarray/data/atlas/experiments/{experiment_id}/"
                    logger.info(
                        f"Downloading experiment {experiment_id} from {ftp_url}{ftp_path}"
                    )
                    try:
                        ftp.cwd(ftp_path)
                    except ftplib.error_perm as e:
                        logger.error(f"Cannot access {ftp_path}: {e}")
                        continue
                        
                    files = ftp.nlst()
                    experiment_config = next(
                        (item for item in config if item["accession"] == experiment_id),
                        None,
                    )
                    if experiment_config:
                        files_to_download = [
                            f["name"] for f in experiment_config["files"]
                        ]
                        logger.info(
                            f"Downloading only the files specified in the config file: {files_to_download}"
                        )
                        files = [f for f in files if f in files_to_download]
                    else:
                        logger.warning(
                            f"No configuration found for accession {experiment_id} in {config_path}. Downloading all files."
                        )

                    # Download each file with retry logic
                    for file in files:
                        file_path = os.path.join(download_path, file)
                        _download_file_with_retry(ftp, file, file_path, ftp_url, ftp_path)
            else:
                logger.error("Either --config_path or --accession must be provided.")
                sys.exit(1)

    except ftplib.error_perm as e:
        logger.error(f"FTP permission error: {e}")
        sys.exit(1)
    except Exception as e:
        logger.exception(f"Error downloading experiments: {e}")
        sys.exit(1)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    download_experiments()