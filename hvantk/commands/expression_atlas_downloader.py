"""
Downloader and transformer for Expression Atlas datasets.
"""

import click
import logging
import os
import sys
import ftplib
import json

from hvantk.datasets.expression_atlas_datasets import ExpressionAtlasDatasetCollection
from hvantk.utils.constants import EXPRESSION_ATLAS_JSON_FILE_PATH

logger = logging.getLogger(__name__)

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
            click.echo("Available datasets:")
            for accession in dataset_accessions:
                click.echo(f"- {accession}")
    except ValueError as e:
        click.echo(f"Error: {e}")


@click.command("expression-atlas-downloader", short_help="Download Expression Atlas dataset")
@click.option(
    "--config_path",
    required=False,
    help="The path to the expression_atlas.json config file.",
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
        _print_dataset_accessions()
        return

    ftp_url = "ftp.ebi.ac.uk"

    # Create download directory if it doesn't exist
    os.makedirs(download_path, exist_ok=True)

    try:
        with ftplib.FTP(ftp_url) as ftp:
            ftp.login()

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

                # Download each file
                for file in files:
                    file_path = os.path.join(download_path, file)
                    with open(file_path, "wb") as f:
                        ftp.retrbinary(f"RETR {file}", f.write)
                    logger.info(f"Downloaded {file} to {file_path}")

            elif config_path:
                with open(config_path, "r") as f:
                    config = json.load(f)

                for experiment in config:
                    experiment_id = experiment["accession"]
                    ftp_path = f"/pub/databases/microarray/data/atlas/experiments/{experiment_id}/"
                    logger.info(
                        f"Downloading experiment {experiment_id} from {ftp_url}{ftp_path}"
                    )
                    ftp.cwd(ftp_path)
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

                    # Download each file
                    for file in files:
                        file_path = os.path.join(download_path, file)
                        with open(file_path, "wb") as f:
                            ftp.retrbinary(f"RETR {file}", f.write)
                        logger.info(f"Downloaded {file} to {file_path}")
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
