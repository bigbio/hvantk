"""
Downloader and transformer for Expression Atlas datasets.
"""

import click
import logging
import os
import sys
import ftplib
import json

logger = logging.getLogger(__name__)


@click.command()
@click.option(
    "--config_path",
    required=False,
    help="The path to the expression_atlas.json config file.",
)
@click.option(
    "--accession", required=False, help="The accession of the experiment to download."
)
@click.option(
    "--download_path", required=True, help="The path to download the data to."
)
def download_experiments(config_path, accession, download_path):
    """Downloads experiment datasets from Expression Atlas using FTP."""
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
