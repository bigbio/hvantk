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
    "--config_path", required=True, help="The path to the expression_atlas.json config file."
)
@click.option(
    "--download_path", required=True, help="The path to download the data to."
)
def download_experiments(config_path, download_path):
    """Downloads experiment datasets from Expression Atlas using FTP and a config file."""
    ftp_url = f"ftp.ebi.ac.uk"

    # Create download directory if it doesn't exist
    if not os.path.exists(download_path):
        os.makedirs(download_path)

    try:
        with open(config_path, "r") as f:
            config = json.load(f)

        for experiment in config:
            experiment_id = experiment["accession"]
            ftp_path = f"/pub/databases/microarray/data/atlas/experiments/{experiment_id}/"
            logger.info(f"Downloading experiment {experiment_id} from {ftp_url}{ftp_path}")

            ftp = ftplib.FTP(ftp_url)
            ftp.login()
            ftp.cwd(ftp_path)

            # Get the list of files in the FTP directory
            files = ftp.nlst()

            # Download each file
            for file in files:
                file_path = os.path.join(download_path, file)
                with open(file_path, "wb") as f:
                    ftp.retrbinary(f"RETR {file}", f.write)
                logger.info(f"Downloaded {file} to {file_path}")

            ftp.quit()

    except Exception as e:
        logger.error(f"Error downloading experiments: {e}")
        sys.exit(1)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    download_experiments()