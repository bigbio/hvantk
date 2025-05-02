import os
import shutil
import zipfile

import requests
from tqdm import tqdm
import os.path as path

import logging

logger = logging.getLogger(__name__)


def download_file(url: str, out_dir: str, file_name: str):
    """
    Download a file from a URL to a local directory.

    :param url: URL of the file to download
    :param out_dir: Local directory to save the file
    :param file_name: Name of the file to save

    :return: None
    """
    # Validate file_name doesn't contain path traversal
    if path.isabs(file_name) or ".." in file_name:
        raise ValueError(
            f"Invalid file_name: {file_name}. Must be a simple filename without path components."
        )

    os.makedirs(out_dir, exist_ok=True)
    local_path = os.path.join(out_dir, file_name)
    logger.info(f"Downloading {url} to {local_path}")
    try:
        response = requests.get(url, stream=True, timeout=30)  # Add timeout
        response.raise_for_status()  # Raises exception for 4XX/5XX responses

        total_size = int(response.headers.get("content-length", 0))

        # Optional: Add file size check
        max_size = 1024 * 1024 * 1024  # 1GB
        if total_size > max_size:
            logger.warning(
                f"File is very large ({total_size/1024/1024:.1f} MB), exceeding recommended size of {max_size/1024/1024:.1f} MB"
            )

        with open(local_path, "wb") as f, tqdm(
            desc=file_name,
            total=total_size,
            unit="B",
            unit_scale=True,
            unit_divisor=1024,
        ) as bar:
            for chunk in response.iter_content(chunk_size=8192):
                size = f.write(chunk)
                bar.update(size)
        logger.info(f"Download completed: {local_path}")

    except requests.exceptions.RequestException as e:
        # Clean up partial download if it exists
        if os.path.exists(local_path):
            os.remove(local_path)
        raise Exception(f"Failed to download {url}: {str(e)}")


# define a function to test if a given url exists
def url_exists(url: str) -> bool:
    """
    Check if a URL exists by sending a HEAD request.

    :param url: URL to check
    :type url: str
    :return: True if the URL exists, False otherwise
    :rtype: bool
    """
    try:
        response = requests.head(url, allow_redirects=True)
        return response.status_code == 200
    except requests.RequestException:
        return False


def compress_files(source_dir: str, output_zip: str, remove_originals: bool = False) -> None:
    """
    Compresses the files in a given source directory into a ZIP archive, maintaining
    the folder structure. Optionally, the original source directory can be removed
    after compression.

    :param source_dir: Path to the source directory to compress.
    :type source_dir: str
    :param output_zip: Path to the resulting ZIP file, including the file name and
        extension.
    :type output_zip: str
    :param remove_originals: Indicates whether the original source directory should
        be removed after compression. Defaults to False.
    :type remove_originals: bool
    :return: This function does not return any value.
    :rtype: None
    """
    with zipfile.ZipFile(output_zip, "w", zipfile.ZIP_DEFLATED) as zipf:
        for root, _, files in os.walk(source_dir):
            for file in files:
                file_path = os.path.join(root, file)
                # Use relative path in the archive to preserve folder structure.
                arcname = os.path.relpath(file_path, source_dir)
                zipf.write(file_path, arcname)
    logging.info(f"Compressed '{source_dir}' into '{output_zip}'")

    if remove_originals:
        shutil.rmtree(source_dir)
        logging.info(f"Removed original directory '{source_dir}' after compression.")


def decompress_files(zip_path: str, extract_to: str, remove_originals: bool = False) -> None:
    """
    Decompresses files from a zip archive to a specified location. Optionally, removes the original zip
    archive after extraction.

    :param zip_path: Path to the zip file to be decompressed.
    :type zip_path: str

    :param extract_to: Directory where the contents of the zip file will be extracted.
    :type extract_to: str

    :param remove_originals: Whether to delete the original zip file after extraction. Defaults to False.
    :type remove_originals: bool

    :return: None
    :raises FileNotFoundError: If the zip file does not exist
    :raises zipfile.BadZipFile: If the file is not a valid zip file
    :raises PermissionError: If there are permission issues with extracting files
    """
    if not os.path.exists(zip_path):
        raise FileNotFoundError(f"Zip file '{zip_path}' does not exist")

    # Create extract directory if it doesn't exist
    os.makedirs(extract_to, exist_ok=True)

    with zipfile.ZipFile(zip_path, "r") as zipf:
        # Check if the zip file is valid
        if zipf.testzip() is not None:
            raise zipfile.BadZipFile(f"'{zip_path}' is corrupted")

        zipf.extractall(extract_to)
    logging.info(f"Extracted '{zip_path}' into '{extract_to}'")

    if remove_originals:
        os.remove(zip_path)
        logging.info(f"Removed archive file '{zip_path}' after decompression.")
