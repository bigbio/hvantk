import os
import shutil
import zipfile

import requests
from tqdm import tqdm
import os.path as path

import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)  # Set default log level to DEBUG


def download_file(url: str, out_dir: str, file_name: str):
    """
    Downloads a file from the specified URL to a local directory with a given file name.
    
    Validates that the file name does not contain path traversal components. Creates the output directory if it does not exist. Raises an exception if the download fails or if the file name is invalid.
    """
    # Validate file_name doesn't contain path traversal
    if path.isabs(file_name) or ".." in file_name:
        raise ValueError(
            f"Invalid file_name: {file_name}. Must be a simple filename without path components."
        )

    os.makedirs(out_dir, exist_ok=True)
    local_path = os.path.join(out_dir, file_name)
    logger.debug(f"Downloading {url} to {local_path}")
    try:
        response = requests.get(url, stream=True, timeout=30)  # Add timeout
        response.raise_for_status()  # Raises exception for 4XX/5XX responses
        logger.debug(f"Download response status code: {response.status_code}")

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
        logger.exception(f"Failed to download {url}: {str(e)}")
        raise Exception(f"Failed to download {url}: {str(e)}")


# define a function to test if a given url exists
def url_exists(url: str) -> bool:
    """
    Checks whether the specified URL exists by performing an HTTP HEAD request.
    
    Returns:
        True if the URL responds with status code 200; otherwise, False.
    """
    try:
        logger.debug(f"Checking if URL exists: {url}")
        response = requests.head(url, allow_redirects=True)
        if response.status_code == 200:
            logger.debug(f"URL {url} exists with status code: {response.status_code}")
            return True
        else:
            logger.warning(
                f"URL {url} does not exist with status code: {response.status_code}"
            )
            return False
    except requests.RequestException as e:
        logger.exception(f"Error checking URL {url}: {e}")
        return False


def compress_files(
    source_dir: str, output_zip: str, remove_originals: bool = False
) -> None:
    """
    Compresses all files in a directory into a ZIP archive, preserving folder structure.
    
    Args:
        source_dir: Path to the directory whose contents will be compressed.
        output_zip: Path to the output ZIP file.
        remove_originals: If True, deletes the source directory after compression.
    
    Raises:
        FileNotFoundError: If the source directory does not exist.
        ValueError: If the source path is not a directory.
        PermissionError: If there are file access permission issues.
    """
    if not os.path.exists(source_dir):
        raise FileNotFoundError(f"Source directory '{source_dir}' does not exist")

    if not os.path.isdir(source_dir):
        raise ValueError(f"Source path '{source_dir}' is not a directory")

    # Create parent directory for output_zip if it doesn't exist
    output_dir = os.path.dirname(output_zip)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    with zipfile.ZipFile(output_zip, "w", zipfile.ZIP_DEFLATED) as zipf:
        for root, _, files in os.walk(source_dir):
            for file in files:
                file_path = os.path.join(root, file)
                # Use relative path in the archive to preserve folder structure.
                arcname = os.path.relpath(file_path, source_dir)
                zipf.write(file_path, arcname)
    logger.info(f"Compressed '{source_dir}' into '{output_zip}'")

    if remove_originals:
        shutil.rmtree(source_dir)
        logger.info(f"Removed original directory '{source_dir}' after compression.")


def decompress_files(
    zip_path: str, extract_to: str, remove_originals: bool = False
) -> None:
    """
    Extracts all files from a ZIP archive to a target directory.
    
    If the ZIP file is corrupted, a BadZipFile exception is raised. Optionally deletes the original ZIP file after extraction if remove_originals is True.
    
    Args:
        zip_path: Path to the ZIP archive.
        extract_to: Directory where files will be extracted.
        remove_originals: If True, deletes the ZIP file after extraction.
    
    Raises:
        FileNotFoundError: If the ZIP file does not exist.
        zipfile.BadZipFile: If the ZIP file is invalid or corrupted.
        PermissionError: If there are file access permission issues.
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
    logger.info(f"Extracted '{zip_path}' into '{extract_to}'")

    if remove_originals:
        os.remove(zip_path)
        logger.info(f"Removed archive file '{zip_path}' after decompression.")
