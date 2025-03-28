import os
import requests
from tqdm import tqdm


def download_file(url: str, out_dir: str, file_name: str):
    """
    Download a file from a URL to a local directory.

    :param url: URL of the file to download
    :param out_dir: Local directory to save the file
    :param file_name: Name of the file to save

    :return: None
    """
    os.makedirs(out_dir, exist_ok=True)
    local_path = os.path.join(out_dir, file_name)
    print(f"Downloading {url} to {local_path}")

    try:
        response = requests.get(url, stream=True, timeout=30)  # Add timeout
        response.raise_for_status()  # Raises exception for 4XX/5XX responses
        
        total_size = int(response.headers.get('content-length', 0))
        
        # Optional: Add file size check
        max_size = 1024 * 1024 * 1024  # 1GB
        if total_size > max_size:
            logger.warning(f"File is very large ({total_size/1024/1024:.1f} MB), exceeding recommended size of {max_size/1024/1024:.1f} MB")
        
        with open(local_path, 'wb') as f, tqdm(
                desc=file_name,
                total=total_size,
                unit='B',
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
