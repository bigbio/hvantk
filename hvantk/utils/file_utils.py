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

    response = requests.get(url, stream=True)
    if response.status_code == 200:
        total_size = int(response.headers.get('content-length', 0))
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
        print(f"Download completed: {local_path}")
    else:
        raise Exception(f"Failed to download {url}: Status code {response.status_code}")
