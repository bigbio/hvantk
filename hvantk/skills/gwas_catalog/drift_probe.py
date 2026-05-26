"""GWAS Catalog drift probe: HEAD probe of the EBI FTP latest-release archive.

The NHGRI-EBI GWAS Catalog publishes its full v1.0 associations as a zipped
TSV at a stable URL under the EBI FTP server
(``https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations-full.zip``).
The legacy ``api/search/downloads/full`` endpoint has been retired (404).

The archive is a binary zip; we cannot cheaply stream the first column line
out of it without unpacking. We instead fingerprint the HEAD metadata --
``Last-Modified``, ``ETag``, and ``Content-Length`` -- which together flip
deterministically whenever EBI publishes a new release (typically weekly).
This keeps the probe small, network-light, and honest about what it tracks.
"""

from __future__ import annotations

import hashlib
from datetime import datetime, timezone

import requests

from hvantk.core.plugin.api import DriftProbeError

PROBE_VERSION = 1
GWAS_CATALOG_FULL_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
    "gwas-catalog-associations-full.zip"
)
_FILENAME = "gwas-catalog-associations-full.zip"
_TIMEOUT_S = 30


def fetch_fingerprint() -> dict:
    """Fingerprint of the live GWAS Catalog full-associations archive.

    Issues a single HEAD request and captures ``Last-Modified``, ``ETag``,
    and ``Content-Length``. The sha256 over the ETag+Content-Length pair
    is the file-level checksum; ``source_version`` is the ``Last-Modified``
    timestamp.
    """
    try:
        head = requests.head(
            GWAS_CATALOG_FULL_URL, timeout=_TIMEOUT_S, allow_redirects=True
        )
        head.raise_for_status()
    except requests.RequestException as exc:
        raise DriftProbeError(f"HTTP failure: {exc}") from exc

    last_modified = head.headers.get("Last-Modified")
    etag = head.headers.get("ETag", "")
    content_length = head.headers.get("Content-Length", "")
    fingerprint_blob = f"{etag}|{content_length}".encode("utf-8")
    checksum = hashlib.sha256(fingerprint_blob).hexdigest()
    return {
        "probe_version": PROBE_VERSION,
        "source_version": last_modified,
        "headers": {_FILENAME: ["ETag", "Content-Length", "Last-Modified"]},
        "checksums": {_FILENAME: checksum},
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
