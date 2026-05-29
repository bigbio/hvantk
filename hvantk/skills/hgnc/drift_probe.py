"""HGNC drift probe: HEAD + first-line fetch against the complete-set TSV."""

from __future__ import annotations

import hashlib
import io
from datetime import datetime, timezone

import requests

from hvantk.core.plugin.api import DriftProbeError
from hvantk.skills.hgnc.shared.constants import HGNC_DOWNLOAD_URL

PROBE_VERSION = 1
_FILENAME = "hgnc_complete_set.txt"
_TIMEOUT_S = 30


def fetch_fingerprint() -> dict:
    """Lightweight fingerprint of the live HGNC source.

    Reads only the HEAD (for Last-Modified) and the first line (for column
    headers). Does not stream the full ~50 MB TSV.
    """
    try:
        head = requests.head(HGNC_DOWNLOAD_URL, timeout=_TIMEOUT_S, allow_redirects=True)
        head.raise_for_status()
        last_modified = head.headers.get("Last-Modified")
        with requests.get(
            HGNC_DOWNLOAD_URL, timeout=_TIMEOUT_S, stream=True, allow_redirects=True
        ) as resp:
            resp.raise_for_status()
            buf = io.StringIO()
            for chunk in resp.iter_content(chunk_size=4096, decode_unicode=True):
                buf.write(chunk)
                if "\n" in buf.getvalue():
                    break
            first_line = buf.getvalue().splitlines()[0]
    except requests.RequestException as exc:
        raise DriftProbeError(f"HTTP failure: {exc}") from exc

    columns = first_line.split("\t")
    checksum = hashlib.sha256(first_line.encode("utf-8")).hexdigest()
    return {
        "probe_version": PROBE_VERSION,
        "source_version": last_modified,
        "headers": {_FILENAME: columns},
        "checksums": {_FILENAME: checksum},
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
