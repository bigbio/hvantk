"""GenCC drift probe: HEAD + first-line fetch against the submissions TSV.

GenCC serves a single rolling submissions TSV at a stable URL
(``GENCC_BASE_URL``); there are no dated archives. The body is a
tab-separated table whose first line is the column-header row (begins
with ``sgc_id``). We only need a fingerprint sensitive enough to surface
upstream column drift, so we stream the response until we have read the
column-header line, hash it, and stop. ``Last-Modified`` (if the server
returns one) is captured as the source version.
"""

from __future__ import annotations

import hashlib
import io
from datetime import datetime, timezone

import requests

from hvantk.core.constants import GENCC_BASE_URL, GENCC_FILE_PREFIX
from hvantk.core.plugin_api import DriftProbeError

PROBE_VERSION = 1
_FILENAME = f"{GENCC_FILE_PREFIX}.tsv"
_TIMEOUT_S = 30
_HEADER_MARKER = "sgc_id"


def fetch_fingerprint() -> dict:
    """Lightweight fingerprint of the live GenCC submissions TSV.

    Reads only the HTTP HEAD (for ``Last-Modified``, when present) and the
    portion of the body up to and including the column-header line (begins
    with ``sgc_id``). Does not stream the full TSV (~MB-scale).
    """
    try:
        head = requests.head(
            GENCC_BASE_URL, timeout=_TIMEOUT_S, allow_redirects=True
        )
        head.raise_for_status()
        last_modified = head.headers.get("Last-Modified")

        with requests.get(
            GENCC_BASE_URL, timeout=_TIMEOUT_S, stream=True, allow_redirects=True
        ) as resp:
            resp.raise_for_status()
            buf = io.StringIO()
            header_line = None
            for chunk in resp.iter_content(chunk_size=4096, decode_unicode=True):
                if not chunk:
                    continue
                buf.write(chunk)
                if _HEADER_MARKER in buf.getvalue():
                    # Scan accumulated buffer line-by-line for the header row.
                    for line in buf.getvalue().splitlines():
                        if _HEADER_MARKER in line:
                            header_line = line
                            break
                    if header_line is not None:
                        break
            if header_line is None:
                raise DriftProbeError(
                    f"GenCC header marker {_HEADER_MARKER!r} not found in "
                    f"response body."
                )
    except requests.RequestException as exc:
        raise DriftProbeError(f"HTTP failure: {exc}") from exc

    # Header is TSV; strip surrounding quotes/whitespace for a stable list.
    columns = [col.strip().strip('"') for col in header_line.split("\t")]
    checksum = hashlib.sha256(header_line.encode("utf-8")).hexdigest()
    return {
        "probe_version": PROBE_VERSION,
        "source_version": last_modified,
        "headers": {_FILENAME: columns},
        "checksums": {_FILENAME: checksum},
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
