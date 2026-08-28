"""GenCC drift probe: HEAD + first-line fetch against the submissions TSV.

GenCC serves a single rolling submissions TSV at a stable URL
(``GENCC_BASE_URL``); there are no dated archives. The body is a
tab-separated table whose first line is the column-header row (begins
with ``sgc_id``). We only need a fingerprint sensitive enough to surface
upstream column drift, so we stream the response until we have read the
column-header line, hash it, and stop.

``Last-Modified`` is recorded under ``informational`` rather than
``source_version``: GenCC re-publishes on a weekly cadence without a schema
or content change, and across all 4 fingerprint commits from 2026-07-27 to
2026-08-23 the checksum ``6f07ac79f9...`` never moved while Last-Modified
moved every time, so treating it as the source version opened a
no-information PR on every run. ``Content-Length`` is the real content
signal, matching the ClinGen and HGNC probes' reasoning.
"""

from __future__ import annotations

import hashlib
import io
from datetime import datetime, timezone

import requests

from hvantk.skills.gencc.shared.constants import GENCC_BASE_URL, GENCC_FILE_PREFIX
from hvantk.core.plugin.api import DriftProbeError
from hvantk.core.utils.http import request_with_retry

PROBE_VERSION = 2
_FILENAME = f"{GENCC_FILE_PREFIX}.tsv"
_TIMEOUT_S = 30
_HEADER_MARKER = "sgc_id"


def fetch_fingerprint() -> dict:
    """Lightweight fingerprint of the live GenCC submissions TSV.

    Reads only the HTTP HEAD (for ``Content-Length`` and ``Last-Modified``,
    when present) and the portion of the body up to and including the
    column-header line (begins with ``sgc_id``). Does not stream the full
    TSV (~MB-scale).
    """
    # Retrying, because GenCC rate-limits: on 2026-08-04 it answered the scheduled
    # drift regeneration with HTTP 429, which failed the probe and then the whole
    # workflow run. request_with_retry honours Retry-After but clamps the wait, so a
    # long backoff request cannot stall CI.
    try:
        head = request_with_retry(
            "HEAD", GENCC_BASE_URL, timeout=_TIMEOUT_S, allow_redirects=True
        )
        head.raise_for_status()
        last_modified = head.headers.get("Last-Modified")
        content_length = head.headers.get("Content-Length")
        if content_length is None:
            # Fail closed -- see the hgnc probe for the full reasoning. The bot
            # regenerates drifted baselines automatically, so a partial fingerprint
            # would become permanent.
            raise DriftProbeError(
                "GenCC response omitted Content-Length; refusing to record a "
                "fingerprint with no content signal."
            )

        with request_with_retry(
            "GET", GENCC_BASE_URL, timeout=_TIMEOUT_S, stream=True, allow_redirects=True
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
        "source_version": None,
        "headers": {_FILENAME: columns},
        "checksums": {_FILENAME: checksum},
        "extras": {"content_length": content_length},
        "informational": {"last_modified": last_modified},
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
