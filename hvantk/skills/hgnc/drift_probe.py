"""HGNC drift probe: HEAD + first-line fetch against the complete-set TSV."""

from __future__ import annotations

import hashlib
import io
from datetime import datetime, timezone

import requests

from hvantk.core.plugin.api import DriftProbeError
from hvantk.skills.hgnc.shared.constants import HGNC_DOWNLOAD_URL

PROBE_VERSION = 2
_FILENAME = "hgnc_complete_set.txt"
_TIMEOUT_S = 30


def fetch_fingerprint() -> dict:
    """Lightweight fingerprint of the live HGNC source.

    Reads only the HEAD (for Content-Length and Last-Modified) and the first line
    (for column headers). Does not stream the full ~50 MB TSV.

    ``Last-Modified`` is recorded under ``informational`` rather than
    ``source_version`` because HGNC re-publishes byte-identical content with a fresh
    timestamp: across all 8 fingerprint commits from 2026-05-16 to 2026-08-27 the
    content checksum never moved while Last-Modified moved every time, and each of
    those opened a PR carrying no information. ``Content-Length`` is the real content
    signal, matching the ClinGen probe's reasoning.
    """
    try:
        head = requests.head(HGNC_DOWNLOAD_URL, timeout=_TIMEOUT_S, allow_redirects=True)
        head.raise_for_status()
        last_modified = head.headers.get("Last-Modified")
        content_length = head.headers.get("Content-Length")
        if content_length is None:
            # Fail closed. Without this the column-header hash is the only signal left,
            # so a row-level HGNC change would read as clean -- and because the bot
            # regenerates drifted baselines automatically, one transient omission would
            # bake that in permanently. probe_failed is loud and recoverable.
            raise DriftProbeError(
                "HGNC response omitted Content-Length; refusing to record a "
                "fingerprint with no content signal."
            )
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
        "source_version": None,
        "headers": {_FILENAME: columns},
        "checksums": {_FILENAME: checksum},
        "extras": {"content_length": content_length},
        "informational": {"last_modified": last_modified},
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
