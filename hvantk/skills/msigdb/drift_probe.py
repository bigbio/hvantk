"""MSigDB drift probe: index-page version scrape.

MSigDB distributes per-release GMT files at versioned URLs (e.g.
``c2.cp.v2026.1.Hs.symbols.gmt``); there is no single "latest" stable
URL. Releases are documented on the GSEA-MSigDB portal index page at
``https://www.gsea-msigdb.org/gsea/msigdb/index.jsp`` which surfaces
the current human + mouse version strings (``v2026.1.Hs`` /
``v2026.1.Mm``).

The probe fetches the index page, regex-matches all ``vYYYY.N.[HM]s``
tokens, and fingerprints the sorted list. ``source_version`` is the
highest version found (lexicographic over the canonical zero-padded
form). Any release roll-over (annual cadence) flips the checksum.
"""

from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime, timezone

import requests

from hvantk.core.plugin.api import DriftProbeError

PROBE_VERSION = 1
MSIGDB_INDEX_URL = "https://www.gsea-msigdb.org/gsea/msigdb/index.jsp"
_VERSION_REGEX = re.compile(r"v(20\d{2}\.\d+(?:\.[A-Za-z]+)?)")
_FILENAME = "index.jsp"
_TIMEOUT_S = 30


def fetch_fingerprint() -> dict:
    """Scrape the MSigDB index page for current release version strings."""
    try:
        resp = requests.get(MSIGDB_INDEX_URL, timeout=_TIMEOUT_S, allow_redirects=True)
        resp.raise_for_status()
        body = resp.text
    except requests.RequestException as exc:
        raise DriftProbeError(f"HTTP failure: {exc}") from exc

    versions = sorted(set(_VERSION_REGEX.findall(body)))
    if not versions:
        raise DriftProbeError(
            "MSigDB index page did not contain any vYYYY.N.[HM]s version tokens"
        )

    canonical = json.dumps(versions, separators=(",", ":")).encode("utf-8")
    checksum = hashlib.sha256(canonical).hexdigest()

    return {
        "probe_version": PROBE_VERSION,
        "source_version": versions[-1],
        "headers": {_FILENAME: ["release_version"]},
        "checksums": {_FILENAME: checksum},
        "fetched_at": datetime.now(timezone.utc).isoformat(),
        "extras": {"versions_found": versions},
    }
