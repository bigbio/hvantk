"""GTEx eQTL drift probe: portal HEAD + version-meta scrape.

GTEx ships per-tissue v11 cis-eQTL parquet files but does not expose a
machine-readable manifest of the release artifacts. The GTEx portal
download page (``https://gtexportal.org/home/downloads/adult-gtex/qtl``)
does, however, return ``Last-Modified``/``ETag`` headers and the HTML
body contains ``<meta name="gtex-version" content="...">`` and
``<meta name="gtex-build-date" content="...">`` tags that flip when the
portal redeploys (which historically tracks data releases).

The probe fingerprints both the HEAD metadata and the version/build-date
meta-tag pair. ``source_version`` is the gtex-version string from the
HTML; the checksum hashes ``(gtex-version, gtex-build-date, ETag)``.

This is the best signal we have without an authenticated portal API or
gs:// directory listing. It is honest about its limitation: it tracks
"portal freshness", which is highly correlated with -- but not strictly
equivalent to -- v11 dataset releases.
"""

from __future__ import annotations

import hashlib
import re
from datetime import datetime, timezone

import requests

from hvantk.core.plugin_api import DriftProbeError

PROBE_VERSION = 1
GTEX_PORTAL_DOWNLOAD_URL = "https://gtexportal.org/home/downloads/adult-gtex/qtl"
_GTEX_VERSION_META_RE = re.compile(
    r'<meta\s+name="gtex-version"\s+content="([^"]+)"', re.IGNORECASE
)
_GTEX_BUILD_DATE_META_RE = re.compile(
    r'<meta\s+name="gtex-build-date"\s+content="([^"]+)"', re.IGNORECASE
)
_FILENAME = "gtex-portal-qtl-downloads"
_TIMEOUT_S = 30


def fetch_fingerprint() -> dict:
    """Fingerprint the GTEx QTL download page + version meta-tags."""
    try:
        resp = requests.get(
            GTEX_PORTAL_DOWNLOAD_URL, timeout=_TIMEOUT_S, allow_redirects=True
        )
        resp.raise_for_status()
        body = resp.text
        etag = resp.headers.get("ETag", "")
        last_modified = resp.headers.get("Last-Modified", "")
    except requests.RequestException as exc:
        raise DriftProbeError(f"HTTP failure: {exc}") from exc

    version_match = _GTEX_VERSION_META_RE.search(body)
    build_date_match = _GTEX_BUILD_DATE_META_RE.search(body)
    if version_match is None:
        raise DriftProbeError(
            "GTEx portal page did not contain a gtex-version meta tag"
        )

    gtex_version = (version_match.group(1) or "").strip()
    gtex_build_date = (
        (build_date_match.group(1) if build_date_match else "") or ""
    ).strip()

    fingerprint_blob = (
        f"{gtex_version}|{gtex_build_date}|{etag}".encode("utf-8")
    )
    checksum = hashlib.sha256(fingerprint_blob).hexdigest()

    return {
        "probe_version": PROBE_VERSION,
        "source_version": gtex_version,
        "headers": {
            _FILENAME: ["gtex-version", "gtex-build-date", "ETag", "Last-Modified"],
        },
        "checksums": {_FILENAME: checksum},
        "fetched_at": datetime.now(timezone.utc).isoformat(),
        "extras": {
            "gtex_version": gtex_version,
            "gtex_build_date": gtex_build_date,
            "etag": etag,
            "last_modified": last_modified,
        },
    }
