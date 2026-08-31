"""dbnsfp drift probe: release-list scrape of the upstream landing page.

dbNSFP has no reachable direct data URL, which issue #177 recorded correctly --
but the conclusion that nothing could be probed does not follow. The landing page
itself is fetchable and advertises the release set, so a release roll-over is
detectable even though the archives are not.

Why not the archives: the landing page links every release at
``https://dbnsfp.s3.amazonaws.com/dbNSFP<version>.zip``, and every one of those
is dead -- the bucket answers ``NoSuchBucket``, so the links 404 rather than
merely being gated. The alternative mirror the page names
(``database.liulab.science``) does not resolve at all. Acquisition therefore
stays manual, and this probe deliberately reports on the *advertised release
set* rather than on any data file.

**The raw page must never be hashed.** Google Sites re-renders per request: two
consecutive fetches returned 352,830 and 352,716 bytes. Hashing the body would
flag drift on literally every probe run and open a nightly no-op pull request.
The extracted release list is stable across the same two fetches (31 entries,
identical digest), so the probe fingerprints that projection instead. This is the
msigdb index-page pattern, and here it is a correctness requirement rather than a
stylistic choice.

What this detects: a new dbNSFP release being published, or the download links
being repaired. What it cannot detect: an in-place change to an archive's
contents. That limit is a property of how dbNSFP is published, and is recorded in
SKILL.md so the coverage claim stays honest.
"""

from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime, timezone

import requests

from hvantk.core.plugin.api import DriftProbeError

PROBE_VERSION = 1
DBNSFP_LANDING_URL = "https://sites.google.com/site/jpopgen/dbNSFP"

# Matches the release archives the page advertises, e.g. dbNSFP4.9a.zip,
# dbNSFP4.0b1a.zip, dbNSFP2.0.zip. Anchored on the "dbNSFP" prefix and the
# ".zip" suffix so page furniture cannot leak into the fingerprint.
_RELEASE_REGEX = re.compile(r"dbNSFP(\d[\w.]*?)\.zip", re.IGNORECASE)

# Academic releases carry the "a" suffix ("c" is the commercial build); hvantk
# builds from the academic one, so that is what `source_version` reports.
_ACADEMIC_REGEX = re.compile(r"^(\d+)\.(\d+)a$")

_FILENAME = "dbNSFP-release-index"
_TIMEOUT_S = 30


def _latest_academic(versions: list[str]) -> str | None:
    """Highest ``<major>.<minor>a`` release, compared numerically."""
    parsed = []
    for v in versions:
        m = _ACADEMIC_REGEX.match(v)
        if m:
            parsed.append(((int(m.group(1)), int(m.group(2))), v))
    if not parsed:
        return None
    return max(parsed)[1]


def fetch_fingerprint() -> dict:
    """Fingerprint the release set advertised on the dbNSFP landing page."""
    try:
        resp = requests.get(
            DBNSFP_LANDING_URL, timeout=_TIMEOUT_S, allow_redirects=True
        )
        resp.raise_for_status()
        body = resp.text
    except requests.RequestException as exc:
        raise DriftProbeError(f"HTTP failure: {exc}") from exc

    versions = sorted({m.group(1) for m in _RELEASE_REGEX.finditer(body)})
    # Fail closed. An empty match set means the page moved or its markup changed,
    # not that dbNSFP shipped zero releases; recording it would bake an empty
    # baseline that every later run compares equal to.
    if not versions:
        raise DriftProbeError(
            "dbNSFP landing page advertised no dbNSFP<version>.zip releases; "
            "the page layout has probably changed."
        )

    canonical = json.dumps(versions, separators=(",", ":")).encode("utf-8")
    checksum = hashlib.sha256(canonical).hexdigest()

    return {
        "probe_version": PROBE_VERSION,
        "source_version": _latest_academic(versions),
        "headers": {_FILENAME: ["release_version"]},
        "checksums": {_FILENAME: checksum},
        "extras": {"releases_found": versions},
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
