"""alphagenome drift probe: SDK release check against the PyPI JSON API.

AlphaGenome is a credentialed live prediction service, so there is no static
artifact to fingerprint and issue #177 shipped a stub sentinel. For a live model
API, though, the meaningful upstream change is not a file but a *model or client
release*: predictions are generated on demand, so what makes a stored artifact
stale is the service behind it moving. The SDK is published openly on PyPI, whose
JSON API needs no credentials, so that release stream is directly probeable.

Compared surface: the current version plus the sorted set of released versions.
Upload timestamps and file digests are excluded -- PyPI can re-host an unchanged
release, and the hgnc precedent is that a validator which moves without the
content changing produces nothing but no-op pull requests.

What this detects: a new AlphaGenome SDK release, which is the signal to re-check
whether predictions still match a stored artifact. What it cannot detect: a
server-side model update shipped without an SDK release, which no unauthenticated
probe can see. That limit is a property of the service and is recorded in
SKILL.md so the coverage claim stays honest.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone

import requests

from hvantk.core.plugin.api import DriftProbeError

PROBE_VERSION = 1
ALPHAGENOME_PYPI_URL = "https://pypi.org/pypi/alphagenome/json"

_FILENAME = "alphagenome-sdk-releases"
_TIMEOUT_S = 30


def fetch_fingerprint() -> dict:
    """Fingerprint the published AlphaGenome SDK release set."""
    try:
        resp = requests.get(
            ALPHAGENOME_PYPI_URL, timeout=_TIMEOUT_S, allow_redirects=True
        )
        resp.raise_for_status()
    except requests.RequestException as exc:
        raise DriftProbeError(f"HTTP failure: {exc}") from exc

    # Parsed in its own block: requests' JSONDecodeError subclasses both
    # ValueError and RequestException, so decoding inside the block above would
    # report a malformed body as an HTTP failure.
    try:
        payload = resp.json()
    except ValueError as exc:
        raise DriftProbeError(f"PyPI returned non-JSON: {exc}") from exc

    current = (payload.get("info") or {}).get("version")
    releases = sorted((payload.get("releases") or {}).keys())
    # Fail closed. No version means the project was renamed or the API changed
    # shape, not that AlphaGenome has no releases; recording it would bake an
    # empty baseline that every later run compares equal to.
    if not current or not releases:
        raise DriftProbeError(
            "PyPI returned no version or no releases for alphagenome; the "
            "project or the API shape has probably changed."
        )

    canonical = json.dumps(
        {"current": current, "releases": releases}, separators=(",", ":")
    ).encode("utf-8")
    checksum = hashlib.sha256(canonical).hexdigest()

    return {
        "probe_version": PROBE_VERSION,
        "source_version": current,
        "headers": {_FILENAME: ["sdk_version"]},
        "checksums": {_FILENAME: checksum},
        "extras": {"current_version": current, "releases_found": releases},
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
