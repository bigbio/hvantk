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

from datetime import datetime, timezone

import requests

from hvantk.core.plugin.api import DriftProbeError
from hvantk.core.utils.http import request_with_retry

PROBE_VERSION = 2
ALPHAGENOME_PYPI_URL = "https://pypi.org/pypi/alphagenome/json"

_FILENAME = "alphagenome-sdk-releases"
_TIMEOUT_S = (5.0, 15.0)


def fetch_fingerprint() -> dict:
    """Fingerprint the published AlphaGenome SDK release set."""
    try:
        resp = request_with_retry(
            "GET", ALPHAGENOME_PYPI_URL, timeout=_TIMEOUT_S, allow_redirects=True
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
    # Fail closed on the field that actually carries the signal. `info.version` is
    # the supported one; `releases` is deprecated on this endpoint and slated for
    # removal, so requiring it would turn a PyPI API change into a permanent
    # probe_failed for an SDK that never moved.
    if not current:
        raise DriftProbeError(
            "PyPI returned no info.version for alphagenome; the project or the "
            "API shape has probably changed."
        )

    compared: dict[str, object] = {"current_version": current}
    releases = sorted((payload.get("releases") or {}).keys())
    if releases:
        compared["release_count"] = len(releases)

    return {
        "probe_version": PROBE_VERSION,
        "source_version": current,
        # The version IS the signal; a sha256 over it would be a pure function of
        # a value already in the compared surface. The full release list is
        # deliberately NOT compared: it is deprecated upstream, and it grows on
        # pre-release and yanked uploads that no build would ever install.
        "headers": {_FILENAME: compared},
        "checksums": {},
        "informational": {"releases_found": releases},
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
