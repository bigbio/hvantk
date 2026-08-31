"""cosmic-cgc drift probe: release-version scrape of the public release notes.

The Cancer Gene Census *data* is login- and licence-gated, which issue #177
recorded correctly: ``cancer.sanger.ac.uk/census`` answers 302 to
``/cosmic/login``, so no direct data URL exists and acquisition stays manual.
The conclusion that nothing could be probed does not follow, though. COSMIC
publishes its release notes without a login, and those name the current release,
so a release roll-over is detectable even though the archive is not.

Note the trailing slash matters: ``/cosmic/release_notes`` returns 200 while
``/cosmic/release_notes/`` redirects to the login page.

Compared surface: the sorted set of ``COSMIC v<N>`` tokens the page names, not
the page body. The body is in fact byte-stable across fetches (unlike dbNSFP's
Google Sites page, which re-renders per request), but fingerprinting the
projection rather than the markup keeps an unrelated editorial edit from opening
a no-op pull request.

What this detects: a new COSMIC release. What it cannot detect: a change to the
Census contents within a release, which no unauthenticated probe can see. That
limit is a property of the licence gate and is recorded in SKILL.md so the
coverage claim stays honest.
"""

from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime, timezone

import requests

from hvantk.core.plugin.api import DriftProbeError

PROBE_VERSION = 1
COSMIC_RELEASE_NOTES_URL = "https://cancer.sanger.ac.uk/cosmic/release_notes"

_RELEASE_REGEX = re.compile(r"COSMIC\s+v(\d+)", re.IGNORECASE)
_FILENAME = "cosmic-release-index"
_TIMEOUT_S = 30


def fetch_fingerprint() -> dict:
    """Fingerprint the release set named on the public COSMIC release notes."""
    try:
        resp = requests.get(
            COSMIC_RELEASE_NOTES_URL, timeout=_TIMEOUT_S, allow_redirects=True
        )
        resp.raise_for_status()
        body = resp.text
    except requests.RequestException as exc:
        raise DriftProbeError(f"HTTP failure: {exc}") from exc

    versions = sorted(
        {int(m.group(1)) for m in _RELEASE_REGEX.finditer(body)}
    )
    # Fail closed. No matches means the page moved or was replaced by the login
    # form, not that COSMIC has no releases; recording it would bake an empty
    # baseline that every later run compares equal to.
    if not versions:
        raise DriftProbeError(
            "COSMIC release notes named no 'COSMIC v<N>' releases; the page has "
            "probably moved or redirected to the login form."
        )

    canonical = json.dumps(versions, separators=(",", ":")).encode("utf-8")
    checksum = hashlib.sha256(canonical).hexdigest()

    return {
        "probe_version": PROBE_VERSION,
        "source_version": f"v{max(versions)}",
        "headers": {_FILENAME: ["release_version"]},
        "checksums": {_FILENAME: checksum},
        "extras": {"releases_found": [f"v{v}" for v in versions]},
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
