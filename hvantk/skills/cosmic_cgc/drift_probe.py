"""cosmic-cgc drift probe: release index from the public release notes.

The Cancer Gene Census *data* is login- and licence-gated, which issue #177
recorded correctly: ``cancer.sanger.ac.uk/census`` answers 302 to
``/cosmic/login``, so no direct data URL exists and acquisition stays manual.
The conclusion that nothing could be probed does not follow, though. COSMIC
publishes its release notes without a login, and those carry a per-release
anchor, so a release roll-over is detectable even though the archive is not.

Note the trailing slash matters: ``/cosmic/release_notes`` returns 200 while
``/cosmic/release_notes/`` redirects to the login page.

**Anchor on the id attributes, never on prose.** A first version matched
``COSMIC\\s+v(\\d+)`` anywhere in the body and produced
``[v16, v18, v20, v101, v102, v103, v104]`` -- a non-contiguous set that is not a
release index at all. v16/v18/v20 come from sentences about the *Actionability*
product ("COSMIC v20 of the Actionability data"), a different product line with
its own version series. Any editorial sentence naming an old release would have
entered the compared surface and opened a no-op pull request, and a
forward-looking "coming in COSMIC v105" would have reported a release that did
not exist. The page instead carries ``id="v101"`` ... ``id="v104"`` anchors, one
per real release, which is what this probe reads.

What this detects: a new COSMIC release. What it cannot detect: a change to the
Census contents within a release, which no unauthenticated probe can see. That
limit is a property of the licence gate and is recorded in SKILL.md.
"""

from __future__ import annotations

import re
from datetime import datetime, timezone

import requests

from hvantk.core.plugin.api import DriftProbeError
from hvantk.core.utils.http import request_with_retry

PROBE_VERSION = 2
COSMIC_RELEASE_NOTES_URL = "https://cancer.sanger.ac.uk/cosmic/release_notes"

# Per-release anchors, e.g. id="v104". Deliberately NOT a prose pattern.
_RELEASE_ANCHOR_RE = re.compile(r'id="v(\d+)"', re.IGNORECASE)

_FILENAME = "cosmic-release-index"
_TIMEOUT_S = (5.0, 15.0)


def fetch_fingerprint() -> dict:
    """Fingerprint the release index on the public COSMIC release notes."""
    try:
        resp = request_with_retry(
            "GET", COSMIC_RELEASE_NOTES_URL, timeout=_TIMEOUT_S, allow_redirects=True
        )
        resp.raise_for_status()
    except requests.RequestException as exc:
        raise DriftProbeError(f"HTTP failure: {exc}") from exc

    # The host is known to redirect to a login form (the trailing-slash path does
    # exactly that), and a login page answers 200. Without this a redirect would
    # be scraped as if it were the release notes.
    if resp.history:
        raise DriftProbeError(
            f"COSMIC release notes redirected to {resp.url!r}; the page has "
            "probably moved or now requires a login."
        )

    # Decoded explicitly: requests falls back to ISO-8859-1 for text/html with no
    # charset, which would mangle a non-breaking space and silently change what
    # the pattern matches.
    body = resp.content.decode("utf-8", errors="replace")
    versions = sorted({int(m.group(1)) for m in _RELEASE_ANCHOR_RE.finditer(body)})
    if not versions:
        raise DriftProbeError(
            "COSMIC release notes carried no 'id=\"v<N>\"' release anchors; the "
            "page layout has probably changed."
        )

    return {
        "probe_version": PROBE_VERSION,
        "source_version": f"v{max(versions)}",
        # The release list IS the signal; a sha256 over it would be a pure
        # function of a value already in the compared surface and would add no
        # detection power.
        "headers": {_FILENAME: [f"v{v}" for v in versions]},
        "checksums": {},
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
