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
stays manual, and this probe reports on the *advertised release set* rather than
on any data file. That the documented download path is broken is tracked
separately as issue #321; it is a documentation defect, not a drift signal.

**The raw page must never be hashed.** Google Sites re-renders per request: two
consecutive fetches returned 352,830 and 352,716 bytes. Hashing the body would
flag drift on every probe run and open a nightly no-op pull request. The
extracted release list was identical across those same two fetches, so the probe
compares that projection instead. This is the msigdb index-page pattern, and here
it is a correctness requirement rather than a stylistic choice.

What this detects: a new dbNSFP release being advertised. What it cannot detect:
an in-place change to an archive's contents, or the download links being
repaired -- the page markup names the same archives either way, so a restored
bucket would not move this fingerprint. Those limits are recorded in SKILL.md.
"""

from __future__ import annotations

import re
from datetime import datetime, timezone

import requests

from hvantk.core.plugin.api import DriftProbeError
from hvantk.core.utils.http import request_with_retry

PROBE_VERSION = 2
DBNSFP_LANDING_URL = "https://sites.google.com/site/jpopgen/dbNSFP"

# Matches the release archives the page advertises. The optional `v` is
# load-bearing: dbNSFP's 2.x and 3.x generations are published as `dbNSFPv3.5a.zip`,
# and a pattern requiring a digit straight after "dbNSFP" silently matched none of
# them -- so a next release named `dbNSFPv5.0a.zip` would have left the projection
# unchanged and reported clean on a real roll-over.
_RELEASE_REGEX = re.compile(r"dbNSFPv?(\d[\w.]*?)\.zip", re.IGNORECASE)

# Academic releases carry the "a" suffix ("c" is the commercial build); hvantk
# builds from the academic one. Case-insensitive to match _RELEASE_REGEX, which
# would otherwise capture `4.9A` into the set while this rejected it, leaving a
# phantom release in the list and a regressed source_version.
_ACADEMIC_REGEX = re.compile(r"^(\d+(?:\.\d+)*)a$", re.IGNORECASE)

_FILENAME = "dbNSFP-release-index"
_TIMEOUT_S = (5.0, 15.0)


def _latest_academic(versions: list[str]) -> str | None:
    """Highest academic release, compared componentwise rather than as text."""
    parsed = []
    for v in versions:
        m = _ACADEMIC_REGEX.match(v)
        if m:
            parts = tuple(int(p) for p in m.group(1).split("."))
            parsed.append((parts, v))
    if not parsed:
        return None
    # Componentwise so 4.10a beats 4.9a; a text sort would not.
    return max(parsed)[1]


def fetch_fingerprint() -> dict:
    """Fingerprint the release set advertised on the dbNSFP landing page."""
    try:
        resp = request_with_retry(
            "GET", DBNSFP_LANDING_URL, timeout=_TIMEOUT_S, allow_redirects=True
        )
        resp.raise_for_status()
    except requests.RequestException as exc:
        raise DriftProbeError(f"HTTP failure: {exc}") from exc

    # Decoded explicitly: requests falls back to ISO-8859-1 for text/html with no
    # charset, and its chardet/charset_normalizer fallback is environment
    # dependent, so the same page could otherwise yield two different projections
    # on two machines.
    body = resp.content.decode("utf-8", errors="replace")

    # Lowercased before entering the set so a cosmetic recasing of one link
    # cannot present as an extra release.
    versions = sorted({m.group(1).lower() for m in _RELEASE_REGEX.finditer(body)})
    # Fail closed. An empty match set means the page moved or its markup changed,
    # not that dbNSFP shipped zero releases; recording it would bake an empty
    # baseline that every later run compares equal to.
    if not versions:
        raise DriftProbeError(
            "dbNSFP landing page advertised no dbNSFP<version>.zip releases; "
            "the page layout has probably changed."
        )

    latest = _latest_academic(versions)
    # Also fail closed here. Silently recording source_version: null would let the
    # bot commit that null as the baseline, after which the probe reports clean
    # forever having quietly stopped identifying a release at all.
    if latest is None:
        raise DriftProbeError(
            f"dbNSFP advertised {len(versions)} releases but none matched the "
            "academic '<version>a' naming; the release scheme has probably "
            "changed."
        )

    return {
        "probe_version": PROBE_VERSION,
        "source_version": latest,
        # The release list IS the signal; a sha256 over it would be a pure
        # function of a value already in the compared surface.
        "headers": {_FILENAME: versions},
        "checksums": {},
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
