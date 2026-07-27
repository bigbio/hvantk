"""Drift probe for ensembl-gene:structure.

Unlike the BioMart ``genes`` dataset -- acquired manually, hence a stub probe -- the GTF
has a stable, predictable release URL, so it can be probed for real: a HEAD against
``ENSEMBL_GTF_URL`` records the validators Ensembl's archive serves (``ETag``,
``Last-Modified``, ``Content-Length``) alongside the pinned release.

What this detects is the **pinned** artifact moving underneath a build -- Ensembl
re-issuing a release file, or a contributor bumping ``ENSEMBL_RELEASE`` without
rebuilding. Both change CDS lengths, MANE Select assignments and coordinates for the same
gene, which is a silent correctness bug rather than a loud one.

What it deliberately does not do is report the mere existence of a newer Ensembl release
as drift. The release is a repo-side pin, not an upstream fact: treating release 114's
publication as movement in our source would leave this dataset permanently drifted
through no change to anything it actually reads. Choosing when to advance the pin is an
upgrade decision, not a drift signal.

Before this probe did any network I/O it returned the two repo constants below and
nothing else, so it compared equal to its own committed baseline forever -- a permanent
false "clean" that the module docstring described as a real fingerprint.
"""
from __future__ import annotations

from datetime import datetime, timezone

import requests

from hvantk.core.plugin.api import DriftProbeError
from hvantk.resources.ensembl_release import (
    ENSEMBL_GTF_FILENAME,
    ENSEMBL_GTF_URL,
    ENSEMBL_RELEASE,
)

PROBE_VERSION = 2
_TIMEOUT_S = 30


def fetch_fingerprint() -> dict:
    """HEAD the pinned Ensembl GTF and fingerprint its HTTP validators.

    Downloads no payload: the GTF is ~64 MB, and its validators are enough to
    tell whether the file behind the pin has changed.
    """
    try:
        resp = requests.head(ENSEMBL_GTF_URL, timeout=_TIMEOUT_S, allow_redirects=True)
        resp.raise_for_status()
    except requests.RequestException as exc:
        raise DriftProbeError(f"HTTP failure: {exc}") from exc

    etag = resp.headers.get("ETag")
    content_length = resp.headers.get("Content-Length")
    if etag is None and content_length is None:
        raise DriftProbeError(
            f"{ENSEMBL_GTF_URL} returned neither ETag nor Content-Length; "
            "nothing stable to fingerprint."
        )

    return {
        "probe_version": PROBE_VERSION,
        "source_version": resp.headers.get("Last-Modified"),
        "headers": {
            ENSEMBL_GTF_FILENAME: {
                "etag": etag,
                "content_length": content_length,
            }
        },
        "checksums": {},
        "extras": {"release": ENSEMBL_RELEASE, "url": ENSEMBL_GTF_URL},
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
