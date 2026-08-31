"""INSIDER drift probe: HEAD against both distributed interactome products.

The Interactome Insider portal (http://interactomeinsider.yulab.org) presents its
files through a download page that carries no links in its markup, which is why
issue #177 recorded the source as needing a page-scraping probe and shipped a
stub sentinel instead. The underlying paths are stable and directly addressable,
so no scraping is required:

* ``/bed/all.bed`` -- ``Whole_Human_Interactome_Interface_hg38.bed`` (~1.17 GB),
  the interval-keyed source for ``insider:variants``;
* ``/downloads/interfacesALL/H_sapiens_interfacesALL.txt`` (~49 MB), the
  protein-pair source for ``insider:interfaces``.

Both datasets share this module, so one probe fingerprints both objects and the
two datasets move together. Only HEAD requests are issued; the 1.17 GB BED is
never transferred.

Compared surface: Content-Length per file. ``ETag`` is recorded but demoted to
``informational``, because nginx derives it from (size, mtime) rather than from
the body -- ``0x2f752a5`` in the interfaces ETag is exactly its 49,762,981-byte
length. It therefore adds nothing to Content-Length while inheriting
Last-Modified's failure mode, where a byte-identical re-upload opens a pull
request carrying no information (the hgnc precedent: 8 of 8 drift PRs moved only
a timestamp). ``Last-Modified`` is demoted for the same reason.
"""

from __future__ import annotations

from datetime import datetime, timezone

import requests

from hvantk.core.plugin.api import DriftProbeError

PROBE_VERSION = 1

INSIDER_BASE_URL = "http://interactomeinsider.yulab.org"
INSIDER_BED_URL = f"{INSIDER_BASE_URL}/bed/all.bed"
INSIDER_INTERFACES_URL = (
    f"{INSIDER_BASE_URL}/downloads/interfacesALL/H_sapiens_interfacesALL.txt"
)

# Keyed by the filename the SKILL.md and catalog use, not by the URL basename
# (`all.bed`), so the fingerprint reads against the documented product names.
_FILES = {
    "Whole_Human_Interactome_Interface_hg38.bed": INSIDER_BED_URL,
    "H_sapiens_interfacesALL.txt": INSIDER_INTERFACES_URL,
}

_TIMEOUT_S = 30

# nginx gzips text/plain on the fly when the client offers it, and a compressed
# response carries no Content-Length at all while its ETag gains a "-gzip"
# suffix. requests sends "gzip, deflate" by default, so the probe must opt out
# explicitly or its only compared signal disappears for the interfaces file.
_HEADERS = {"Accept-Encoding": "identity"}


def fetch_fingerprint() -> dict:
    """Lightweight fingerprint of both live INSIDER products (HEAD only)."""
    extras: dict[str, str] = {}
    informational: dict[str, dict[str, str | None]] = {}

    for filename, url in _FILES.items():
        try:
            head = requests.head(
                url, timeout=_TIMEOUT_S, allow_redirects=True, headers=_HEADERS
            )
            head.raise_for_status()
        except requests.RequestException as exc:
            raise DriftProbeError(f"HTTP failure for {filename}: {exc}") from exc

        content_length = head.headers.get("Content-Length")
        # Fail closed: Content-Length is the only compared signal here, so
        # recording a fingerprint without it would compare equal forever.
        if content_length is None:
            raise DriftProbeError(
                f"INSIDER response for {filename} omitted Content-Length; refusing "
                "to record a fingerprint with no content signal."
            )

        etag = head.headers.get("ETag")
        extras[filename] = content_length
        informational[filename] = {
            "etag": etag.strip('"') if etag is not None else None,
            "last_modified": head.headers.get("Last-Modified"),
        }

    return {
        "probe_version": PROBE_VERSION,
        "source_version": None,
        "headers": {},
        # Deliberately empty: this probe fingerprints HTTP validators rather than
        # bodies, the same shape peptideatlas ships. An empty map is not a
        # hand-seeded-baseline marker (an empty *value* inside it would be).
        "checksums": {},
        "extras": extras,
        "informational": informational,
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
