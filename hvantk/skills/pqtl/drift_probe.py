"""pqtl drift probe: preprint-version check against the medRxiv API.

The Fang et al. pQTL summary statistics are published as supplementary material
to a preprint, which issue #177 recorded correctly: there is no direct data URL,
and acquisition stays manual. The conclusion that nothing could be probed does
not follow. For a publication-only source, the publication *is* the upstream, and
medRxiv exposes its metadata through a public JSON API. A new preprint version,
or the preprint being published in a journal, is precisely the event after which
the supplementary data may no longer match what an artifact was built from.

Source: Fang et al. (2025), "Regulation of protein abundance in normal human
tissues", medRxiv, doi:10.1101/2025.01.10.25320181 -- 10,841 proteins across
>700 GTEx samples in five tissues, the cis-pQTL allpairs the builder consumes.

Compared surface: the version number, the posting date, and the ``published``
field (``"NA"`` until a journal version exists, then the journal DOI). The title
and abstract are deliberately excluded -- an editorial typo fix would otherwise
open a pull request carrying no information, the failure mode the hgnc probe was
corrected for.

What this detects: a new preprint version, or journal publication. What it cannot
detect: an in-place replacement of a supplementary file under an unchanged
version. That limit is a property of how the source is published and is recorded
in SKILL.md so the coverage claim stays honest.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone

import requests

from hvantk.core.plugin.api import DriftProbeError

PROBE_VERSION = 1

PQTL_SOURCE_DOI = "10.1101/2025.01.10.25320181"
MEDRXIV_API_URL = f"https://api.biorxiv.org/details/medrxiv/{PQTL_SOURCE_DOI}"

# Only these fields form the compared surface; see module docstring.
_COMPARED_FIELDS = ("version", "date", "published", "doi")

_FILENAME = "medrxiv-preprint-metadata"
_TIMEOUT_S = 30


def fetch_fingerprint() -> dict:
    """Fingerprint the upstream preprint's version metadata."""
    try:
        resp = requests.get(MEDRXIV_API_URL, timeout=_TIMEOUT_S, allow_redirects=True)
        resp.raise_for_status()
    except requests.RequestException as exc:
        raise DriftProbeError(f"HTTP failure: {exc}") from exc

    # Parsed in its own block: requests' JSONDecodeError subclasses both
    # ValueError and RequestException, so decoding inside the block above would
    # report a malformed body as an HTTP failure.
    try:
        payload = resp.json()
    except ValueError as exc:
        raise DriftProbeError(f"medRxiv API returned non-JSON: {exc}") from exc

    collection = payload.get("collection") or []
    # Fail closed. An empty collection means the DOI stopped resolving or the API
    # changed shape, not that the preprint has no versions; recording it would
    # bake an empty baseline that every later run compares equal to.
    if not collection:
        raise DriftProbeError(
            f"medRxiv API returned no records for {PQTL_SOURCE_DOI}; the DOI or "
            "the API shape has probably changed."
        )

    # The API lists one record per version, oldest first; the newest is what a
    # rebuild would pick up.
    latest = collection[-1]
    compared = {field: latest.get(field) for field in _COMPARED_FIELDS}

    canonical = json.dumps(compared, sort_keys=True, separators=(",", ":")).encode(
        "utf-8"
    )
    checksum = hashlib.sha256(canonical).hexdigest()

    return {
        "probe_version": PROBE_VERSION,
        "source_version": str(latest.get("version")),
        "headers": {_FILENAME: list(_COMPARED_FIELDS)},
        "checksums": {_FILENAME: checksum},
        "extras": compared,
        "informational": {
            "title": latest.get("title"),
            "versions_listed": len(collection),
        },
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
