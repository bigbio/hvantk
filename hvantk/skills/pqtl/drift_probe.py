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

from datetime import datetime, timezone

import requests

from hvantk.core.plugin.api import DriftProbeError
from hvantk.core.utils.http import request_with_retry

PROBE_VERSION = 2

PQTL_SOURCE_DOI = "10.1101/2025.01.10.25320181"
MEDRXIV_API_URL = f"https://api.biorxiv.org/details/medrxiv/{PQTL_SOURCE_DOI}"

# Only these fields form the compared surface; see module docstring. `doi` is
# deliberately absent: the request URL is built FROM the DOI, so echoing it back
# is a constant that can never drift.
_COMPARED_FIELDS = ("version", "date", "published")

_FILENAME = "medrxiv-preprint-metadata"
_TIMEOUT_S = (5.0, 15.0)


def fetch_fingerprint() -> dict:
    """Fingerprint the upstream preprint's version metadata."""
    try:
        resp = request_with_retry(
            "GET", MEDRXIV_API_URL, timeout=_TIMEOUT_S, allow_redirects=True
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

    # Selected by version rather than by position. The API's ordering is not
    # documented, and `collection[-1]` on a newest-first response would pin the
    # OLDEST record -- so a v2 posting, the single event this probe exists to
    # detect, would compare equal to the baseline and report clean forever.
    def _version_of(record: dict) -> int:
        raw = record.get("version")
        try:
            return int(raw)
        except (TypeError, ValueError):
            return -1

    latest = max(collection, key=_version_of)
    version = _version_of(latest)
    # Fail closed rather than stringifying a missing field. `str(None)` yields the
    # literal "None", which is truthy and indistinguishable from a real version
    # label, and the bot would commit it as the baseline.
    if version < 0:
        raise DriftProbeError(
            f"medRxiv record for {PQTL_SOURCE_DOI} carried no usable version "
            f"field (got {latest.get('version')!r}); the API shape has probably "
            "changed."
        )

    compared = {field: latest.get(field) for field in _COMPARED_FIELDS}
    # Normalised so an API that switches "1" to 1 does not read as drift.
    compared["version"] = str(version)

    return {
        "probe_version": PROBE_VERSION,
        "source_version": str(version),
        # The metadata IS the signal; a sha256 over it would be a pure function of
        # values already in the compared surface.
        "headers": {_FILENAME: compared},
        "checksums": {},
        "informational": {
            "title": latest.get("title"),
            "versions_listed": len(collection),
        },
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
