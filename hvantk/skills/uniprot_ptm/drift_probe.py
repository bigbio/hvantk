"""UniProt PTM drift probe: HEAD + minimal page-of-1 GET against the REST API.

UniProt's REST search endpoint (``UNIPROT_API_URL``) is a live JSON API,
not a static TSV. There is no ``Last-Modified`` header on a query URL and
no archival versioning. To surface upstream schema drift we issue a HEAD
(reachability) and a small ``size=1`` GET, then fingerprint the keys
present in the first result entry. The keys come directly from the
UniProt JSON schema; if UniProt removes or renames a key the fingerprint
changes.

A second list of column names — the local TSV column order written by
:mod:`hvantk.skills.uniprot_ptm.shared.datasets` — is also captured so
that downstream consumers can compare against the output file shape.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone

import requests

from hvantk.core.plugin.api import DriftProbeError
from hvantk.skills.uniprot_ptm.shared.constants import (
    UNIPROT_API_FIELDS,
    UNIPROT_API_URL,
    UNIPROT_HUMAN_PTM_QUERY,
)
from hvantk.skills.uniprot_ptm.shared.datasets import _TSV_COLUMNS, _build_search_url

PROBE_VERSION = 1
_FILENAME = "uniprot-ptm-human.tsv"
_TIMEOUT_S = 30


def fetch_fingerprint() -> dict:
    """Lightweight fingerprint of the live UniProt PTM REST endpoint.

    Performs a HEAD for reachability and a ``size=1`` GET to capture the
    JSON keys of a single result entry. Hashes the canonicalised key list
    so column drift in the upstream response surfaces as a checksum
    change.
    """
    probe_url = _build_search_url(
        UNIPROT_API_URL, UNIPROT_HUMAN_PTM_QUERY, UNIPROT_API_FIELDS, 1
    )

    try:
        head = requests.head(probe_url, timeout=_TIMEOUT_S, allow_redirects=True)
        head.raise_for_status()
        last_modified = head.headers.get("Last-Modified")

        resp = requests.get(
            probe_url,
            timeout=_TIMEOUT_S,
            headers={"Accept": "application/json"},
            allow_redirects=True,
        )
        resp.raise_for_status()
        payload = resp.json()
    except requests.RequestException as exc:
        raise DriftProbeError(f"HTTP failure: {exc}") from exc
    except ValueError as exc:
        raise DriftProbeError(f"Non-JSON response: {exc}") from exc

    results = payload.get("results", [])
    if not results:
        raise DriftProbeError(
            "UniProt size=1 probe returned zero results; cannot fingerprint."
        )

    entry_keys = sorted(results[0].keys())
    canonical = json.dumps(entry_keys, separators=(",", ":"), sort_keys=True)
    checksum = hashlib.sha256(canonical.encode("utf-8")).hexdigest()

    return {
        "probe_version": PROBE_VERSION,
        "source_version": last_modified,
        "headers": {
            _FILENAME: list(_TSV_COLUMNS),
            "uniprot_entry_keys": entry_keys,
        },
        "checksums": {_FILENAME: checksum},
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
