"""UniProt PTM drift probe: the release UniProt itself reports, plus the query's result count.

The REST search endpoint is a live JSON API: no ``Last-Modified``, no archive. Probe
version 1 therefore recorded ``source_version: null`` and hashed only the KEYS of the
first result -- the shape of the API response, which changes when UniProt changes its
API and not when it releases new annotations. A UniProt release (roughly every eight
weeks, routinely revising PTM annotations) left that fingerprint byte-identical (#271).

UniProt exposes what the probe needs on every response: ``X-UniProt-Release`` (the
release tag, ``2026_03``), ``X-UniProt-Release-Date``, and ``X-Total-Results`` (the
number of entries matching the exact query the downloader runs). The release is the
version; the count is a content signal in ``extras``; the date is informational, since
it is a function of the release. The key list and the local TSV column order stay as
the schema signal under ``headers``/``checksums``.

Fail-closed: a response without the release or the count is not fingerprinted. The
drift bot regenerates drifted baselines unattended, so one transient omission would
otherwise be committed as the new baseline and retire content detection for good --
the reasoning clingen, hgnc and gencc already follow.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone

import requests

from hvantk.core.plugin.api import DriftProbeError
from hvantk.core.utils.http import request_with_retry
from hvantk.skills.uniprot_ptm.shared.constants import (
    UNIPROT_API_FIELDS,
    UNIPROT_API_URL,
    UNIPROT_HUMAN_PTM_QUERY,
)
from hvantk.skills.uniprot_ptm.shared.datasets import _TSV_COLUMNS, _build_search_url

PROBE_VERSION = 2
_FILENAME = "uniprot-ptm-human.tsv"

# drift_runner wraps every probe in a 60s SIGALRM by default
# (`run_drift_checks(timeout=60)`; the drift workflows pass no `--timeout`), so the
# retry budget below must fit under it with room to spare. Worst case -- every
# attempt exhausts both timeouts and every retry sleeps the full backoff:
# 2+4 s backoff (attempts 1->2, 2->3) + 3x15 s worst-case timeout (5 s connect +
# 10 s read, per attempt) = 51 s < 60 s. The 51 s bound covers the exponential-backoff
# path; a Retry-After sleep is clamped to max_sleep_s and can exceed it, and the drift
# runner's SIGALRM is the backstop in that case.
_ATTEMPTS = 3
_TIMEOUT_S = (5.0, 10.0)

_RELEASE_HEADER = "X-UniProt-Release"
_RELEASE_DATE_HEADER = "X-UniProt-Release-Date"
_TOTAL_HEADER = "X-Total-Results"


def fetch_fingerprint() -> dict:
    """Fingerprint the live UniProt PTM query: release tag, result count, response keys."""
    probe_url = _build_search_url(
        UNIPROT_API_URL, UNIPROT_HUMAN_PTM_QUERY, UNIPROT_API_FIELDS, 1
    )
    try:
        resp = request_with_retry(
            "GET",
            probe_url,
            timeout=_TIMEOUT_S,
            headers={"Accept": "application/json"},
            allow_redirects=True,
            attempts=_ATTEMPTS,
            # UniProt is a JSON API; an empty 200 is the same invisible transient
            # fault #352 documents for medRxiv, not a legitimate empty response.
            retry_on_empty_body=True,
        )
        resp.raise_for_status()
    except requests.RequestException as exc:
        raise DriftProbeError(f"HTTP failure: {exc}") from exc
    try:
        payload = resp.json()
    except ValueError as exc:
        raise DriftProbeError(f"Non-JSON response: {exc}") from exc

    release = (resp.headers.get(_RELEASE_HEADER) or "").strip()
    if not release:
        raise DriftProbeError(
            f"UniProt response carried no {_RELEASE_HEADER}; refusing to record a "
            "fingerprint with no version signal."
        )
    total_raw = resp.headers.get(_TOTAL_HEADER)
    try:
        total_results = int(total_raw)
    except (TypeError, ValueError):
        raise DriftProbeError(
            f"UniProt response carried no usable {_TOTAL_HEADER} (got {total_raw!r}); "
            "refusing to record a fingerprint with no content signal."
        ) from None

    results = payload.get("results", [])
    if not results:
        raise DriftProbeError(
            "UniProt size=1 probe returned zero results; cannot fingerprint."
        )
    if total_results < len(results):
        raise DriftProbeError(
            f"UniProt reported {_TOTAL_HEADER}={total_results}, impossible for "
            f"{len(results)} returned result(s); refusing to record a fingerprint "
            "with an inconsistent content signal."
        )
    entry_keys = sorted(results[0].keys())
    canonical = json.dumps(entry_keys, separators=(",", ":"), sort_keys=True)
    checksum = hashlib.sha256(canonical.encode("utf-8")).hexdigest()

    return {
        "probe_version": PROBE_VERSION,
        "source_version": release,
        "headers": {
            _FILENAME: list(_TSV_COLUMNS),
            "uniprot_entry_keys": entry_keys,
        },
        "checksums": {_FILENAME: checksum},
        "extras": {"total_results": total_results},
        "informational": {"release_date": resp.headers.get(_RELEASE_DATE_HEADER)},
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
