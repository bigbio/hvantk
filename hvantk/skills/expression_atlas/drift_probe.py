"""Expression Atlas drift probe: experiment-index JSON fingerprint.

Expression Atlas hosts many per-accession experiments, each with its
own release cadence. There is no single "Expression Atlas version"
header. The EBI gxa portal does expose a JSON index of every
experiment at ``https://www.ebi.ac.uk/gxa/json/experiments`` listing
``experimentAccession`` and ``lastUpdate`` for each. That index
provides a useful (if noisy) global drift signal: any experiment
load/update anywhere in the atlas flips the fingerprint.

The probe fetches the index JSON, projects each experiment down to
``(experimentAccession, lastUpdate)``, sorts the resulting list for
determinism, and hashes the sorted JSON-serialized tuples. The
``source_version`` is the count of experiments, which is the most
human-meaningful single value we can extract.

Note: this is a "global atlas freshness" probe, not a per-experiment
probe. Per-accession drift detection should be handled by the builder
when materializing each experiment (the registry/transcriptomics
catalog already pins per-accession versions for reproducible builds).
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone

import requests

from hvantk.core.plugin_api import DriftProbeError

PROBE_VERSION = 1
EXPRESSION_ATLAS_INDEX_URL = "https://www.ebi.ac.uk/gxa/json/experiments"
_FILENAME = "experiments.json"
_TIMEOUT_S = 60


def fetch_fingerprint() -> dict:
    """Fingerprint of the Expression Atlas experiment index."""
    try:
        resp = requests.get(
            EXPRESSION_ATLAS_INDEX_URL, timeout=_TIMEOUT_S, allow_redirects=True
        )
        resp.raise_for_status()
        payload = resp.json()
    except requests.RequestException as exc:
        raise DriftProbeError(f"HTTP failure: {exc}") from exc
    except ValueError as exc:
        raise DriftProbeError(f"Expression Atlas JSON parse failure: {exc}") from exc

    experiments = payload.get("experiments", []) or []
    projection = sorted(
        (
            (
                exp.get("experimentAccession", ""),
                exp.get("lastUpdate", ""),
            )
            for exp in experiments
            if isinstance(exp, dict)
        )
    )
    canonical = json.dumps(projection, separators=(",", ":")).encode("utf-8")
    checksum = hashlib.sha256(canonical).hexdigest()

    return {
        "probe_version": PROBE_VERSION,
        "source_version": f"{len(projection)} experiments",
        "headers": {
            _FILENAME: ["experimentAccession", "lastUpdate"],
        },
        "checksums": {_FILENAME: checksum},
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }
