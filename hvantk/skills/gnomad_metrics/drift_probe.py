"""Drift probe for gnomad-metrics — documentation-only source (stub).

gnomAD constraint metrics are acquired manually here and ship no committed
baseline fingerprint. A HEAD probe of the public v2.1.1 constraint file on GCS is
marginally feasible but deferred (see issue #177, Option 2). For now this probe
returns a structured stub sentinel so ``hvantk drift`` reports a visible WARNING
(status="stub") rather than a silent false-green.
"""
from __future__ import annotations

from hvantk.core.plugin.api import stub_fingerprint

_REASON = (
    "gnomAD constraint metrics acquired manually; no committed baseline "
    "(a GCS HEAD probe of the public v2.1.1 file is feasible but deferred — "
    "issue #177)"
)


def fetch_fingerprint() -> dict:
    return stub_fingerprint(_REASON)
