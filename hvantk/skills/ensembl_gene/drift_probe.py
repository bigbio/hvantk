"""Drift probe for ensembl-gene — documentation-only source (stub).

Ensembl gene annotations are acquired manually here and ship no committed
baseline fingerprint. A release-endpoint version probe (Ensembl REST) is
marginally feasible but deferred (see issue #177, Option 2). For now this probe
returns a structured stub sentinel so ``hvantk drift`` reports a visible WARNING
(status="stub") rather than a silent false-green.
"""
from __future__ import annotations

from hvantk.core.plugin.api import stub_fingerprint

_REASON = (
    "Ensembl gene annotations acquired manually; no committed baseline "
    "(a release-endpoint version probe is feasible but deferred — issue #177)"
)


def fetch_fingerprint() -> dict:
    return stub_fingerprint(_REASON)
