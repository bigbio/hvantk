"""Drift probe for pqtl — documentation-only source (stub).

pQTL summary statistics are publication-only (Fang et al. 2025 / GTEx); there is
no stable direct data URL to fingerprint. This probe returns a structured stub
sentinel so ``hvantk drift`` reports a visible WARNING (status="stub") rather
than a silent false-green. Replace with a real probe if a direct data URL becomes
available. See issue #177.
"""
from __future__ import annotations

from hvantk.core.plugin.api import stub_fingerprint

_REASON = (
    "pQTL summary statistics are publication-only (Fang et al. 2025 / GTEx); "
    "no stable direct data URL to fingerprint"
)


def fetch_fingerprint() -> dict:
    return stub_fingerprint(_REASON)
