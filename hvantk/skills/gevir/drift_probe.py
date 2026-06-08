"""Drift probe for gevir — documentation-only source (stub).

GeVIR metrics are published as supplementary data (PMID 31873297); there is no
programmatic data URL to fingerprint. This probe returns a structured stub
sentinel so ``hvantk drift`` reports a visible WARNING (status="stub") rather
than a silent false-green. Replace with a real probe if a direct data URL becomes
available. See issue #177.
"""
from __future__ import annotations

from hvantk.core.plugin.api import stub_fingerprint

_REASON = (
    "GeVIR metrics are published as supplementary data (PMID 31873297); "
    "no programmatic data URL to fingerprint"
)


def fetch_fingerprint() -> dict:
    return stub_fingerprint(_REASON)
