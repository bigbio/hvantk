"""Drift probe for dbnsfp — documentation-only source (stub).

dbNSFP is distributed from a landing page (sites.google.com/site/jpopgen/dbNSFP)
with no stable direct data URL to fingerprint. This probe returns a structured
stub sentinel so ``hvantk drift`` reports a visible WARNING (status="stub")
rather than a silent false-green. Replace with a real probe if a direct data URL
becomes available. See issue #177.
"""
from __future__ import annotations

from hvantk.core.plugin.api import stub_fingerprint

_REASON = (
    "dbNSFP is distributed from a landing page "
    "(sites.google.com/site/jpopgen/dbNSFP); no stable direct data URL"
)


def fetch_fingerprint() -> dict:
    return stub_fingerprint(_REASON)
