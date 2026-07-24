"""Sanity test for the INSIDER stub drift probe.

The INSIDER source has no programmatically-probeable URL, so the probe returns
a structured stub sentinel (status="stub" + reason) instead of a false-green
constant. This asserts that shape so a regression back to the silent
false-green cannot slip through. See issue #177.
"""

from hvantk.core.plugin.api import PROBE_STATUS_STUB, STUB_FINGERPRINT_TOKEN
from hvantk.skills.insider.drift_probe import (
    NOT_IMPLEMENTED_REASON,
    fetch_fingerprint,
)


def test_stub_fingerprint_shape():
    fp = fetch_fingerprint()
    assert fp["probe_status"] == PROBE_STATUS_STUB
    assert fp["fingerprint"] == STUB_FINGERPRINT_TOKEN
    assert fp["reason"] == NOT_IMPLEMENTED_REASON
    assert "not implemented" in fp["reason"]
