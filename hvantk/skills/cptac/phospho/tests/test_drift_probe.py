"""Sanity test for the CPTAC phospho drift probe.

The probe reports the installed ``cptac`` package version. We mock
``importlib.metadata.version`` to assert the fingerprint shape regardless
of whether ``cptac`` is installed in the dev environment.
"""

from unittest import mock


def test_fingerprint_shape():
    from hvantk.skills.cptac.phospho import drift_probe

    with mock.patch.object(drift_probe, "_cptac_version", return_value="1.5.13"):
        fp = drift_probe.fetch_fingerprint()

    assert fp["probe_version"] == drift_probe.PROBE_VERSION
    assert fp["source_version"] == "1.5.13"
    assert fp["headers"] == {"cptac": {"package_version": "1.5.13"}}
    assert "fetched_at" in fp


def test_fingerprint_raises_when_package_missing():
    import pytest

    from hvantk.core.plugin_api import DriftProbeError
    from hvantk.skills.cptac.phospho import drift_probe

    with mock.patch.object(drift_probe, "_cptac_version", return_value=None):
        with pytest.raises(DriftProbeError):
            drift_probe.fetch_fingerprint()
