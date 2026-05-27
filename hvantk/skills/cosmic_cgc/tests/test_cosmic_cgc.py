"""Skipped: no fixture for cosmic-cgc; loader-only test (Phase K)."""
from __future__ import annotations

import pytest

from hvantk.core.models import AnnotationTable
from hvantk.core.plugin import loader as plugin_loader


def test_cosmic_cgc_submissions_registered():
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("cosmic-cgc:submissions")
    assert spec.artifact_type is AnnotationTable
    assert spec.schema_id == "cosmic-cgc-v1"


@pytest.mark.skip(
    reason="No fixture available for cosmic-cgc (COSMIC requires account login); manual smoke-test only"
)
def test_cosmic_cgc_submissions_round_trip():
    pass
