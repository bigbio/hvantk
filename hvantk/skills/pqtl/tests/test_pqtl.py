"""Skipped: no fixture for pqtl; loader-only test (Phase K)."""
from __future__ import annotations

import pytest

from hvantk.core.models import AnnotationTable
from hvantk.core.plugin import loader as plugin_loader


def test_pqtl_metrics_registered():
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("pqtl:metrics")
    assert spec.artifact_type is AnnotationTable
    assert spec.schema_id == "pqtl-v1"


@pytest.mark.skip(
    reason="No fixture available for pqtl; manual smoke-test only"
)
def test_pqtl_metrics_round_trip():
    pass
