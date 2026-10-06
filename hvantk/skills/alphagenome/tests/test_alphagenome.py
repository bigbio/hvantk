"""Loader registration test for alphagenome:predictions (no Hail needed)."""

from __future__ import annotations

from hvantk.core.models import AnnotationTable
from hvantk.core.plugin import loader as plugin_loader


def test_alphagenome_predictions_registered():
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("alphagenome:predictions")
    assert spec.artifact_type is AnnotationTable
    assert spec.schema_id == "alphagenome-v2"
