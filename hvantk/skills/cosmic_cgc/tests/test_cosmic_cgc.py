"""Loader registration test for cosmic-cgc. The round-trip test lives in
test_builder.py, against the committed synthetic fixture."""

from __future__ import annotations

from hvantk.core.models import AnnotationTable
from hvantk.core.plugin import loader as plugin_loader


def test_cosmic_cgc_submissions_registered():
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("cosmic-cgc:submissions")
    assert spec.artifact_type is AnnotationTable
    assert spec.schema_id == "cosmic-cgc-v2"
    assert spec.schema_ids == ("cosmic-cgc-legacy-v1",)
