"""Phase B invariant: every loaded dataset must declare artifact_type and schema_id."""
from __future__ import annotations

from hvantk.core.plugin import loader as plugin_loader


def test_every_dataset_has_artifact_type_and_schema_id():
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    bad = []
    for ds in reg.list_datasets():
        if ds.artifact_type is None or ds.schema_id is None:
            bad.append(ds.name)
    assert not bad, (
        "Every dataset must declare artifact_type + schema_id in plugin.yaml. "
        f"Missing: {bad}"
    )
