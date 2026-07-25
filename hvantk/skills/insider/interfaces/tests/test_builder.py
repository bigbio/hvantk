"""Round-trip and snapshot tests for the INSIDER per-protein interface builder."""

from __future__ import annotations

from pathlib import Path

import pytest

from hvantk.tests._snapshot_utils import (
    collect_sample_rows,
    hail_schema_to_dict,
    load_snapshot,
    phase_b_snapshot_adapter,
)

# Aliased to avoid shadowing the fixture name `regenerate_snapshots` in the test signature
from hvantk.tests._snapshot_utils import regenerate_snapshots as regenerate_snapshots_fn

FIXTURE = (
    "hvantk/skills/insider/interfaces/tests/testdata/raw/interfaces/"
    "H_sapiens_interfacesALL.txt"
)
SNAPSHOT_DIR = Path("hvantk/skills/insider/interfaces/tests/snapshots")

# The fixture's three proteins; the table is keyed on uniprot_id.
SAMPLE_KEYS = [
    {"uniprot_id": "Q00001"},
    {"uniprot_id": "Q00002"},
    {"uniprot_id": "Q00003"},
]


@pytest.mark.hail
def test_insider_interfaces_round_trip(hail_session, tmp_path, regenerate_snapshots):
    """Build the interface reduction from fixture; assert schema + sample-row stability."""
    import hail as hl
    from hvantk.skills.insider.interfaces.builder import build_insider_interfaces

    builder = phase_b_snapshot_adapter(build_insider_interfaces, "insider:interfaces")

    if regenerate_snapshots:
        regenerate_snapshots_fn(
            builder_fn=builder,
            fixture_path=FIXTURE,
            snapshot_dir=SNAPSHOT_DIR,
            keys=SAMPLE_KEYS,
            builder_kwargs={},
        )
        pytest.skip(
            "Snapshots regenerated; rerun without --regenerate-snapshots to assert."
        )

    output_path = str(tmp_path / "insider_interfaces.ht")
    fixture_uri = Path(FIXTURE).resolve().as_uri()
    builder(input_path=fixture_uri, output_path=output_path)
    ht = hl.read_table(output_path)

    expected_schema = load_snapshot(SNAPSHOT_DIR / "schema.json")
    assert hail_schema_to_dict(ht) == expected_schema, "interfaces schema drifted"

    expected_rows = load_snapshot(SNAPSHOT_DIR / "sample_rows.json")
    assert (
        collect_sample_rows(ht, keys=SAMPLE_KEYS) == expected_rows
    ), "interfaces sample rows drifted"


@pytest.mark.hail
def test_empty_interface_table_is_rejected(hail_session, tmp_path):
    """A header-only file means the parse found nothing -- fail loud, don't ship 0 genes."""
    from hvantk.skills.insider.interfaces.builder import build_insider_interfaces

    src = tmp_path / "empty.txt"
    src.write_text("P1\tP2\tSource\tP1_IRES\tP2_IRES\n")

    class _Ctx:
        def provenance(self, **kw):
            from hvantk.core.models import Provenance

            return Provenance(schema_id=kw.get("schema_id", "insider-interfaces-v1"))

    with pytest.raises(ValueError, match="0 proteins"):
        build_insider_interfaces(str(src), _Ctx())
