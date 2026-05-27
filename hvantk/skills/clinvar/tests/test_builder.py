"""Round-trip and update tests for the ClinVar builder skill."""

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

FIXTURE = "hvantk/skills/clinvar/tests/testdata/raw/clinvar/clinvar_20220403_chr20.vcf.bgz"
SNAPSHOT_DIR = Path("hvantk/skills/clinvar/tests/snapshots")

# (locus, alleles) is unique-in-table for ClinVar, so keys are inlined here
# rather than maintained in a separate sample_keys.json (per conventions §9).
SAMPLE_KEYS = [
    {"locus": "chr20:289563", "alleles": ["T", "C"]},
    {"locus": "chr20:8132666", "alleles": ["A", "T"]},
    {"locus": "chr20:21215630", "alleles": ["G", "T"]},
    {"locus": "chr20:34989865", "alleles": ["T", "C"]},
    {"locus": "chr20:44626468", "alleles": ["C", "T"]},
    {"locus": "chr20:50894815", "alleles": ["C", "T"]},
    {"locus": "chr20:63202482", "alleles": ["A", "G"]},
    {"locus": "chr20:63692857", "alleles": ["C", "G"]},
]


@pytest.mark.hail
def test_clinvar_round_trip(hail_session, tmp_path, regenerate_snapshots):
    """Build ClinVar from fixture; assert schema and sample-row stability."""
    import hail as hl
    from hvantk.skills.clinvar.builder import build_clinvar

    builder = phase_b_snapshot_adapter(build_clinvar, "clinvar:variants")
    builder_kwargs = {"reference_genome": "GRCh38"}

    if regenerate_snapshots:
        regenerate_snapshots_fn(
            builder_fn=builder,
            fixture_path=FIXTURE,
            snapshot_dir=SNAPSHOT_DIR,
            keys=SAMPLE_KEYS,
            builder_kwargs=builder_kwargs,
        )
        pytest.skip("Snapshots regenerated; rerun without --regenerate-snapshots to assert.")

    output_path = str(tmp_path / "clinvar.ht")
    builder(input_path=FIXTURE, output_path=output_path, **builder_kwargs)
    ht = hl.read_table(output_path)

    expected_schema = load_snapshot(SNAPSHOT_DIR / "schema.json")
    actual_schema = hail_schema_to_dict(ht)
    assert actual_schema == expected_schema, "ClinVar schema drifted from snapshot"

    expected_rows = load_snapshot(SNAPSHOT_DIR / "sample_rows.json")
    actual_rows = collect_sample_rows(ht, keys=SAMPLE_KEYS)
    assert actual_rows == expected_rows, "ClinVar sample rows drifted from snapshot"
