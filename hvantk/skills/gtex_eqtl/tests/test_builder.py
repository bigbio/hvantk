"""Round-trip and update tests for the GTEx v11 eQTL builder skill."""

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

_TESTS_DIR = Path(__file__).parent
FIXTURE_FILE = _TESTS_DIR / "testdata/raw/gtex-eqtl/Liver.v11.eQTLs.signif_pairs.parquet"
FIXTURE_DIR = str(FIXTURE_FILE.parent)
SNAPSHOT_DIR = _TESTS_DIR / "snapshots"

# (locus, alleles, gene_id) is unique-in-table for v11 signif_pairs (per skill §5),
# so keys are inlined here rather than maintained in a separate sample_keys.json
# (per conventions §9). Keys picked from observed fixture rows after gene-version
# stripping (ENSG00000268903.1 -> ENSG00000268903).
SAMPLE_KEYS = [
    {"locus": "chr1:63671", "alleles": ["G", "A"], "gene_id": "ENSG00000268903"},
    {"locus": "chr1:282564", "alleles": ["C", "A"], "gene_id": "ENSG00000268903"},
    {"locus": "chr1:64649", "alleles": ["A", "C"], "gene_id": "ENSG00000308579"},
]


@pytest.mark.hail
def test_gtex_eqtl_round_trip(hail_session, tmp_path, regenerate_snapshots):
    """Build GTEx v11 eQTL from Liver fixture; assert schema and sample-row stability."""
    import hail as hl
    from hvantk.skills.gtex_eqtl.builder import build_eqtl_associations

    builder = phase_b_snapshot_adapter(build_eqtl_associations, "gtex_eqtl:eqtls")
    builder_kwargs = {
        "reference_genome": "GRCh38",
        "source": "gtex_v11",
        "tissue": "Liver",
        "p_threshold": 0,  # retain all fixture rows (signif_pairs is pre-filtered)
    }

    if regenerate_snapshots:
        regenerate_snapshots_fn(
            builder_fn=builder,
            fixture_path=FIXTURE_DIR,
            snapshot_dir=SNAPSHOT_DIR,
            keys=SAMPLE_KEYS,
            builder_kwargs=builder_kwargs,
        )
        pytest.skip(
            "Snapshots regenerated; rerun without --regenerate-snapshots to assert."
        )

    output_path = str(tmp_path / "gtex_eqtl.ht")
    builder(input_path=FIXTURE_DIR, output_path=output_path, **builder_kwargs)
    ht = hl.read_table(output_path)

    expected_schema = load_snapshot(SNAPSHOT_DIR / "schema.json")
    actual_schema = hail_schema_to_dict(ht)
    assert actual_schema == expected_schema, "GTEx eQTL schema drifted from snapshot"

    expected_rows = load_snapshot(SNAPSHOT_DIR / "sample_rows.json")
    actual_rows = collect_sample_rows(ht, keys=SAMPLE_KEYS)
    assert actual_rows == expected_rows, "GTEx eQTL sample rows drifted from snapshot"
