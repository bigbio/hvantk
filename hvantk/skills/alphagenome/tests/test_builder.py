"""Round-trip test for the AlphaGenome tidy-scores builder.

The fixture is a 150-row real subset of an AlphaGenome ``tidy_scores()`` shard,
licensed under the AlphaGenome Output Terms (non-commercial; see
``testdata/raw/alphagenome/NOTICE.md``).
"""

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
FIXTURE_DIR = str(_TESTS_DIR / "testdata/raw/alphagenome")
SNAPSHOT_DIR = _TESTS_DIR / "snapshots"

# The fixture's three variants; the output has exactly one row per variant.
SAMPLE_KEYS = [
    {"locus": "chr3:39408741", "alleles": ["T", "C"]},
    {"locus": "chr6:112216367", "alleles": ["C", "A"]},
    {"locus": "chrX:153694448", "alleles": ["T", "G"]},
]

# The oracle values below (AG_SPLICE_SITES, HEART_SPLICE_SITE_USAGE) are AlphaGenome
# outputs: AlphaGenome Output Terms of Use, non-commercial, not MIT-licensed (see
# testdata/raw/alphagenome/NOTICE.md).
#
# Oracle: ag_splice_sites = max |raw_score| over a variant's SPLICE_SITES rows,
# from an independent per-variant ClinVar benchmark computed by the maintainer,
# read at float32 on 2026-10-06. Over the full 411-variant shard the builder
# reproduces it for all 339 variants that benchmark covers.
AG_SPLICE_SITES = {
    "chr3:39408741:T>C": 0.71823883,
    "chr6:112216367:C>A": 0.9982147,
    "chrX:153694448:T>G": 0.9948883,
}

# Heart build: splice_site_usage.max_abs_raw with ontology_curies=HEART_CURIES.
# These values are curie-defined. UBERON:0006566 + UBERON:0006631 select three
# SPLICE_SITE_USAGE tracks: the two GTEx heart tracks and the ENCODE
# `usage_UBERON:0006631 total RNA-seq` track. Cross-checked pandas vs Hail on the
# full 411-variant shard (2026-10-06; every summary of the heart build equal).
# They intentionally differ from that benchmark's GTEx-only ag_heart_ssu for 2 of 3
# variants (chr6: 0.9678116, chrX: 0.91789865), where the ENCODE track scores
# higher; the builder has no gtex_tissue filter to reproduce that feature.
HEART_CURIES = ["UBERON:0006566", "UBERON:0006631"]
HEART_SPLICE_SITE_USAGE = {
    "chr3:39408741:T>C": 0.00390625,
    "chr6:112216367:C>A": 0.98340607,
    "chrX:153694448:T>G": 0.93353176,
}


@pytest.mark.hail
def test_alphagenome_round_trip(hail_session, tmp_path, regenerate_snapshots):
    """Build the real-subset fixture; check snapshots, the oracle and a heart build."""
    import hail as hl
    from hvantk.skills.alphagenome.builder import build_alphagenome_predictions

    builder = phase_b_snapshot_adapter(
        build_alphagenome_predictions, "alphagenome:predictions"
    )

    if regenerate_snapshots:
        regenerate_snapshots_fn(
            builder_fn=builder,
            fixture_path=FIXTURE_DIR,
            snapshot_dir=SNAPSHOT_DIR,
            keys=SAMPLE_KEYS,
        )
        pytest.skip(
            "Snapshots regenerated; rerun without --regenerate-snapshots to assert."
        )

    output_path = str(tmp_path / "alphagenome.ht")
    builder(input_path=FIXTURE_DIR, output_path=output_path)
    ht = hl.read_table(output_path)
    from hvantk.core.io import load

    assert load(output_path).provenance.build_parameters == {
        "output_types": None,
        "ontology_curies": None,
    }

    expected_schema = load_snapshot(SNAPSHOT_DIR / "schema.json")
    assert hail_schema_to_dict(ht) == expected_schema, (
        "AlphaGenome schema drifted from snapshot"
    )
    expected_rows = load_snapshot(SNAPSHOT_DIR / "sample_rows.json")
    assert collect_sample_rows(ht, keys=SAMPLE_KEYS) == expected_rows, (
        "AlphaGenome sample rows drifted from snapshot"
    )

    rows = {r.variant_id: r for r in ht.collect()}
    splice_sites = {v: rows[v].splice_sites.max_abs_raw for v in AG_SPLICE_SITES}
    assert splice_sites == pytest.approx(AG_SPLICE_SITES, rel=1e-6)

    heart_path = str(tmp_path / "alphagenome_heart.ht")
    # Given in reverse order: provenance records the filter values sorted.
    builder(
        input_path=FIXTURE_DIR,
        output_path=heart_path,
        ontology_curies=HEART_CURIES[::-1],
    )
    assert load(heart_path).provenance.build_parameters == {
        "output_types": None,
        "ontology_curies": HEART_CURIES,
    }
    heart = {
        r.variant_id: r.splice_site_usage.max_abs_raw
        for r in hl.read_table(heart_path).collect()
    }
    assert heart == pytest.approx(HEART_SPLICE_SITE_USAGE, rel=1e-6)

    filtered_path = str(tmp_path / "alphagenome_rna.ht")
    builder(input_path=FIXTURE_DIR, output_path=filtered_path, output_types="RNA_SEQ")
    assert load(filtered_path).provenance.build_parameters == {
        "output_types": ["RNA_SEQ"],
        "ontology_curies": None,
    }
