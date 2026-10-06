"""Round-trip test for the pQTL (Fang et al.) builder skill."""

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
FIXTURE_DIR = str(_TESTS_DIR / "testdata/raw/pqtl")
SNAPSHOT_DIR = _TESTS_DIR / "snapshots"

# (locus, alleles, gene_id) keys for all 7 surviving fixture rows (the 8th, with
# STAT == 0, is dropped) -- see testdata/raw/pqtl/README.md for the full
# row-purpose table. Every surviving row is pinned here, not just a subset, so the
# snapshot actually proves each edge case transformed correctly rather than merely
# surviving: a mapped symbol (BRCA1); a negative-BETA/STAT mapped symbol (BRCA2, so
# SE must still come out positive); an unmapped symbol that must fall back to its
# raw gene_name (SYNTHGENE1); a second gene (TP53) at the *same* locus/alleles as
# the BRCA1 row -- the key stays unique only because gene_id differs; a chrX
# variant; an indel id (multi-base REF/ALT); and an extreme p-value row.
SAMPLE_KEYS = [
    {"locus": "chr17:43094687", "alleles": ["A", "G"], "gene_id": "ENSG00000012048"},
    {"locus": "chr13:32340073", "alleles": ["G", "C"], "gene_id": "ENSG00000139618"},
    {"locus": "chr5:100000", "alleles": ["A", "T"], "gene_id": "SYNTHGENE1"},
    {"locus": "chr17:43094687", "alleles": ["A", "G"], "gene_id": "ENSG00000141510"},
    {"locus": "chrX:71130000", "alleles": ["C", "T"], "gene_id": "ENSG00000139618"},
    {"locus": "chr1:1000000", "alleles": ["AT", "A"], "gene_id": "ENSG00000141510"},
    {"locus": "chr17:43110000", "alleles": ["T", "C"], "gene_id": "ENSG00000012048"},
]


class _StubGeneCatalog:
    """Minimal gene_catalog double: maps only the symbols this fixture knows.

    Mirrors the ``GeneCatalogStreamer.map_ids`` contract via duck typing -- the
    builder only ever calls ``.map_ids``, and ``GeneCatalogStreamer`` is imported
    solely under ``TYPE_CHECKING`` in builder.py, so no subclassing is required.
    Per that contract (``hvantk/core/streamers/gene_catalog.py``), the returned
    dict has the *same keys as the input* -- an unmapped symbol (``SYNTHGENE1``)
    maps to ``None`` rather than being omitted. This does not change which
    ``gene_id`` any row gets: the builder's
    ``hl.or_else(mapping_literal.get(ht.gene_symbol), ht.gene_symbol)`` falls back
    to the raw symbol identically whether the key is missing or present with a
    ``None``/missing value.
    """

    _KNOWN = {
        "BRCA1": "ENSG00000012048",
        "BRCA2": "ENSG00000139618",
        "TP53": "ENSG00000141510",
    }

    def map_ids(self, symbols, source_type, target_type):
        assert source_type == "gene_symbol"
        assert target_type == "ensembl_gene_id"
        return {s: self._KNOWN.get(s) for s in symbols}


@pytest.mark.hail
def test_pqtl_metrics_round_trip(hail_session, tmp_path, regenerate_snapshots):
    """Build pQTL metrics from the synthetic Liver fixture; assert schema/row stability.

    Also pins two behaviors that have no dedicated test elsewhere, folded into this
    test rather than split into separate ones so they reuse this test's single Hail
    build instead of each paying for their own:
      - the fixture's STAT == 0 row (gene_name=BRCA1, chr17:43095211) is dropped,
        since SE is undefined for it;
      - every surviving row's derived SE (= |beta / stat|) is strictly positive.
    """
    import hail as hl
    from hvantk.skills.pqtl.builder import build_pqtl_metrics

    builder = phase_b_snapshot_adapter(build_pqtl_metrics, "pqtl:metrics")
    builder_kwargs = {
        "reference_genome": "GRCh38",
        "tissue": "Liver",
        "gene_catalog": _StubGeneCatalog(),
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

    output_path = str(tmp_path / "pqtl_metrics.ht")
    builder(input_path=FIXTURE_DIR, output_path=output_path, **builder_kwargs)
    ht = hl.read_table(output_path)

    expected_schema = load_snapshot(SNAPSHOT_DIR / "schema.json")
    actual_schema = hail_schema_to_dict(ht)
    assert actual_schema == expected_schema, "pQTL schema drifted from snapshot"

    expected_rows = load_snapshot(SNAPSHOT_DIR / "sample_rows.json")
    actual_rows = collect_sample_rows(ht, keys=SAMPLE_KEYS)
    assert actual_rows == expected_rows, "pQTL sample rows drifted from snapshot"

    # The fixture has 8 rows and exactly one (STAT == 0) must be dropped; catches a
    # row silently dropped or duplicated among rows the snapshot does not sample.
    assert ht.count() == 7, "expected 8 fixture rows minus 1 dropped STAT==0 row"

    # The fixture's STAT == 0 row (BRCA1 @ chr17:43095211) must be dropped by the
    # builder before SE derivation (|beta / stat|  is undefined at stat == 0).
    dropped_row_count = ht.filter(
        (ht.locus == hl.locus("chr17", 43095211, reference_genome="GRCh38"))
        & (ht.alleles == ["A", "G"])
        & (ht.gene_id == "ENSG00000012048")
    ).count()
    assert dropped_row_count == 0, "STAT == 0 row must be dropped by the builder"

    # Every surviving row's derived SE must be strictly positive.
    min_se = ht.aggregate(hl.agg.min(ht.se))
    assert min_se is not None and min_se > 0, f"expected every se > 0, min was {min_se}"
