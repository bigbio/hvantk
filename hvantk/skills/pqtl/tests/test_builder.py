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

# (locus, alleles, gene_id) keys for the synthetic fixture rows -- see
# testdata/raw/pqtl/README.md for the full row-purpose table. Picked to cover a
# mapped symbol (BRCA1), a negative-BETA/STAT mapped symbol (BRCA2, so SE must still
# come out positive), and an unmapped symbol that must fall back to its raw
# gene_name (SYNTHGENE1).
SAMPLE_KEYS = [
    {"locus": "chr17:43094687", "alleles": ["A", "G"], "gene_id": "ENSG00000012048"},
    {"locus": "chr13:32340073", "alleles": ["G", "C"], "gene_id": "ENSG00000139618"},
    {"locus": "chr5:100000", "alleles": ["A", "T"], "gene_id": "SYNTHGENE1"},
]


class _StubGeneCatalog:
    """Minimal gene_catalog double: maps only the symbols this fixture uses.

    Mirrors the ``GeneCatalogStreamer.map_ids`` contract via duck typing -- the
    builder only ever calls ``.map_ids``, and ``GeneCatalogStreamer`` is imported
    solely under ``TYPE_CHECKING`` in builder.py, so no subclassing is required.
    Returns a mapping restricted to the known symbols among those passed in; an
    unmapped symbol (``SYNTHGENE1``) is simply absent from the result, which is what
    drives the builder's ``hl.or_else(..., ht.gene_symbol)`` fallback-to-raw-symbol
    path.
    """

    _KNOWN = {
        "BRCA1": "ENSG00000012048",
        "BRCA2": "ENSG00000139618",
        "TP53": "ENSG00000141510",
    }

    def map_ids(self, symbols, source_type, target_type):
        assert source_type == "gene_symbol"
        assert target_type == "ensembl_gene_id"
        return {s: self._KNOWN[s] for s in symbols if s in self._KNOWN}


@pytest.mark.hail
def test_pqtl_metrics_round_trip(hail_session, tmp_path, regenerate_snapshots):
    """Build pQTL metrics from the synthetic Liver fixture; assert schema/row stability.

    Also pins two behaviors that have no dedicated test elsewhere (folded in here
    per the project's test-minimization policy rather than added as separate test
    functions):
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
