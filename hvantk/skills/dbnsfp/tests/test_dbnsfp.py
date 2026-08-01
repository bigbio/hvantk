"""Conformance test for the dbnsfp plugin (Phase K)."""
from __future__ import annotations

from pathlib import Path

import pytest

from hvantk.core.models import AnnotationTable
from hvantk.core.plugin import loader as plugin_loader
from hvantk.core.plugin.run_builder import run_builder_for_spec


def test_dbnsfp_variants_registered():
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("dbnsfp:variants")
    assert spec.artifact_type is AnnotationTable
    assert spec.schema_id == "dbnsfp-v1"


@pytest.mark.hail
def test_dbnsfp_variants_round_trip(tmp_path):
    plugin_loader.reset_registry_for_tests()
    reg = plugin_loader.get_registry()
    spec = reg.get_dataset("dbnsfp:variants")

    assert spec.artifact_type is AnnotationTable
    assert spec.schema_id == "dbnsfp-v1"

    fixture = Path(
        "hvantk/tests/testdata/raw/dbnsfp/dbNSFP4_v49a_example_variants.bgz"
    )
    assert fixture.exists(), f"Fixture not found: {fixture}"

    out = tmp_path / "variants.ht"
    prov = run_builder_for_spec(
        spec,
        parsed_input=fixture,
        output_path=out,
        plugin_version=spec.plugin_version,
    )

    assert prov.plugin == "dbnsfp"
    assert prov.schema_id == "dbnsfp-v1"

    from hvantk.core import io as core_io
    loaded = core_io.load(out)
    assert isinstance(loaded, AnnotationTable)
    assert loaded.backend == "hail"
    assert loaded.count() > 0


@pytest.mark.hail
def test_rankscores_are_numeric_not_strings(tmp_path):
    """dbNSFP's 57 ``*_rankscore`` columns must land as float64, not str.

    A rankscore is ONE value per variant (rank-normalised across the whole database),
    unlike the ``*_score`` columns, which are per-transcript dicts. Leaving them as raw
    strings makes them unusable as features: ``hl.agg.mean`` on a str fails, and any
    downstream aggregation would need a per-consumer cast.

    Rankscores matter because they are the only dbNSFP surface where every predictor is
    comparable -- rank-normalised to 0-1 AND direction-normalised, which is what the
    ``_converted_`` prefix means (SIFT and FATHMM are inverted so that, for every
    rankscore without exception, HIGHER = more damaging). A wide predictor axis built on
    raw scores would be silently backwards for those tools.
    """
    import hail as hl

    plugin_loader.reset_registry_for_tests()
    spec = plugin_loader.get_registry().get_dataset("dbnsfp:variants")
    fixture = Path("hvantk/tests/testdata/raw/dbnsfp/dbNSFP4_v49a_example_variants.bgz")
    out = tmp_path / "variants.ht"
    run_builder_for_spec(spec, parsed_input=fixture, output_path=out,
                         plugin_version=spec.plugin_version)

    ht = hl.read_table(str(out))
    rank = [f for f in ht.row for _ in (0,) if f.endswith("_rankscore")]
    assert len(rank) >= 50, f"expected dbNSFP's ~57 rankscore fields, found {len(rank)}"
    wrong = [f for f in rank if ht[f].dtype != hl.tfloat64]
    assert not wrong, f"rankscore fields still non-float: {wrong[:5]} ({len(wrong)} total)"

    # "." must become missing, not 0.0 -- a missing predictor is not a benign one.
    row = ht.filter(hl.is_defined(ht.MutationTaster_converted_rankscore)).head(1).collect()
    assert row, "no row with a defined MutationTaster_converted_rankscore"
    assert 0.0 <= row[0].MutationTaster_converted_rankscore <= 1.0
