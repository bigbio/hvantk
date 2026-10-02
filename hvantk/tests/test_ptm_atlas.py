"""``build_atlas`` reports the sources that actually reached the pipeline.

The pipeline core is stubbed out: these tests check only how ``build_atlas`` translates
its config, so they need no GTF, no network and no Hail.
"""

import pytest

from hvantk.algorithms.ptm import atlas
from hvantk.algorithms.ptm.atlas import PTMAtlasConfig, build_atlas
from hvantk.algorithms.ptm.pipeline import PTMBuildResult


@pytest.fixture
def captured(monkeypatch):
    seen = {}

    def fake_core(build_cfg):
        seen["cfg"] = build_cfg
        return PTMBuildResult(n_mapped=3, mapped_tsv_path="mapped.tsv.bgz")

    monkeypatch.setattr(atlas, "ptm_build_pipeline_core", fake_core)
    return seen


def _config(tmp_path, **kw):
    uniprot = tmp_path / "uniprot.tsv"
    uniprot.write_text("x\n")
    return PTMAtlasConfig(
        output_dir=str(tmp_path / "out"),
        output_ht=str(tmp_path / "out" / "sites.ht"),
        uniprot_tsv=str(uniprot),
        **kw,
    )


def test_a_requested_source_without_its_tsv_is_not_reported_as_used(tmp_path, captured):
    # peptideatlas is in the default sources, but no TSV is passed, so the core gets
    # None for it and skips it. Reporting it as used claimed data the atlas never held.
    result = build_atlas(_config(tmp_path, sources=["uniprot", "peptideatlas"]))
    assert captured["cfg"].peptideatlas_tsv is None
    assert result.sources_used == ["uniprot"]


def test_a_source_with_its_tsv_is_reported_as_used(tmp_path, captured):
    pa = tmp_path / "peptideatlas.tsv"
    pa.write_text("x\n")
    result = build_atlas(
        _config(tmp_path, sources=["uniprot", "peptideatlas"], peptideatlas_tsv=str(pa))
    )
    assert captured["cfg"].peptideatlas_tsv == str(pa)
    assert result.sources_used == ["uniprot", "peptideatlas"]


def test_an_unselected_source_is_not_used_even_when_its_tsv_is_passed(
    tmp_path, captured
):
    cptac = tmp_path / "cptac.tsv"
    cptac.write_text("x\n")
    result = build_atlas(_config(tmp_path, sources=["uniprot"], cptac_tsv=str(cptac)))
    assert captured["cfg"].cptac_tsv is None
    assert result.sources_used == ["uniprot"]
