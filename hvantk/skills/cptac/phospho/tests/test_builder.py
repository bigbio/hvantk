"""Unit tests for the CPTAC phospho Phase B builder (#198).

These exercise the builder without the upstream ``cptac`` package: they feed
matrix/metadata CSVs in the same layout ``CPTACPhosphoDataset.download`` writes
(matrix index name ``Site``; metadata index name ``SampleID``).
"""
from __future__ import annotations

import pandas as pd
import pytest


def _ctx():
    from hvantk.core.models import BuildContext

    return BuildContext(
        plugin="cptac",
        dataset="cptac:phospho",
        plugin_version="0.1.0",
        source_fingerprint="test",
        builder_commit=None,
    )


def _write_fixture(raw, ct="brca"):
    pd.DataFrame(
        {"S1": [1.0, 2.0], "S2": [3.0, 4.0]},
        index=pd.Index(["TP53_S15", "EGFR_Y1068"], name="Site"),
    ).to_csv(raw / f"cptac-phospho-{ct}-matrix.csv")
    pd.DataFrame(
        {"stage": ["I", "II"]},
        index=pd.Index(["S1", "S2"], name="SampleID"),
    ).to_csv(raw / f"cptac-phospho-{ct}-metadata.csv")


def test_build_from_raw_dir(tmp_path):
    from hvantk.core.models import ExpressionMatrix
    from hvantk.skills.cptac.phospho.builder import build_cptac_phospho

    raw = tmp_path / "raw"
    raw.mkdir()
    _write_fixture(raw, "brca")

    art = build_cptac_phospho(str(raw), _ctx(), cancer_type="brca")

    assert isinstance(art, ExpressionMatrix)
    adata = art.to_anndata()
    assert adata.shape == (2, 2)  # 2 samples x 2 sites
    assert set(adata.obs.index) == {"S1", "S2"}
    assert set(adata.var.index) == {"TP53_S15", "EGFR_Y1068"}
    assert "stage" in adata.obs.columns


def test_build_from_dict(tmp_path):
    from hvantk.core.models import ExpressionMatrix
    from hvantk.skills.cptac.phospho.builder import build_cptac_phospho

    raw = tmp_path / "raw"
    raw.mkdir()
    _write_fixture(raw, "brca")

    art = build_cptac_phospho(
        {
            "expression": str(raw / "cptac-phospho-brca-matrix.csv"),
            "metadata": str(raw / "cptac-phospho-brca-metadata.csv"),
        },
        _ctx(),
    )
    assert isinstance(art, ExpressionMatrix)
    assert art.to_anndata().shape == (2, 2)


def test_raw_dir_requires_cancer_type(tmp_path):
    from hvantk.skills.cptac.phospho.builder import build_cptac_phospho

    raw = tmp_path / "raw"
    raw.mkdir()
    _write_fixture(raw, "brca")

    with pytest.raises(ValueError, match="cancer_type"):
        build_cptac_phospho(str(raw), _ctx())


def test_missing_files_raise(tmp_path):
    from hvantk.skills.cptac.phospho.builder import build_cptac_phospho

    raw = tmp_path / "raw"
    raw.mkdir()
    with pytest.raises(FileNotFoundError):
        build_cptac_phospho(str(raw), _ctx(), cancer_type="brca")


def test_ignores_stray_kwargs(tmp_path):
    from hvantk.skills.cptac.phospho.builder import build_cptac_phospho

    raw = tmp_path / "raw"
    raw.mkdir()
    _write_fixture(raw, "brca")

    # reprocess forwards download-only --plugin-arg values (e.g. overwrite) to the
    # builder too; they must not break the build.
    art = build_cptac_phospho(str(raw), _ctx(), cancer_type="brca", overwrite=True)
    assert art.to_anndata().shape == (2, 2)
