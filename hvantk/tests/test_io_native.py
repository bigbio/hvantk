"""Tests for core/io.load_native and save_native — Hail-native passthrough."""
from __future__ import annotations

import dataclasses
from datetime import datetime, timezone

import pandas as pd
import pytest

from hvantk.core import io as core_io
from hvantk.core.models import (
    AnnotationTable,
    ExpressionMatrix,
    GeneSet,
    Provenance,
)


def _prov(schema_id: str = "t-rows-v1") -> Provenance:
    return Provenance(
        plugin="t",
        dataset="t:rows",
        plugin_version="0",
        source_fingerprint="sha256:x",
        schema_id=schema_id,
        build_timestamp=datetime(2026, 5, 21, tzinfo=timezone.utc),
        builder_commit=None,
    )


def test_load_native_parquet_returns_dataframe_and_provenance(tmp_path):
    df = pd.DataFrame({"x": [1, 2, 3]})
    ann = AnnotationTable.from_pandas(df, provenance=_prov())
    out = tmp_path / "rows.parquet"
    core_io.save(ann, out)

    native, prov = core_io.load_native(out)
    assert isinstance(native, pd.DataFrame)
    assert native["x"].tolist() == [1, 2, 3]
    assert prov == ann.provenance


def test_save_native_parquet_roundtrip(tmp_path):
    """A raw pandas DataFrame can be saved + loaded via the native API."""
    df = pd.DataFrame({"gene": ["BRCA1", "TP53"], "score": [0.7, 0.9]})
    out = tmp_path / "rows.parquet"
    core_io.save_native(df, out, provenance=_prov())

    # Sidecar exists; the artifact-typed load() still works on the same file.
    assert (out.with_name(out.name + ".provenance.json")).exists()
    loaded_ann = core_io.load(out)
    assert isinstance(loaded_ann, AnnotationTable)
    out_df = loaded_ann.to_pandas()
    assert out_df["gene"].tolist() == ["BRCA1", "TP53"]


def test_load_native_geneset_returns_list(tmp_path):
    gs = GeneSet(
        name="brca", provenance=_prov(), _members=frozenset({"BRCA1", "BRCA2"})
    )
    out = tmp_path / "brca.geneset.json"
    core_io.save(gs, out)

    native, prov = core_io.load_native(out)
    assert isinstance(native, list)
    assert set(native) == {"BRCA1", "BRCA2"}
    assert prov == gs.provenance


def test_save_native_geneset_from_iterable(tmp_path):
    out = tmp_path / "kinases.geneset.json"
    core_io.save_native(["AURKA", "AURKB", "PLK1"], out, provenance=_prov())

    loaded = core_io.load(out)
    assert isinstance(loaded, GeneSet)
    assert loaded.to_set() == {"AURKA", "AURKB", "PLK1"}


def test_load_save_native_chains_provenance(tmp_path):
    """A Hail-native algorithm pattern: load native, derive, save with parents."""
    src_df = pd.DataFrame({"x": [1, 2, 3, 4]})
    src_ann = AnnotationTable.from_pandas(src_df, provenance=_prov())
    src_path = tmp_path / "src.parquet"
    core_io.save(src_ann, src_path)

    df, src_prov = core_io.load_native(src_path)
    derived = df[df["x"] > 2].reset_index(drop=True)

    derived_prov = dataclasses.replace(
        _prov(schema_id="derived-v1"), parents=(src_prov,)
    )
    derived_path = tmp_path / "derived.parquet"
    core_io.save_native(derived, derived_path, provenance=derived_prov)

    loaded_derived = core_io.load(derived_path)
    assert src_prov in loaded_derived.provenance.parents


def test_save_native_unrecognized_extension_raises(tmp_path):
    from hvantk.core.io._errors import ArtifactTypeError

    df = pd.DataFrame({"x": [1]})
    with pytest.raises(ArtifactTypeError, match="unrecognized extension"):
        core_io.save_native(df, tmp_path / "rows.csv", provenance=_prov())


def test_load_native_h5ad_returns_anndata(tmp_path):
    import anndata as ad
    import numpy as np

    obs = pd.DataFrame({"tissue": ["liver"]}, index=["s1"])
    var = pd.DataFrame({"gene": ["BRCA1"]}, index=["g1"])
    adata = ad.AnnData(X=np.array([[1.0]]), obs=obs, var=var)
    em = ExpressionMatrix.from_anndata(adata, provenance=_prov("expr-v1"))
    out = tmp_path / "expr.h5ad"
    core_io.save(em, out)

    native, prov = core_io.load_native(out)
    assert isinstance(native, ad.AnnData)
    assert native.n_obs == 1
    assert native.n_vars == 1
    assert prov.schema_id == "expr-v1"


@pytest.mark.hail
def test_load_native_ht_returns_hail_table(tmp_path):
    """Zero-cost native passthrough for Hail-backed artifacts."""
    import hail as hl

    ht = hl.Table.parallelize(
        [{"gene": "BRCA1", "score": 0.7}],
        hl.tstruct(gene=hl.tstr, score=hl.tfloat64),
    )
    ann = AnnotationTable.from_hail(ht, provenance=_prov())
    out = tmp_path / "rows.ht"
    core_io.save(ann, out)

    native, prov = core_io.load_native(out)
    assert isinstance(native, hl.Table)
    assert native.count() == 1
    assert prov == ann.provenance


@pytest.mark.hail
def test_save_native_ht_roundtrip(tmp_path):
    import hail as hl

    ht = hl.Table.parallelize(
        [{"x": 1}, {"x": 2}],
        hl.tstruct(x=hl.tint32),
    )
    out = tmp_path / "rows.ht"
    core_io.save_native(ht, out, provenance=_prov(schema_id="hail-rows-v1"))

    loaded = core_io.load(out)
    assert isinstance(loaded, AnnotationTable)
    assert loaded.backend == "hail"
    assert loaded.count() == 2
    assert loaded.provenance.schema_id == "hail-rows-v1"
