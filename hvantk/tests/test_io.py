"""Tests for core/io: save / load round-trip, sidecar manifest, dispatch by ext."""
from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from hvantk.core import io as core_io
from hvantk.core.io._errors import ArtifactTypeError, SchemaIdMismatchError
from hvantk.core.models._expr import col
from hvantk.core.models.annotation_table import AnnotationTable
from hvantk.core.models.expression_matrix import ExpressionMatrix
from hvantk.core.models.provenance import Provenance


def _prov(schema_id="t-rows-v1"):
    return Provenance(
        plugin="t",
        dataset="t:rows",
        plugin_version="0.0",
        source_fingerprint="sha256:abc",
        schema_id=schema_id,
        build_timestamp=datetime(2026, 5, 20, tzinfo=timezone.utc),
        builder_commit=None,
    )


def test_save_load_pandas_round_trip(tmp_path):
    df = pd.DataFrame({"gene": ["BRCA1"], "score": [0.7]})
    ann = AnnotationTable.from_pandas(df, provenance=_prov())
    out = tmp_path / "rows.parquet"
    core_io.save(ann, out)

    assert out.exists()
    assert (out.with_suffix(".parquet.provenance.json")).exists()

    loaded = core_io.load(out)
    assert isinstance(loaded, AnnotationTable)
    assert loaded.backend == "pandas"
    assert loaded.provenance == ann.provenance
    pd.testing.assert_frame_equal(
        loaded.to_pandas().reset_index(drop=True),
        df.reset_index(drop=True),
    )


def test_save_load_filter_chain_round_trip(tmp_path):
    """Save, load, filter on the loaded artifact — provenance survives."""
    df = pd.DataFrame({"gene": ["BRCA1", "BRCA2"], "score": [0.7, 0.4]})
    ann = AnnotationTable.from_pandas(df, provenance=_prov())
    out = tmp_path / "rows.parquet"
    core_io.save(ann, out)

    loaded = core_io.load(out)
    filtered = loaded.filter(col("score") > 0.5).collect()
    assert filtered == [{"gene": "BRCA1", "score": 0.7}]
    assert loaded.provenance == ann.provenance


def test_save_load_anndata_round_trip(tmp_path):
    obs = pd.DataFrame({"tissue": ["liver"]}, index=["s1"])
    var = pd.DataFrame({"gene": ["BRCA1"]}, index=["g1"])
    adata = ad.AnnData(X=np.array([[1.0]]), obs=obs, var=var)
    em = ExpressionMatrix.from_anndata(adata, provenance=_prov("expr-v1"))
    out = tmp_path / "expr.h5ad"
    core_io.save(em, out)

    loaded = core_io.load(out)
    assert isinstance(loaded, ExpressionMatrix)
    assert loaded.n_obs == 1
    assert loaded.n_vars == 1
    assert loaded.provenance == em.provenance
