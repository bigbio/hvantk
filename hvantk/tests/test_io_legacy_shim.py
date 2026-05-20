"""Legacy raw-file shim: core/io.load wraps unmanifested files with Provenance.unknown."""
from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from hvantk.core import io as core_io
from hvantk.core.models.annotation_table import AnnotationTable
from hvantk.core.models.expression_matrix import ExpressionMatrix


def test_legacy_parquet_loaded_with_unknown_provenance(tmp_path):
    df = pd.DataFrame({"x": [1, 2, 3]})
    legacy_path = tmp_path / "legacy.parquet"
    df.to_parquet(legacy_path, index=False)

    ann = core_io.load(legacy_path)
    assert isinstance(ann, AnnotationTable)
    assert ann.provenance.plugin == "<unknown>"
    assert "legacy file" in ann.provenance.source_fingerprint


def test_legacy_h5ad_loaded_with_unknown_provenance(tmp_path):
    adata = ad.AnnData(X=np.array([[1.0]]),
                       obs=pd.DataFrame({"t": ["a"]}, index=["s1"]),
                       var=pd.DataFrame({"g": ["g1"]}, index=["v1"]))
    legacy_path = tmp_path / "legacy.h5ad"
    adata.write_h5ad(str(legacy_path))

    em = core_io.load(legacy_path)
    assert isinstance(em, ExpressionMatrix)
    assert em.provenance.plugin == "<unknown>"


def test_legacy_loaded_artifact_is_usable(tmp_path):
    """Legacy-loaded artifacts must still support the portable API."""
    from hvantk.core.models._expr import col

    df = pd.DataFrame({"x": [1, 2, 3]})
    legacy_path = tmp_path / "legacy.parquet"
    df.to_parquet(legacy_path, index=False)

    ann = core_io.load(legacy_path)
    out = ann.filter(col("x") > 1).collect()
    assert out == [{"x": 2}, {"x": 3}]
