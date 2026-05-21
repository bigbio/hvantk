"""Phase P: artifact-typed wrapper for compute_specificity."""
from __future__ import annotations

from datetime import datetime, timezone

import pandas as pd

from hvantk.algorithms.expression.tissue_specificity import (
    compute_specificity_artifact,
)
from hvantk.core.models import AnnotationTable, Provenance


def _prov() -> Provenance:
    return Provenance(
        plugin="t", dataset="t:expr", plugin_version="0",
        source_fingerprint="sha256:x", schema_id="t-expr-v1",
        build_timestamp=datetime(2026, 5, 21, tzinfo=timezone.utc),
        builder_commit=None,
    )


def test_specificity_artifact_returns_annotation_table():
    df = pd.DataFrame({
        "gene_id": ["G1", "G2", "G3"],
        "liver": [10.0, 5.0, 1.0],
        "brain": [1.0, 5.0, 10.0],
        "heart": [1.0, 5.0, 1.0],
    })
    ann = AnnotationTable.from_pandas(df, provenance=_prov())

    result = compute_specificity_artifact(ann)
    assert isinstance(result, AnnotationTable)
    out_df = result.to_pandas()
    assert set(out_df.columns) == {"gene_id", "specificity"}
    assert len(out_df) == 3


def test_specificity_artifact_chains_provenance():
    """The @algorithm decorator should chain input provenance to output."""
    df = pd.DataFrame({
        "gene_id": ["G1", "G2"],
        "liver": [10.0, 5.0],
        "brain": [1.0, 5.0],
    })
    src_prov = _prov()
    ann = AnnotationTable.from_pandas(df, provenance=src_prov)

    result = compute_specificity_artifact(ann)
    # Provenance chaining: the input prov becomes a parent of the output prov.
    assert src_prov in result.provenance.parents
