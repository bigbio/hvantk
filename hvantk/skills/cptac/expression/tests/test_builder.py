"""Snapshot round-trip test for the CPTAC expression builder.

Builds from the committed fixture and asserts the AnnData schema and a sample of rows
against committed snapshots.

Uses the AnnData wrapper pattern rather than ``phase_b_snapshot_adapter``: that adapter
ends with ``artifact.to_hail()``, which is wrong for an AnnData-backed ExpressionMatrix.
Returning the AnnData directly makes ``regenerate_snapshots`` take its anndata branch.
See hvantk/skills/ucsc_cellbrowser/tests/test_builder.py for the same shape.

Needs no Hail, so this runs in the default pytest selection.

Regenerate after an intentional change:
    pytest hvantk/skills/cptac/expression/tests/test_builder.py --regenerate-snapshots
"""

from __future__ import annotations

from pathlib import Path

import pytest

from hvantk.tests._snapshot_utils import (
    anndata_sample_rows,
    anndata_schema_to_dict,
    load_snapshot,
)
from hvantk.tests._snapshot_utils import regenerate_snapshots as regenerate_snapshots_fn

FIXTURE_DIR = Path("hvantk/skills/cptac/expression/tests/testdata/raw/cptac-expression")
EXPRESSION = str(FIXTURE_DIR / "expression.tsv")
METADATA = str(FIXTURE_DIR / "metadata.tsv")
SNAPSHOT_DIR = Path("hvantk/skills/cptac/expression/tests/snapshots")

# The AnnData snapshot path samples the first N rows rather than selecting by key
# (anndata_sample_rows(adata, n=5)), so there is no key list to maintain here.


def _fake_ctx():
    """Deterministic BuildContext, so provenance never perturbs a snapshot."""
    from hvantk.core.models.build_context import BuildContext

    return BuildContext(
        plugin="cptac",
        dataset="cptac:expression",
        plugin_version="test",
        source_fingerprint="sha256:test",
        builder_commit=None,
    )


def _build_for_snapshot(expression_path, **call_kwargs):
    """Adapt the Phase B builder to the snapshot helper's calling convention."""
    from hvantk.skills.cptac.expression.builder import build_cptac_expression

    metadata_path = call_kwargs.pop("metadata_path", METADATA)
    call_kwargs.pop("output_path", None)
    call_kwargs.pop("overwrite", None)

    artifact = build_cptac_expression(
        parsed_input={"expression": expression_path, "metadata": metadata_path},
        ctx=_fake_ctx(),
        **call_kwargs,
    )
    return artifact.to_anndata()


def test_cptac_expression_snapshot_round_trip(tmp_path, regenerate_snapshots):
    """Build CPTAC expression from fixture; assert schema and sample-row stability."""
    if regenerate_snapshots:
        regenerate_snapshots_fn(
            builder_fn=_build_for_snapshot,
            fixture_path=EXPRESSION,
            snapshot_dir=SNAPSHOT_DIR,
            builder_kwargs={"metadata_path": METADATA},
            input_path_kwarg="expression_path",
        )
        pytest.skip("Snapshots regenerated; rerun without --regenerate-snapshots to assert.")

    adata = _build_for_snapshot(EXPRESSION, metadata_path=METADATA)

    assert adata.n_obs == 2, "two samples in the fixture"
    assert adata.n_vars == 2, "two genes in the fixture"

    expected_schema = load_snapshot(SNAPSHOT_DIR / "schema.json")
    assert anndata_schema_to_dict(adata) == expected_schema, \
        "CPTAC expression schema drifted from snapshot"

    expected_rows = load_snapshot(SNAPSHOT_DIR / "sample_rows.json")
    assert anndata_sample_rows(adata) == expected_rows, \
        "CPTAC expression sample rows drifted from snapshot"
