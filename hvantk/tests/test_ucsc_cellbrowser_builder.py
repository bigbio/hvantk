"""Round-trip and update tests for the UCSC Cell Browser builder skill."""

from __future__ import annotations

from pathlib import Path

import pytest

from hvantk.tests._snapshot_utils import (
    anndata_sample_rows,
    anndata_schema_to_dict,
    load_snapshot,
)
from hvantk.tests._snapshot_utils import regenerate_snapshots as regenerate_snapshots_fn

EXPR_FIXTURE = "hvantk/tests/testdata/raw/ucsc-cellbrowser/expression_matrix.tsv"
META_FIXTURE = "hvantk/tests/testdata/raw/ucsc-cellbrowser/metadata.tsv"
SNAPSHOT_DIR = Path("hvantk/tests/snapshots/ucsc-cellbrowser")


def test_ucsc_cellbrowser_round_trip(tmp_path, regenerate_snapshots):
    """Build UCSC sc AnnData from fixture; assert schema and head sample stability."""
    from hvantk.tables.matrix_builders import build_ucsc_ad

    if regenerate_snapshots:
        regenerate_snapshots_fn(
            builder_fn=build_ucsc_ad,
            fixture_path=EXPR_FIXTURE,
            snapshot_dir=SNAPSHOT_DIR,
            builder_kwargs={
                "metadata_path": META_FIXTURE,
                "output_path": str(tmp_path / "ucsc.h5ad"),
                "overwrite": True,
            },
            input_path_kwarg="expression_matrix_path",
        )
        pytest.skip("Snapshots regenerated; rerun without --regenerate-snapshots to assert.")

    adata = build_ucsc_ad(
        expression_matrix_path=EXPR_FIXTURE,
        metadata_path=META_FIXTURE,
        output_path=str(tmp_path / "ucsc.h5ad"),
        overwrite=True,
    )

    expected_schema = load_snapshot(SNAPSHOT_DIR / "schema.json")
    actual_schema = anndata_schema_to_dict(adata)
    assert actual_schema == expected_schema, "UCSC sc schema drifted from snapshot"

    expected_rows = load_snapshot(SNAPSHOT_DIR / "sample_rows.json")
    actual_rows = anndata_sample_rows(adata)
    assert actual_rows == expected_rows, "UCSC sc sample rows drifted from snapshot"
