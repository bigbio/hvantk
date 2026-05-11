"""Round-trip and update tests for the UCSC Cell Browser builder skill.

Parametrized over per-collection fixtures so adding a new UCSC collection
that ships a structurally distinct ``obs`` schema only requires:

  1. Dropping raw TSVs under ``hvantk/tests/testdata/raw/<source>/``.
  2. Appending a case to ``_UCSC_CASES`` below.
  3. Running ``pytest <this file> --regenerate-snapshots`` once.

The builder, snapshot helpers, and registry entry stay shared.
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


# Per-collection cases. ``builder_kwargs`` carries registry-derived overrides
# (e.g. ``gene_column``, ``split_gene_field``) when a dataset deviates from
# ``build_ucsc_ad`` defaults; the default values produce empty dicts.
_UCSC_CASES = [
    pytest.param(
        "hvantk/tests/testdata/raw/ucsc-cellbrowser/expression_matrix.tsv",
        "hvantk/tests/testdata/raw/ucsc-cellbrowser/metadata.tsv",
        Path("hvantk/tests/snapshots/ucsc-cellbrowser"),
        {},
        id="asp_2019-celltype-summary",
    ),
    pytest.param(
        "hvantk/tests/testdata/raw/ucsc-cellbrowser-adult-ctx/expression_matrix.tsv",
        "hvantk/tests/testdata/raw/ucsc-cellbrowser-adult-ctx/metadata.tsv",
        Path("hvantk/tests/snapshots/ucsc-cellbrowser-adult-ctx"),
        {},
        id="adult-ctx-meta-atlas-class-summary",
    ),
    pytest.param(
        "hvantk/tests/testdata/raw/ucsc-cellbrowser-dev-ctx/expression_matrix.tsv",
        "hvantk/tests/testdata/raw/ucsc-cellbrowser-dev-ctx/metadata.tsv",
        Path("hvantk/tests/snapshots/ucsc-cellbrowser-dev-ctx"),
        {},
        id="dev-ctx-meta-atlas-type-v2-summary",
    ),
]


@pytest.mark.parametrize(
    ("expr_fixture", "meta_fixture", "snapshot_dir", "builder_kwargs"),
    _UCSC_CASES,
)
def test_ucsc_cellbrowser_round_trip(
    tmp_path,
    regenerate_snapshots,
    expr_fixture,
    meta_fixture,
    snapshot_dir,
    builder_kwargs,
):
    """Build UCSC sc AnnData from fixture; assert schema and head sample stability."""
    from hvantk.tables.matrix_builders import build_ucsc_ad

    call_kwargs = {
        "metadata_path": meta_fixture,
        "output_path": str(tmp_path / "ucsc.h5ad"),
        "overwrite": True,
        **builder_kwargs,
    }

    if regenerate_snapshots:
        regenerate_snapshots_fn(
            builder_fn=build_ucsc_ad,
            fixture_path=expr_fixture,
            snapshot_dir=snapshot_dir,
            builder_kwargs=call_kwargs,
            input_path_kwarg="expression_matrix_path",
        )
        pytest.skip("Snapshots regenerated; rerun without --regenerate-snapshots to assert.")

    adata = build_ucsc_ad(
        expression_matrix_path=expr_fixture,
        **call_kwargs,
    )

    expected_schema = load_snapshot(snapshot_dir / "schema.json")
    actual_schema = anndata_schema_to_dict(adata)
    assert actual_schema == expected_schema, f"UCSC sc schema drifted from snapshot ({snapshot_dir})"

    expected_rows = load_snapshot(snapshot_dir / "sample_rows.json")
    actual_rows = anndata_sample_rows(adata)
    assert actual_rows == expected_rows, f"UCSC sc sample rows drifted from snapshot ({snapshot_dir})"
