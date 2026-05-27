"""Round-trip and update tests for the UCSC Cell Browser builder skill.

Parametrized over per-collection fixtures so adding a new UCSC collection
that ships a structurally distinct ``obs`` schema only requires:

  1. Dropping raw TSVs under
     ``hvantk/skills/ucsc_cellbrowser/tests/testdata/raw/<source>/``.
  2. Appending a case to ``_UCSC_CASES`` below.
  3. Running ``pytest <this file> --regenerate-snapshots`` once.

The builder, snapshot helpers, and plugin manifest entry stay shared.
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


_FIXTURE_ROOT = Path("hvantk/skills/ucsc_cellbrowser/tests/testdata/raw")
_SNAPSHOT_ROOT = Path("hvantk/skills/ucsc_cellbrowser/tests/snapshots")


# Per-collection cases. ``builder_kwargs`` carries plugin-derived overrides
# (e.g. ``gene_column``, ``split_gene_field``) when a dataset deviates from
# Phase B ``build_ucsc_cellbrowser`` defaults; the default values produce empty
# dicts. ``dataset_name`` selects the per-dataset schema_id in the builder's
# ``_SCHEMA_IDS`` table.
_UCSC_CASES = [
    pytest.param(
        str(_FIXTURE_ROOT / "ucsc-cellbrowser" / "expression_matrix.tsv"),
        str(_FIXTURE_ROOT / "ucsc-cellbrowser" / "metadata.tsv"),
        _SNAPSHOT_ROOT / "ucsc-cellbrowser",
        "ucsc-cellbrowser:default",
        {},
        id="asp_2019-celltype-summary",
    ),
    pytest.param(
        str(_FIXTURE_ROOT / "ucsc-cellbrowser-adult-ctx" / "expression_matrix.tsv"),
        str(_FIXTURE_ROOT / "ucsc-cellbrowser-adult-ctx" / "metadata.tsv"),
        _SNAPSHOT_ROOT / "ucsc-cellbrowser-adult-ctx",
        "ucsc-cellbrowser:adult-ctx",
        {},
        id="adult-ctx-meta-atlas-class-summary",
    ),
    pytest.param(
        str(_FIXTURE_ROOT / "ucsc-cellbrowser-dev-ctx" / "expression_matrix.tsv"),
        str(_FIXTURE_ROOT / "ucsc-cellbrowser-dev-ctx" / "metadata.tsv"),
        _SNAPSHOT_ROOT / "ucsc-cellbrowser-dev-ctx",
        "ucsc-cellbrowser:dev-ctx",
        {},
        id="dev-ctx-meta-atlas-type-v2-summary",
    ),
]


def _fake_ctx(dataset_name: str):
    """Construct a deterministic BuildContext for snapshot tests."""
    from hvantk.core.models.build_context import BuildContext

    plugin = dataset_name.split(":", 1)[0]
    return BuildContext(
        plugin=plugin,
        dataset=dataset_name,
        plugin_version="test",
        source_fingerprint="sha256:test",
        builder_commit=None,
    )


def _build_for_snapshot(expression_matrix_path, **call_kwargs):
    """Phase B build wrapper for snapshot regeneration / assertion.

    Accepts the (input_path, metadata_path, ...) signature the snapshot
    helper expects plus an injected ``dataset_name``; returns the underlying
    ``AnnData`` so the existing snapshot helpers apply without modification.
    """
    from hvantk.skills.ucsc_cellbrowser.builder import build_ucsc_cellbrowser

    dataset_name = call_kwargs.pop("dataset_name")
    metadata_path = call_kwargs.pop("metadata_path")
    # Strip args that only existed on the Phase A signature.
    call_kwargs.pop("output_path", None)
    call_kwargs.pop("overwrite", None)

    artifact = build_ucsc_cellbrowser(
        parsed_input={
            "expression_matrix": expression_matrix_path,
            "metadata": metadata_path,
        },
        ctx=_fake_ctx(dataset_name),
        **call_kwargs,
    )
    return artifact.to_anndata()


@pytest.mark.parametrize(
    ("expr_fixture", "meta_fixture", "snapshot_dir", "dataset_name", "builder_kwargs"),
    _UCSC_CASES,
)
def test_ucsc_cellbrowser_round_trip(
    tmp_path,
    regenerate_snapshots,
    expr_fixture,
    meta_fixture,
    snapshot_dir,
    dataset_name,
    builder_kwargs,
):
    """Build UCSC sc AnnData from fixture; assert schema and head sample stability."""
    call_kwargs = {
        "metadata_path": meta_fixture,
        "dataset_name": dataset_name,
        **builder_kwargs,
    }

    if regenerate_snapshots:
        regenerate_snapshots_fn(
            builder_fn=_build_for_snapshot,
            fixture_path=expr_fixture,
            snapshot_dir=snapshot_dir,
            builder_kwargs=call_kwargs,
            input_path_kwarg="expression_matrix_path",
        )
        pytest.skip("Snapshots regenerated; rerun without --regenerate-snapshots to assert.")

    adata = _build_for_snapshot(
        expression_matrix_path=expr_fixture,
        **call_kwargs,
    )

    expected_schema = load_snapshot(snapshot_dir / "schema.json")
    actual_schema = anndata_schema_to_dict(adata)
    assert actual_schema == expected_schema, f"UCSC sc schema drifted from snapshot ({snapshot_dir})"

    expected_rows = load_snapshot(snapshot_dir / "sample_rows.json")
    actual_rows = anndata_sample_rows(adata)
    assert actual_rows == expected_rows, f"UCSC sc sample rows drifted from snapshot ({snapshot_dir})"
