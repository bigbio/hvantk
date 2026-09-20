"""Snapshot round-trip test for the Expression Atlas builder.

Builds from the committed fixture and asserts the AnnData schema and a sample
of rows against committed snapshots.

Uses the AnnData wrapper pattern rather than ``phase_b_snapshot_adapter``:
that adapter ends with ``artifact.to_hail()``, which is wrong for an
AnnData-backed ExpressionMatrix. Returning the AnnData directly makes
``regenerate_snapshots`` take its anndata branch. See
hvantk/skills/cptac/expression/tests/test_builder.py and
hvantk/skills/ucsc_cellbrowser/tests/test_builder.py for the same shape.

Needs no Hail, so this runs in the default pytest selection.

Regenerate after an intentional change:
    pytest hvantk/skills/expression_atlas/tests/test_builder.py --regenerate-snapshots

Fixture derivation recipe (source: the real upstream files already committed
at hvantk/tests/testdata/raw/expression_atlas/, ~116k transcripts x 320
samples for the expression matrix -- far too large to commit as a fixture):

  1. Expression matrix (E-MTAB-6798-transcripts-tpms.tsv.bgz -> .tsv):
     read the bgzipped TSV. Its real header is
     "Gene ID\tGene Name\tGeneID\t<320 ERR* sample columns>" -- note the
     THIRD column "GeneID" is the per-row TRANSCRIPT id (e.g.
     "ENSMUST00000000001"), distinct from "Gene ID" (gene id) by a
     one-space difference. create_anndata_from_expression_atlas() only
     strips {"Gene ID", "Gene Name"} as non-sample columns
     (gene_column/gene_name_column defaults), so keeping "GeneID" verbatim
     makes the builder try to cast transcript-id strings into the float32
     expression matrix and crash:
         ValueError: could not convert string to float: 'ENSMUST...'
     Reproduced independently on a 1-row/2-sample toy frame. This is a
     latent builder bug present in the real upstream file shape, not an
     artifact of truncation; fixing it is out of scope here (builder.py /
     shared/expression_atlas.py are not in the allowed edit list for this
     fixture-seeding task -- see local report to the ledger owner).
     The fixture therefore DROPS the "GeneID" transcript column and keeps
     the first 20 DISTINCT genes (by "Gene ID", first transcript row seen,
     via drop_duplicates) so var_names come out unique, then keeps 4 sample
     columns (ERR2588382, ERR2588384, ERR2588383, ERR2588399) -- the first
     four in header order. Written back out as plain uncompressed TSV.
  2. SDRF (E-MTAB-6798.condensed-sdrf.tsv): long-format, one row per
     (sample, characteristic|factor). Kept every line whose 3rd
     tab-separated field (sample_id) is one of the 4 retained sample IDs
     above, preserving the original ragged field layout (some rows have a
     trailing ontology-URI column, some don't) so
     convert_sdrf_to_dataframe()'s usecols-inference keeps working
     unmodified.

  Script: /private/tmp/claude-502/-Users-enrique-projects-github-pyvatk/
  3a2337e9-30fa-4bab-8821-edd370a160ba/scratchpad/derive_fixture.py
  (not committed; rerun the recipe above against the real files to
  reproduce).
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

FIXTURE_DIR = Path("hvantk/skills/expression_atlas/tests/testdata/raw/expression-atlas")
EXPRESSION = str(FIXTURE_DIR / "E-MTAB-6798-transcripts-tpms.tsv")
SDRF = str(FIXTURE_DIR / "E-MTAB-6798.condensed-sdrf.tsv")
SNAPSHOT_DIR = Path("hvantk/skills/expression_atlas/tests/snapshots")

# The AnnData snapshot path samples the first N rows rather than selecting by
# key (anndata_sample_rows(adata, n=5)), so there is no key list to maintain
# here.


def _fake_ctx():
    """Deterministic BuildContext, so provenance never perturbs a snapshot."""
    from hvantk.core.models.build_context import BuildContext

    return BuildContext(
        plugin="expression-atlas",
        dataset="expression-atlas:dataset",
        plugin_version="test",
        source_fingerprint="sha256:test",
        builder_commit=None,
    )


def _build_for_snapshot(expression_matrix_path, **call_kwargs):
    """Adapt the Phase B builder to the snapshot helper's calling convention."""
    from hvantk.skills.expression_atlas.builder import build_expression_atlas

    sdrf_path = call_kwargs.pop("sdrf_path", SDRF)
    call_kwargs.pop("output_path", None)
    call_kwargs.pop("overwrite", None)

    artifact = build_expression_atlas(
        parsed_input={"expression_matrix": expression_matrix_path, "sdrf": sdrf_path},
        ctx=_fake_ctx(),
        **call_kwargs,
    )
    return artifact.to_anndata()


def test_expression_atlas_snapshot_round_trip(tmp_path, regenerate_snapshots):
    """Build Expression Atlas from fixture; assert schema and sample-row stability."""
    if regenerate_snapshots:
        regenerate_snapshots_fn(
            builder_fn=_build_for_snapshot,
            fixture_path=EXPRESSION,
            snapshot_dir=SNAPSHOT_DIR,
            builder_kwargs={"sdrf_path": SDRF},
            input_path_kwarg="expression_matrix_path",
        )
        pytest.skip("Snapshots regenerated; rerun without --regenerate-snapshots to assert.")

    adata = _build_for_snapshot(EXPRESSION, sdrf_path=SDRF)

    assert adata.n_obs == 4, "four samples in the fixture"
    assert adata.n_vars == 20, "twenty genes in the fixture"

    expected_schema = load_snapshot(SNAPSHOT_DIR / "schema.json")
    assert anndata_schema_to_dict(adata) == expected_schema, \
        "Expression Atlas schema drifted from snapshot"

    expected_rows = load_snapshot(SNAPSHOT_DIR / "sample_rows.json")
    assert anndata_sample_rows(adata) == expected_rows, \
        "Expression Atlas sample rows drifted from snapshot"
