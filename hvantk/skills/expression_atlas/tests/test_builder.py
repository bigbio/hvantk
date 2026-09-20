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
     one-space difference.

     That third column is KEPT, deliberately: the fixture reproduces the real
     production header rather than a sanitised one. It used to be dropped,
     because create_anndata_from_expression_atlas() classified it as a sample
     and died casting transcript strings to float32 (issue #342) -- so the
     passing test was exercising a file shape no download ever produces. The
     builder now sets non-numeric leading columns aside into .var, and this
     fixture is what proves it.

     Kept the first 20 DISTINCT genes (by "Gene ID", first transcript row
     seen) so var_names come out unique, then 4 sample columns
     (ERR2588382, ERR2588384, ERR2588383, ERR2588399) -- the first four in
     header order. Written back out as plain uncompressed TSV.
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
        pytest.skip(
            "Snapshots regenerated; rerun without --regenerate-snapshots to assert."
        )

    adata = _build_for_snapshot(EXPRESSION, sdrf_path=SDRF)

    assert adata.n_obs == 4, "four samples in the fixture"
    assert adata.n_vars == 20, "twenty genes in the fixture"

    expected_schema = load_snapshot(SNAPSHOT_DIR / "schema.json")
    assert (
        anndata_schema_to_dict(adata) == expected_schema
    ), "Expression Atlas schema drifted from snapshot"

    expected_rows = load_snapshot(SNAPSHOT_DIR / "sample_rows.json")
    assert (
        anndata_sample_rows(adata) == expected_rows
    ), "Expression Atlas sample rows drifted from snapshot"


# --- Regression: issue #342 ---------------------------------------------------
# The real upstream header carries THREE leading metadata columns, the third being
# `GeneID` (transcript id, distinct from `Gene ID` by one space). Classifying it as a
# sample sent transcript strings into the float32 cast, so every build from an
# unmodified download died. The committed fixture now reproduces that header, but these
# pin the behaviour directly -- a fixture can be re-derived, an assertion cannot drift.


def test_transcript_id_column_is_annotation_not_sample(tmp_path):
    """The third metadata column must reach .var, never the expression matrix."""
    import pandas as pd
    from hvantk.skills.expression_atlas.shared.expression_atlas import (
        create_anndata_from_expression_atlas,
    )

    path = tmp_path / "tpms.tsv"
    pd.DataFrame(
        {
            "Gene ID": ["ENSMUSG00000000001", "ENSMUSG00000000002"],
            "Gene Name": ["Gnai3", "Cdc45"],
            "GeneID": ["ENSMUST00000000001", "ENSMUST00000000002"],
            "ERR1": [16, 3],
            "ERR2": [7, 1],
        }
    ).to_csv(path, sep="\t", index=False)

    adata = create_anndata_from_expression_atlas(expression_matrix_path=str(path))

    assert adata.shape == (
        2,
        2,
    ), "2 samples x 2 genes -- GeneID must not become a third sample"
    assert list(adata.obs_names) == ["ERR1", "ERR2"]
    # Preserved, not dropped: it is the only thing disambiguating repeated gene ids.
    assert "GeneID" in adata.var.columns
    assert list(adata.var["GeneID"]) == ["ENSMUST00000000001", "ENSMUST00000000002"]


def test_all_missing_sample_column_stays_a_sample(tmp_path):
    """Inference must key on 'non-numeric', not 'not obviously a number'.

    A sample with no measurements reads back as all-NaN. Treating that as metadata
    would silently drop a real sample from the matrix -- the failure mode that makes
    content-based inference risky, so it is pinned here.
    """
    import numpy as np
    import pandas as pd
    from hvantk.skills.expression_atlas.shared.expression_atlas import (
        create_anndata_from_expression_atlas,
    )

    path = tmp_path / "tpms.tsv"
    pd.DataFrame(
        {
            "Gene ID": ["ENSMUSG00000000001"],
            "Gene Name": ["Gnai3"],
            "ERR1": [16.0],
            "ERR_empty": [np.nan],
        }
    ).to_csv(path, sep="\t", index=False)

    adata = create_anndata_from_expression_atlas(expression_matrix_path=str(path))
    assert list(adata.obs_names) == ["ERR1", "ERR_empty"]


def test_malformed_sample_named_in_sdrf_raises_rather_than_vanishing(tmp_path):
    """A declared sample carrying a bad sentinel must fail loudly, not move to .var.

    Inferring annotations from content alone made this silent: an all-text sample was
    reclassified as metadata, dropped from the matrix, and the build SUCCEEDED with
    that sample missing -- strictly worse than the crash the inference was added to
    prevent. `metadata_df.index` is the authoritative sample list, so it outranks the
    heuristic.
    """
    import pandas as pd
    import pytest
    from hvantk.skills.expression_atlas.shared.expression_atlas import (
        create_anndata_from_expression_atlas,
    )

    path = tmp_path / "tpms.tsv"
    pd.DataFrame(
        {
            "Gene ID": ["ENSMUSG00000000001", "ENSMUSG00000000002"],
            "Gene Name": ["Gnai3", "Cdc45"],
            "GeneID": ["ENSMUST00000000001", "ENSMUST00000000002"],
            "ERR1": [16, 3],
            # "-" is NOT one of pandas' recognised NA strings, so it survives read_csv
            # as text and would otherwise look exactly like an annotation column.
            "ERR_broken": ["-", "-"],
        }
    ).to_csv(path, sep="\t", index=False)
    metadata = pd.DataFrame(
        index=["ERR1", "ERR_broken"], data={"organism": ["Mus"] * 2}
    )

    with pytest.raises(ValueError, match="could not convert string to float"):
        create_anndata_from_expression_atlas(
            expression_matrix_path=str(path), metadata_df=metadata
        )


def test_malformed_sample_without_sdrf_still_raises(tmp_path):
    """The positional guard covers the case where no sample list is supplied.

    Annotations are a LEADING block in this format, so a non-numeric column appearing
    after the first sample is malformed data, not metadata -- regardless of whether a
    caller passed metadata_df.
    """
    import pandas as pd
    import pytest
    from hvantk.skills.expression_atlas.shared.expression_atlas import (
        create_anndata_from_expression_atlas,
    )

    path = tmp_path / "tpms.tsv"
    pd.DataFrame(
        {
            "Gene ID": ["ENSMUSG00000000001"],
            "Gene Name": ["Gnai3"],
            "GeneID": ["ENSMUST00000000001"],
            "ERR1": [16],
            "ERR_broken": ["-"],
        }
    ).to_csv(path, sep="\t", index=False)

    with pytest.raises(ValueError, match="could not convert string to float"):
        create_anndata_from_expression_atlas(expression_matrix_path=str(path))
