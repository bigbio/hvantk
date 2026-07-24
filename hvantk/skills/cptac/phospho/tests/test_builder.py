"""Unit tests for the CPTAC phospho Phase B builder (#198).

These exercise the builder without the upstream ``cptac`` package: they feed
matrix/metadata CSVs in the same layout ``CPTACPhosphoDataset.download`` writes
(matrix index name ``Site``; metadata index name ``SampleID``).
"""
from __future__ import annotations

import pandas as pd
import pytest


def _ctx():
    from hvantk.core.models import BuildContext

    return BuildContext(
        plugin="cptac",
        dataset="cptac:phospho",
        plugin_version="0.1.0",
        source_fingerprint="test",
        builder_commit=None,
    )


def _write_fixture(raw, ct="brca"):
    pd.DataFrame(
        {"S1": [1.0, 2.0], "S2": [3.0, 4.0]},
        index=pd.Index(["TP53_S15", "EGFR_Y1068"], name="Site"),
    ).to_csv(raw / f"cptac-phospho-{ct}-matrix.csv")
    pd.DataFrame(
        {"stage": ["I", "II"]},
        index=pd.Index(["S1", "S2"], name="SampleID"),
    ).to_csv(raw / f"cptac-phospho-{ct}-metadata.csv")


def test_build_from_raw_dir(tmp_path):
    from hvantk.core.models import ExpressionMatrix
    from hvantk.skills.cptac.phospho.builder import build_cptac_phospho

    raw = tmp_path / "raw"
    raw.mkdir()
    _write_fixture(raw, "brca")

    art = build_cptac_phospho(str(raw), _ctx(), cancer_type="brca")

    assert isinstance(art, ExpressionMatrix)
    adata = art.to_anndata()
    assert adata.shape == (2, 2)  # 2 samples x 2 sites
    assert set(adata.obs.index) == {"S1", "S2"}
    assert set(adata.var.index) == {"TP53_S15", "EGFR_Y1068"}
    assert "stage" in adata.obs.columns


def test_build_from_dict(tmp_path):
    from hvantk.core.models import ExpressionMatrix
    from hvantk.skills.cptac.phospho.builder import build_cptac_phospho

    raw = tmp_path / "raw"
    raw.mkdir()
    _write_fixture(raw, "brca")

    art = build_cptac_phospho(
        {
            "expression": str(raw / "cptac-phospho-brca-matrix.csv"),
            "metadata": str(raw / "cptac-phospho-brca-metadata.csv"),
        },
        _ctx(),
    )
    assert isinstance(art, ExpressionMatrix)
    assert art.to_anndata().shape == (2, 2)


def test_raw_dir_requires_cancer_type(tmp_path):
    from hvantk.skills.cptac.phospho.builder import build_cptac_phospho

    raw = tmp_path / "raw"
    raw.mkdir()
    _write_fixture(raw, "brca")

    with pytest.raises(ValueError, match="cancer_type"):
        build_cptac_phospho(str(raw), _ctx())


def test_missing_files_raise(tmp_path):
    from hvantk.skills.cptac.phospho.builder import build_cptac_phospho

    raw = tmp_path / "raw"
    raw.mkdir()
    with pytest.raises(FileNotFoundError):
        build_cptac_phospho(str(raw), _ctx(), cancer_type="brca")


def test_ignores_stray_kwargs(tmp_path):
    from hvantk.skills.cptac.phospho.builder import build_cptac_phospho

    raw = tmp_path / "raw"
    raw.mkdir()
    _write_fixture(raw, "brca")

    # reprocess forwards download-only --plugin-arg values (e.g. overwrite) to the
    # builder too; they must not break the build.
    art = build_cptac_phospho(str(raw), _ctx(), cancer_type="brca", overwrite=True)
    assert art.to_anndata().shape == (2, 2)


# ---------------------------------------------------------------------------
# Snapshot round-trip
# ---------------------------------------------------------------------------
# Appended rather than put in its own file so the round-trip test lives in
# test_builder.py, per hvantk/skills/_conventions/SKILL.md section 8.
#
# Uses the AnnData wrapper pattern rather than phase_b_snapshot_adapter: that adapter ends
# with artifact.to_hail(), which is wrong for an AnnData-backed ExpressionMatrix. Returning
# the AnnData makes regenerate_snapshots take its anndata branch. Needs no Hail.
#
# Regenerate after an intentional change:
#     pytest hvantk/skills/cptac/phospho/tests/test_builder.py --regenerate-snapshots

from pathlib import Path  # noqa: E402

from hvantk.tests._snapshot_utils import (  # noqa: E402
    anndata_sample_rows,
    anndata_schema_to_dict,
    load_snapshot,
)
from hvantk.tests._snapshot_utils import (  # noqa: E402
    regenerate_snapshots as regenerate_snapshots_fn,
)

_FIXTURE_DIR = Path("hvantk/skills/cptac/phospho/tests/testdata/raw/cptac-phospho")
_PHOSPHO = str(_FIXTURE_DIR / "phospho.tsv")
_METADATA = str(_FIXTURE_DIR / "metadata.tsv")
_SNAPSHOT_DIR = Path("hvantk/skills/cptac/phospho/tests/snapshots")


def _build_for_snapshot(expression_path, **call_kwargs):
    """Adapt the Phase B builder to the snapshot helper's calling convention."""
    from hvantk.skills.cptac.phospho.builder import build_cptac_phospho

    metadata_path = call_kwargs.pop("metadata_path", _METADATA)
    call_kwargs.pop("output_path", None)
    call_kwargs.pop("overwrite", None)

    artifact = build_cptac_phospho(
        parsed_input={"expression": expression_path, "metadata": metadata_path},
        ctx=_ctx(),
        **call_kwargs,
    )
    return artifact.to_anndata()


def test_cptac_phospho_snapshot_round_trip(regenerate_snapshots):
    """Build CPTAC phospho from the committed fixture; assert schema/row stability."""
    if regenerate_snapshots:
        regenerate_snapshots_fn(
            builder_fn=_build_for_snapshot,
            fixture_path=_PHOSPHO,
            snapshot_dir=_SNAPSHOT_DIR,
            builder_kwargs={"metadata_path": _METADATA},
            input_path_kwarg="expression_path",
        )
        pytest.skip("Snapshots regenerated; rerun without --regenerate-snapshots to assert.")

    adata = _build_for_snapshot(_PHOSPHO, metadata_path=_METADATA)

    assert adata.n_obs == 2, "two samples in the fixture"
    assert adata.n_vars == 2, "two phosphosites in the fixture"

    expected_schema = load_snapshot(_SNAPSHOT_DIR / "schema.json")
    assert anndata_schema_to_dict(adata) == expected_schema, \
        "CPTAC phospho schema drifted from snapshot"

    expected_rows = load_snapshot(_SNAPSHOT_DIR / "sample_rows.json")
    assert anndata_sample_rows(adata) == expected_rows, \
        "CPTAC phospho sample rows drifted from snapshot"
