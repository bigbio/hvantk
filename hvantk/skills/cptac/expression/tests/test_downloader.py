"""Unit tests for ``CPTACExpressionDataset.download()`` (#386, PR #451 review).

Covers what ``test_builder.py::test_build_from_raw_dir`` does not: that path
starts from an already-long CSV someone else wrote, so nothing exercised
``download()`` itself or ``melt_proteomics_wide_to_long``'s wide-to-long
transform. It also covers the atomic-publish guarantee -- a clinical-metadata
fetch failure must leave an existing matrix+metadata pair completely
untouched, with no leftover temp file, rather than publishing a fresh matrix
beside stale (or missing) metadata.

No ``cptac`` package or network is needed: ``_load_cptac_dataset`` is
monkeypatched to return a small fake dataset object, mirroring the live-
verified ``get_proteomics(source=..., tissue_type=...)`` / ``get_clinical
(source=...)`` call shape documented in ``expression/SKILL.md`` s 4.
"""

from __future__ import annotations

import os
import stat

import pandas as pd
import pytest

import hvantk.skills.cptac.shared.datasets as shared_datasets
from hvantk.core.plugin.drift_runner import _current_umask

# Wide proteomics frame get_proteomics() would return: 2 samples x 3 genes, with
# one blank Database_ID (MYC) so melt_proteomics_wide_to_long's fall-back-to-Name
# path is exercised, not just the common Database_ID-present case.
_PROTEOMICS_COLUMNS = pd.MultiIndex.from_tuples(
    [("TP53", "ENSP00000269305.4"), ("EGFR", "ENSP00000275493.2"), ("MYC", "")],
    names=["Name", "Database_ID"],
)


def _proteomics_df():
    return pd.DataFrame(
        [[10.5, 5.2, 1.0], [20.3, 8.1, 2.0]],
        index=pd.Index(["S1", "S2.N"], name="Patient_ID"),
        columns=_PROTEOMICS_COLUMNS,
    )


def _clinical_df():
    return pd.DataFrame(
        {"stage": ["I", "II"]}, index=pd.Index(["S1", "S2.N"], name="Patient_ID")
    )


class _FakeCptacDataset:
    """Stand-in for a ``cptac.Brca()``-like object."""

    def __init__(self, fail_clinical: bool = False):
        self.fail_clinical = fail_clinical

    def get_proteomics(self, source=None, tissue_type="both"):
        assert source == "umich" and tissue_type == "both"
        return _proteomics_df()

    def get_clinical(self, source=None):
        if self.fail_clinical:
            raise RuntimeError("Zenodo 504")
        assert source == "mssm"
        return _clinical_df()


def _temp_files(raw_dir):
    """Leftover files matching the ``.<name>.<random>`` mkstemp convention."""
    return [n for n in os.listdir(raw_dir) if n.startswith(".")]


def _ctx():
    from hvantk.core.models import BuildContext

    return BuildContext(
        plugin="cptac",
        dataset="cptac:expression",
        plugin_version="test",
        source_fingerprint="sha256:test",
        builder_commit=None,
    )


@pytest.mark.parametrize("case", ["ok", "clinical_fails"])
def test_download_writes_atomically(case, tmp_path, monkeypatch):
    raw = tmp_path
    matrix_path = raw / "cptac-expression-brca-matrix.csv"
    metadata_path = raw / "cptac-expression-brca-metadata.csv"

    if case == "clinical_fails":
        # Pre-seed an old (sentinel) pair, as if a previous run had completed.
        matrix_path.write_text("SENTINEL-MATRIX\n")
        metadata_path.write_text("SENTINEL-METADATA\n")
        old_matrix = matrix_path.read_bytes()
        old_metadata = metadata_path.read_bytes()

        monkeypatch.setattr(
            shared_datasets,
            "_load_cptac_dataset",
            lambda ct: _FakeCptacDataset(fail_clinical=True),
        )

        with pytest.raises(RuntimeError, match="clinical"):
            shared_datasets.CPTACExpressionDataset(cancer_type="brca").download(
                str(raw), overwrite=True
            )

        assert matrix_path.read_bytes() == old_matrix, "old matrix must be untouched"
        assert metadata_path.read_bytes() == old_metadata, (
            "old metadata must be untouched"
        )
        assert not _temp_files(raw), "no temp file should be left behind"
        return

    # case == "ok"
    monkeypatch.setattr(
        shared_datasets, "_load_cptac_dataset", lambda ct: _FakeCptacDataset()
    )

    paths = shared_datasets.CPTACExpressionDataset(cancer_type="brca").download(
        str(raw)
    )
    assert paths == {"matrix": str(matrix_path), "metadata": str(metadata_path)}
    assert not _temp_files(raw), "no temp file should be left behind"

    # mkstemp always creates mode 0600; the published files must not keep that --
    # they should land at whatever a normal file create gets (0666 masked by the
    # process umask), so a shared raw dir stays group-readable.
    expected_mode = 0o666 & ~_current_umask()
    for p in (matrix_path, metadata_path):
        assert stat.S_IMODE(os.stat(p).st_mode) == expected_mode

    matrix_df = pd.read_csv(matrix_path)
    assert list(matrix_df.columns) == [
        "SampleID",
        "GeneID",
        "Gene Name",
        "Expression",
    ]
    assert len(matrix_df) == 6  # 2 samples x 3 genes
    assert set(matrix_df["GeneID"]) == {
        "ENSP00000269305.4",
        "ENSP00000275493.2",
        "MYC",  # blank Database_ID falls back to the gene symbol
    }

    metadata_df = pd.read_csv(metadata_path)
    assert "tissue_type" in metadata_df.columns

    from hvantk.core.models import ExpressionMatrix
    from hvantk.skills.cptac.expression.builder import build_cptac_expression

    artifact = build_cptac_expression(str(raw), _ctx(), cancer_type="brca")
    assert isinstance(artifact, ExpressionMatrix)
    assert artifact.to_anndata().shape == (2, 3)  # 2 samples x 3 genes
