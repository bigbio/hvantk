"""#198: `hvantk download cptac-phospho --all` must skip-and-continue per cancer.

coad + ov fail *inside* cptac 1.5.14 (upstream); one bad type must not abort the
whole batch. These tests stub CPTACPhosphoDataset so no `cptac` package or network
is needed.
"""
from __future__ import annotations

import pytest


def test_download_dataset_skips_failing_cancer(monkeypatch, tmp_path):
    import hvantk.skills.cptac.phospho.cli as cli

    class _FakeDS:
        def __init__(self, cancer_type):
            self.ct = cancer_type

        def download(self, raw_dir, overwrite=False):
            if self.ct == "coad":
                raise RuntimeError("upstream cptac boom")
            return {
                "tsv": f"{raw_dir}/{self.ct}.tsv",
                "matrix": f"{raw_dir}/{self.ct}-matrix.csv",
                "metadata": f"{raw_dir}/{self.ct}-metadata.csv",
            }

    monkeypatch.setattr(
        "hvantk.skills.cptac.shared.datasets.CPTACPhosphoDataset", _FakeDS
    )

    res = cli.download_dataset(str(tmp_path), cancer_type=None)  # all cancers

    assert "coad" in res["_failures"]          # failure recorded, not raised
    assert "brca" in res                        # others still succeeded
    assert res["brca"]["tsv"].endswith("brca.tsv")


def test_download_dataset_all_fail_raises(monkeypatch, tmp_path):
    import hvantk.skills.cptac.phospho.cli as cli

    class _BoomDS:
        def __init__(self, cancer_type):
            pass

        def download(self, *a, **k):
            raise RuntimeError("boom")

    monkeypatch.setattr(
        "hvantk.skills.cptac.shared.datasets.CPTACPhosphoDataset", _BoomDS
    )

    with pytest.raises(RuntimeError):
        cli.download_dataset(str(tmp_path), cancer_type="brca")
