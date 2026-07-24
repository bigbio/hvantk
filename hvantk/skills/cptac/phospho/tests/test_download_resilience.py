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

    # Assert the NEW wrapped message (not just any RuntimeError — the pre-fix code
    # also propagated a RuntimeError, so an unmatched raises would pass either way).
    with pytest.raises(RuntimeError, match="All CPTAC phospho downloads failed"):
        cli.download_dataset(str(tmp_path), cancer_type="brca")


def test_download_cmd_all_skips_failing_cancer(monkeypatch, tmp_path):
    """The Click `download cptac-phospho --all` path: one bad cancer is skipped,
    survivors are merged, a summary prints, and the command exits 0."""
    import os

    from click.testing import CliRunner

    import hvantk.skills.cptac.phospho.cli as cli
    from hvantk.skills.cptac.shared.datasets import _TSV_COLUMNS

    class _FakeDS:
        def __init__(self, cancer_type):
            self.ct = cancer_type

        def download(self, output_dir, overwrite=False):
            if self.ct in ("coad", "ov"):
                raise RuntimeError(f"upstream cptac boom for {self.ct}")
            tsv = os.path.join(output_dir, f"cptac-phospho-{self.ct}.tsv")
            with open(tsv, "w") as fh:
                fh.write("\t".join(_TSV_COLUMNS) + "\n")  # header, no rows
            return {
                "tsv": tsv,
                "matrix": os.path.join(output_dir, f"cptac-phospho-{self.ct}-matrix.csv"),
                "metadata": os.path.join(output_dir, f"cptac-phospho-{self.ct}-metadata.csv"),
            }

    monkeypatch.setattr(
        "hvantk.skills.cptac.shared.datasets.CPTACPhosphoDataset", _FakeDS
    )

    result = CliRunner().invoke(cli.download_cmd, ["-o", str(tmp_path), "--all"])

    assert result.exit_code == 0, result.output
    assert "coad" in result.output and "FAILED" in result.output
    assert "Summary" in result.output
    # Pan-cancer merge ran over the survivors only.
    assert (tmp_path / "cptac-phospho-pancancer.tsv").exists()


def test_download_cmd_all_fail_exits_nonzero(monkeypatch, tmp_path):
    """When every cancer fails, `download --all` exits non-zero."""
    from click.testing import CliRunner

    import hvantk.skills.cptac.phospho.cli as cli

    class _BoomDS:
        def __init__(self, cancer_type):
            pass

        def download(self, *a, **k):
            raise RuntimeError("boom")

    monkeypatch.setattr(
        "hvantk.skills.cptac.shared.datasets.CPTACPhosphoDataset", _BoomDS
    )

    result = CliRunner().invoke(cli.download_cmd, ["-o", str(tmp_path), "--all"])
    assert result.exit_code == 1
