"""Tests for the gnomad-metrics downloader (URL map, lifecycle, CLI)."""

from __future__ import annotations

import urllib.request
from pathlib import Path

import pytest

from hvantk.skills.gnomad_metrics import cli as gnomad_cli
from hvantk.skills.gnomad_metrics.shared import constants as C

# --------------------------------------------------------------------------- #
# URL map (pure, no network)
# --------------------------------------------------------------------------- #


def test_default_is_v211_by_gene():
    assert C.DEFAULT_VERSION == "v2.1.1"
    assert C.DEFAULT_TABLE["v2.1.1"] == "by_gene"
    url = C.constraint_url(C.DEFAULT_VERSION)
    assert url.endswith("2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz")
    assert url.startswith(
        "https://storage.googleapis.com/gcp-public-data--gnomad/release/"
    )


def test_v211_by_transcript_and_v40_urls():
    assert C.constraint_url("v2.1.1", "by_transcript").endswith(
        "gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz"
    )
    v40 = C.constraint_url("v4.0")
    assert v40.endswith("v4.0/constraint/gnomad.v4.0.constraint_metrics.tsv")
    assert C.constraint_filename("v4.0") == "gnomad.v4.0.constraint_metrics.tsv"


@pytest.mark.parametrize(
    "version,table",
    [("v9.9", None), ("v2.1.1", "nope"), ("v4.0", "by_gene")],
)
def test_unknown_version_or_table_raises(version, table):
    with pytest.raises(ValueError):
        C.constraint_url(version, table)


# --------------------------------------------------------------------------- #
# Lifecycle download_dataset (mocked network)
# --------------------------------------------------------------------------- #


def _fake_urlretrieve(recorder):
    def _inner(url, filename):
        Path(filename).write_bytes(b"fake gnomad constraint payload")
        recorder.append((url, str(filename)))
        return filename, None

    return _inner


def test_download_dataset_default_writes_by_gene(tmp_path, monkeypatch):
    calls: list = []
    monkeypatch.setattr(urllib.request, "urlretrieve", _fake_urlretrieve(calls))

    out = gnomad_cli.download_dataset(str(tmp_path))

    assert Path(out).name == "gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"
    assert Path(out).exists()
    assert len(calls) == 1
    requested_url, _ = calls[0]
    assert requested_url == C.constraint_url("v2.1.1", "by_gene")


def test_download_dataset_v40(tmp_path, monkeypatch):
    calls: list = []
    monkeypatch.setattr(urllib.request, "urlretrieve", _fake_urlretrieve(calls))

    out = gnomad_cli.download_dataset(str(tmp_path), version="v4.0")

    assert Path(out).name == "gnomad.v4.0.constraint_metrics.tsv"
    assert calls[0][0] == C.constraint_url("v4.0")


def test_download_dataset_no_overwrite(tmp_path, monkeypatch):
    calls: list = []
    monkeypatch.setattr(urllib.request, "urlretrieve", _fake_urlretrieve(calls))
    gnomad_cli.download_dataset(str(tmp_path))  # first write
    with pytest.raises(FileExistsError):
        gnomad_cli.download_dataset(str(tmp_path))  # second without overwrite
    # overwrite=True succeeds
    gnomad_cli.download_dataset(str(tmp_path), overwrite=True)


# --------------------------------------------------------------------------- #
# Click command (mocked network)
# --------------------------------------------------------------------------- #


def test_cli_download_cmd(tmp_path, monkeypatch):
    from click.testing import CliRunner

    calls: list = []
    monkeypatch.setattr(urllib.request, "urlretrieve", _fake_urlretrieve(calls))

    target = tmp_path / "out.bgz"
    result = CliRunner().invoke(gnomad_cli.download_cmd, ["--output", str(target)])
    assert result.exit_code == 0, result.output
    assert target.exists()
    assert calls[0][0] == C.constraint_url("v2.1.1", "by_gene")


# --------------------------------------------------------------------------- #
# Real reachability (network-gated, deselected by default)
# --------------------------------------------------------------------------- #


@pytest.mark.network
@pytest.mark.parametrize("version", ["v2.1.1", "v4.0"])
def test_gnomad_constraint_urls_reachable(version):
    import urllib.request as _r

    url = C.constraint_url(version)
    req = _r.Request(url, method="HEAD")
    with _r.urlopen(req, timeout=30) as resp:
        assert resp.status == 200
        assert int(resp.headers.get("Content-Length", "0")) > 0
