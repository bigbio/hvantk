"""The builder's raw-input resolution -- the reprocess directory path.

``hvantk reprocess ensembl-gene:structure`` passes the builder the raw *directory* (no
parse stage is declared), while the plain-Python GTF parser needs the file. These tests
pin the resolver that bridges the two. Hail-free on purpose, so they run in the default
suite rather than only on a Hail node.
"""
from __future__ import annotations

from hvantk.resources.ensembl_release import ENSEMBL_GTF_FILENAME
from hvantk.skills.ensembl_gene.structure.builder import _resolve_gtf_path


def test_a_directory_resolves_to_the_committed_gtf_filename(tmp_path):
    (tmp_path / ENSEMBL_GTF_FILENAME).write_text("dummy")
    resolved = _resolve_gtf_path(tmp_path)
    assert resolved == str(tmp_path / ENSEMBL_GTF_FILENAME)


def test_a_file_path_is_returned_unchanged(tmp_path):
    gtf = tmp_path / "mini.gtf"
    gtf.write_text("dummy")
    assert _resolve_gtf_path(gtf) == str(gtf)
    assert _resolve_gtf_path(str(gtf)) == str(gtf)


def test_reprocess_directory_convention_reaches_a_real_file(tmp_path):
    """End-to-end path shape: reprocess hands over raw_dir; the resolver must land on the
    file the downloader wrote there (raw_dir/ENSEMBL_GTF_FILENAME)."""
    raw_dir = tmp_path / "raw"
    raw_dir.mkdir()
    downloaded = raw_dir / ENSEMBL_GTF_FILENAME
    downloaded.write_text("dummy")

    resolved = _resolve_gtf_path(raw_dir)  # what reprocess passes: the directory

    assert resolved == str(downloaded)
    import os

    assert os.path.isfile(resolved)
