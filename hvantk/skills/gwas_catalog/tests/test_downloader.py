"""
Unit test for the GWAS Catalog downloader.

Mocks the network call to avoid an actual download, standing a tiny
in-memory zip (built with ``zipfile`` + ``io.BytesIO``) in for the real
~70 MB EBI release archive. One happy-path test only, per review: the brief
called for a single mocked downloader test, mirroring
``hvantk/skills/hgnc/tests/test_downloader.py``.
"""

from __future__ import annotations

import io
import os
import stat
import zipfile
from unittest.mock import MagicMock, patch

from hvantk.core.plugin.drift_runner import _current_umask
from hvantk.skills.gwas_catalog.cli import (
    GWAS_CATALOG_ASSOCIATIONS_FILENAME,
    download_dataset,
)

_TSV_CONTENT = b"DATE ADDED TO CATALOG\tPUBMEDID\n2026-01-01\t12345678\n"


def test_download_dataset_extracts_the_tsv_under_the_canonical_name(tmp_path):
    """Happy path: a mocked zip response is streamed, the single member it
    contains is extracted to ``<raw_dir>/<canonical filename>``, and no temp
    files are left behind in raw_dir once the call returns."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("gwas-catalog-download-associations-v1.0-full.tsv", _TSV_CONTENT)

    response = MagicMock()
    response.raise_for_status.return_value = None
    response.iter_content.return_value = [buf.getvalue()]
    response.__enter__.return_value = response
    response.__exit__.return_value = False

    with patch("requests.get", return_value=response):
        result = download_dataset(str(tmp_path))

    expected = tmp_path / GWAS_CATALOG_ASSOCIATIONS_FILENAME
    assert result == str(expected)
    assert expected.read_bytes() == _TSV_CONTENT
    assert list(tmp_path.iterdir()) == [expected]  # no leftover temp files

    # mkstemp always creates mode 0600; the published file must not keep that --
    # it should land at whatever a normal file create gets (0666 masked by the
    # process umask), so a shared raw dir stays group-readable.
    assert stat.S_IMODE(os.stat(expected).st_mode) == (0o666 & ~_current_umask())
