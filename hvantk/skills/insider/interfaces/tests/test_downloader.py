"""Unit test for the INSIDER interfaces downloader.

Mirrors ``hvantk/skills/hgnc/tests/test_downloader.py``'s mocked-network approach
(no real HTTP calls), kept to the one scenario specific to this downloader: a
streamed, atomically-written download that still honours the same
skip-if-present/overwrite semantics hgnc's downloader uses.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from hvantk.skills.insider.interfaces.builder import DEFAULT_FILENAME
from hvantk.skills.insider.interfaces.cli import download_dataset
from hvantk.skills.insider.interfaces.drift_probe import INSIDER_INTERFACES_URL


def test_download_dataset_streams_atomically_and_respects_overwrite(tmp_path):
    """Writes via temp-then-rename to the filename the builder expects, and
    refuses to clobber an existing file unless overwrite=True."""
    raw_dir = tmp_path / "raw"
    content = b"P1\tP2\tSource\tP1_IRES\tP2_IRES\nQ1\tQ2\tECLAIR\t[1]\t[]\n"

    def _fake_get(url, stream=True, timeout=None, headers=None):
        assert url == INSIDER_INTERFACES_URL
        response = MagicMock()
        response.raise_for_status.return_value = None
        # Split across chunks to exercise the streaming write loop.
        response.iter_content.return_value = [content[:10], content[10:]]
        response.headers = {}
        # download_dataset delegates to hvantk.core.utils.http.stream_to_file, which
        # uses the response as a context manager (`with requests.get(...) as
        # response:`); MagicMock's auto-created __enter__ returns a DIFFERENT,
        # unconfigured mock unless wired to return self.
        response.__enter__.return_value = response
        response.__exit__.return_value = False
        return response

    with patch("requests.get", side_effect=_fake_get) as mock_get:
        target = download_dataset(raw_dir=str(raw_dir))
        mock_get.assert_called_once()

    target_path = Path(target)
    assert target_path == raw_dir / DEFAULT_FILENAME
    assert target_path.read_bytes() == content
    # Atomic write: no leftover temp file of any name once the download completes.
    assert list(raw_dir.iterdir()) == [target_path]

    # Re-downloading without overwrite must not clobber the existing file (and
    # must not even attempt a network call).
    with pytest.raises(FileExistsError):
        download_dataset(raw_dir=str(raw_dir))
