"""INSIDER drift probe — documentation-only source (stub).

The Interactome Insider portal (http://interactomeinsider.yulab.org) hosts
the ``Whole_Human_Interactome_Interface_hg38.bed`` file behind a manual
download page with no machine-readable manifest or release feed; a real
drift probe would need to track the page contents or the file's
Last-Modified/ETag header. Writing that probe is out of scope for the
Phase 1 migration.

Until a real probe lands, ``fetch_fingerprint`` returns a structured stub
sentinel (via :func:`hvantk.core.plugin.api.stub_fingerprint`) so
``hvantk drift insider:variants`` reports a visible WARNING (status="stub")
instead of a silent false-green "clean". See issue #177.
"""

from __future__ import annotations

from hvantk.core.plugin.api import stub_fingerprint

NOT_IMPLEMENTED_REASON = (
    "insider drift probe not implemented; track upstream release manually"
)


def fetch_fingerprint() -> dict:
    return stub_fingerprint(NOT_IMPLEMENTED_REASON)
