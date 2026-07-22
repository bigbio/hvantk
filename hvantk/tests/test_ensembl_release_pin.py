"""The PTM algorithm must reuse the single Ensembl release pin, not define its own.

This assertion lives at the top-level test layer rather than under
``hvantk/skills/ensembl_gene/tests/`` because it imports ``hvantk.algorithms.ptm``, and a
file under ``skills/`` may not import from ``algorithms/`` (test_dependency_directions).
The pin itself lives in ``hvantk/resources/ensembl_release.py`` -- substrate that both the
skill and the algorithm may depend on.
"""
from __future__ import annotations

from hvantk.resources.ensembl_release import (
    ENSEMBL_GTF_FILENAME,
    ENSEMBL_GTF_URL,
    ENSEMBL_RELEASE,
)


def test_ptm_constants_reuse_the_same_pin():
    from hvantk.algorithms.ptm import constants

    assert constants.ENSEMBL_RELEASE == ENSEMBL_RELEASE
    assert constants.ENSEMBL_GTF_URL == ENSEMBL_GTF_URL
    assert constants.ENSEMBL_GTF_FILENAME == ENSEMBL_GTF_FILENAME
