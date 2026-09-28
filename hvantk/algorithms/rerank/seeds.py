"""The one place the default random seed is written down.

Every stochastic step in this package carried its own literal ``42``: the CV partition in
``evaluator._raw_oof`` and again in ``ReRanker.score``, the estimator's own
``random_state`` in ``_gbm``, the bootstrap generator in ``_boot_ci``, and
``SelectionPolicy.seed``. Five copies of one decision is five places to forget it, and the
consequence is not cosmetic: the CV partition is one of at least three variance components
in the reported interval (gene resampling, CV partition, best-of-N axis selection) and it
was the only one no caller could vary. Measured on the CHD lattice, a single-seed run
scored 0.765 against a 25-seed mean of 0.784
(``analysis/baseline-lattice/nested_check.py:161-175``) -- the published figure was the
unluckiest of 25 draws, and nothing in the API let a user notice.

``hvantk/tests/rerank/test_seeds.py`` greps the package and fails if a bare ``42``
reappears outside this file.
"""
from __future__ import annotations

import numpy as np

DEFAULT_SEED = 42
"""The historic value, kept so every result produced before the seed became a parameter
reproduces exactly. It carries no other meaning; callers should feel free to change it."""


def rng_for(seed: int, offset: int = 0) -> np.random.Generator:
    """``default_rng(seed + offset)`` -- deliberately NOT ``SeedSequence.spawn``.

    The permutation null is chunked across array tasks, and a chunk has to be reproducible
    on its own: permutation ``i`` is seeded by ``seed + i`` whichever chunk computes it, so
    rerunning chunk 3 reproduces exactly permutations ``[lo, hi)``
    (``analysis/rerank-homogenised/perm_null_arm.py:77-86``). Spawning from a parent
    sequence makes draw ``i`` depend on how the work was divided, which is the one property
    a chunkable null cannot have.
    """
    return np.random.default_rng(int(seed) + int(offset))
