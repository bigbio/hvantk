"""The one place the default random seed is written down.

Every stochastic step in this package carried its own literal ``42``: the CV partition in
``evaluator._raw_oof`` and again in ``ReRanker.score``, the estimator's own
``random_state`` in ``_gbm``, the bootstrap generator in ``_boot_ci``, and
``SelectionPolicy.seed``. Five copies of one decision is five places to forget it, and the
consequence is not cosmetic: the CV partition is one of at least three variance components
in the reported interval (gene resampling, CV partition, best-of-N axis selection) and it
was the only one no caller could vary.

``hvantk/tests/rerank/test_seeds.py`` greps the package and fails if a bare ``42``
reappears outside this file.
"""
from __future__ import annotations

import numpy as np

DEFAULT_SEED = 42
"""The historic value, kept so every result produced before the seed became a parameter
reproduces exactly. It carries no other meaning; callers should feel free to change it."""


def rng_for(seed: int, offset: int = 0) -> np.random.Generator:
    """``default_rng(seed + offset)``.

    ``seed + i`` is the historic scheme for seeding permutation ``i``, kept so nulls
    computed before this module existed reproduce exactly. It is chunk-independent --
    rerunning one chunk reproduces exactly the permutations the whole run would have drawn
    for it -- and its one cost is the seed-proximity caveat below.

    Caveat: nulls from seeds closer together than ``n_perm`` share permutations (``seed=42``
    and ``seed=43`` share ``n_perm - 1`` of their draws), so an independent replicate null
    needs its seed moved by at least ``n_perm``.
    """
    return np.random.default_rng(int(seed) + int(offset))
