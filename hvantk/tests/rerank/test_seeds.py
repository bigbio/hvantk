"""Pass C of #247: the seed is a parameter, and there is exactly one literal 42 left.

The CV partition is one of at least three variance components in the reported interval
(gene resampling, CV partition, best-of-N axis selection) and it was the only one no caller
could vary.
"""
from __future__ import annotations

import ast
from pathlib import Path

import numpy as np
import pytest

import hvantk.algorithms.rerank as rerank_pkg
from hvantk.algorithms.rerank.evaluator import _boot_ci, _raw_oof
from hvantk.algorithms.rerank.reranker import ReRanker, _gbm
from hvantk.algorithms.rerank.seeds import DEFAULT_SEED, rng_for
from hvantk.tests.rerank._synth import planted_signal

PACKAGE = Path(rerank_pkg.__file__).parent


def test_no_bare_42_outside_seeds_py():
    """The structural form of the plan's grep gate. Five literals existed before #247:
    reranker.py:10 and :32, evaluator.py:17 and :38, selection.py:177 -- five copies of one
    decision, which is five places to forget it."""
    offenders = {}
    for path in sorted(PACKAGE.rglob("*.py")):
        if path.name == "seeds.py":
            continue
        tree = ast.parse(path.read_text(), filename=str(path))
        hits = [
            node.lineno
            for node in ast.walk(tree)
            if isinstance(node, ast.Constant) and node.value == 42 and type(node.value) is int
        ]
        if hits:
            offenders[str(path.relative_to(PACKAGE))] = hits
    assert not offenders, offenders


def test_seeds_py_holds_exactly_one():
    tree = ast.parse((PACKAGE / "seeds.py").read_text())
    hits = [
        n for n in ast.walk(tree)
        if isinstance(n, ast.Constant) and n.value == 42 and type(n.value) is int
    ]
    assert len(hits) == 1
    assert DEFAULT_SEED == 42


def test_rng_for_is_offset_not_spawned():
    """A chunk must reproduce on its own: permutation i is seeded by seed + i however the
    work was divided."""
    assert rng_for(5, 3).random() == np.random.default_rng(8).random()


# --- the seed actually changes the partition ------------------------------------------------


def test_two_seeds_give_different_cv_partitions():
    matrix, y, baseline, _ = planted_signal(n=160, n_noise=0)
    a = _raw_oof(matrix, baseline, y, seed=DEFAULT_SEED)
    b = _raw_oof(matrix, baseline, y, seed=DEFAULT_SEED + 1)
    assert not np.allclose(a, b), "the seed did not reach the partition"


def test_the_default_seed_reproduces_the_historic_partition():
    matrix, y, baseline, _ = planted_signal(n=160, n_noise=0)
    assert np.allclose(_raw_oof(matrix, baseline, y), _raw_oof(matrix, baseline, y, seed=42))


def test_boot_ci_takes_a_seed():
    matrix, y, baseline, axes = planted_signal(n=160, n_noise=0)
    p0 = _raw_oof(matrix, baseline, y)
    p1 = _raw_oof(matrix, baseline + axes["axis0"], y)
    assert _boot_ci(y, p1, p0, n=200, seed=1) != _boot_ci(y, p1, p0, n=200, seed=2)
    assert _boot_ci(y, p1, p0, n=200, seed=1) == _boot_ci(y, p1, p0, n=200, seed=1)


def test_gbm_takes_a_seed():
    assert _gbm().random_state == DEFAULT_SEED
    assert _gbm(seed=7).random_state == 7


def test_reranker_takes_a_seed_and_uses_it():
    matrix, y, baseline, _ = planted_signal(n=160, n_noise=0)
    a = ReRanker(folds=5, seed=DEFAULT_SEED).score(matrix, baseline, y)
    b = ReRanker(folds=5, seed=DEFAULT_SEED + 1).score(matrix, baseline, y)
    assert not np.allclose(a, b)


def test_selection_policy_seed_defaults_to_the_shared_constant():
    from hvantk.algorithms.rerank.selection import SelectionPolicy

    assert SelectionPolicy().seed == DEFAULT_SEED


# --- config wiring ----------------------------------------------------------------------------


def test_config_seed_defaults_to_the_shared_constant():
    from hvantk.algorithms.rerank.config import Config

    assert Config.__dataclass_fields__["seed"].default == DEFAULT_SEED


@pytest.mark.parametrize("bad", ["42", 4.5, None, -1])
def test_config_rejects_a_non_integer_seed(bad):
    from hvantk.algorithms.rerank.config import Config

    with pytest.raises((TypeError, ValueError), match="seed"):
        Config(name="x", features=[], labels=None, seed=bad).__post_init__()
