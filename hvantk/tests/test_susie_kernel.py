"""Hermetic regression for the pure-Python SuSiE-RSS + coloc.susie kernel (issue #193).

Network-free and R-free: validates the CONFIRMED-vs-REFUTED discrimination — the
whole point of the fine-mapping layer — on synthetic z-scores over a simulated
two-block LD matrix. Runs in the default ``pytest -q``.

- Shared causal variant    -> a credible set in each trait, high PP4 (CONFIRM).
- Distinct causal variants  -> a credible set in each trait, PP4 ~ 0 (REFUTE: distinct).
- One trait with no signal  -> no credible set, PP4 = 0 (REFUTE: artifact).

The end-to-end controls on real data (AF->MYOZ1 CONFIRMED, CHD 17q21/NSF REFUTED)
live in the network-marked ``test_gwas_coloc_*_network.py``.
"""

import numpy as np

from hvantk.algorithms.qtlcascade.susie import susie_rss, coloc_susie, SusieFit


HALF = 20          # variants per LD block
N_REF = 600        # reference panel size
N_STUDY = 50_000   # study sample size for the z-score regularisation
C_A, C_B = 5, 25   # causal indices in block A and block B


def _sim_block(rng, n, k, rho, af=0.3):
    """k variants sharing a latent dosage with prob ``rho`` -> within-block LD."""
    base = rng.binomial(2, af, size=n).astype(float)
    cols = []
    for _ in range(k):
        ind = rng.binomial(2, af, size=n).astype(float)
        cols.append(np.where(rng.rand(n) < rho, base, ind))
    return np.array(cols)


def _ld_matrix():
    """Two independent LD blocks (within-block |r| high, cross-block ~0)."""
    rng = np.random.RandomState(0)
    A = _sim_block(rng, N_REF, HALF, 0.9)
    B = _sim_block(rng, N_REF, HALF, 0.9)
    return np.corrcoef(np.vstack([A, B]))


def _z_from(R, causal, z0):
    """Marginal z-scores for a single causal variant propagated through LD."""
    return R[:, causal] * z0


def test_ld_blocks_are_separable():
    R = _ld_matrix()
    within = np.abs(R[:HALF, :HALF][~np.eye(HALF, dtype=bool)]).min()
    cross = np.abs(R[:HALF, HALF:]).max()
    assert within >= 0.5      # purity filter (min_abs_corr) will pass
    assert cross < 0.3        # blocks are effectively independent


def test_shared_causal_confirms():
    R = _ld_matrix()
    sg = susie_rss(_z_from(R, C_A, 9.0), R, N_STUDY)
    se = susie_rss(_z_from(R, C_A, 8.0), R, N_STUDY)
    assert len(sg.cs) == 1
    assert len(se.cs) == 1
    # Each credible set covers the shared causal variant.
    assert C_A in sg.cs[0] and C_A in se.cs[0]
    assert coloc_susie(sg, se) > 0.9      # CONFIRMED analog


def test_distinct_causals_refute():
    R = _ld_matrix()
    sg = susie_rss(_z_from(R, C_A, 9.0), R, N_STUDY)   # causal in block A
    se = susie_rss(_z_from(R, C_B, 8.0), R, N_STUDY)   # causal in block B
    assert len(sg.cs) == 1 and len(se.cs) == 1
    assert coloc_susie(sg, se) < 0.1      # REFUTED (distinct causal variants)


def test_no_signal_yields_no_credible_set():
    R = _ld_matrix()
    sg = susie_rss(_z_from(R, C_A, 9.0), R, N_STUDY)
    se = susie_rss(np.full(R.shape[0], 0.1), R, N_STUDY)  # flat: no signal
    assert len(se.cs) == 0
    assert coloc_susie(sg, se) == 0.0     # REFUTED (single-variant artifact)


def test_coloc_susie_zero_when_no_credible_set():
    empty = SusieFit(alpha=np.zeros((0, 0)), lbf_variable=np.zeros((0, 0)),
                     V=np.zeros(0))
    other = SusieFit(alpha=np.zeros((1, 3)), lbf_variable=np.zeros((1, 3)),
                     V=np.array([0.1]), cs=[np.array([0])], cs_effect=[0])
    assert coloc_susie(empty, other) == 0.0
    assert coloc_susie(other, empty) == 0.0


def test_susie_rss_handles_tiny_input():
    R = np.array([[1.0, 0.2], [0.2, 1.0]])
    fit = susie_rss(np.array([0.3, -0.1]), R, 1000)
    assert fit.alpha.shape == (2, 2)      # L capped at p
    assert len(fit.cs) == 0               # no real signal -> no credible set
