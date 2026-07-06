"""
Deterministic tests for the Part 2 statistical machinery (permutation-z on
tie-corrected Spearman, and rank-partial correlation).

Run:  python -m pytest test_part2_stats.py -q
"""
import numpy as np

import part2_sweep as p2


def test_perm_z_perfect_negative_monotone():
    x = list(range(1, 12))
    y = list(range(11, 0, -1))            # perfect negative rank order
    rho, z, p = p2.perm_z(x, y, n_perm=2000, seed=1)
    assert rho < -0.999
    assert z < -3.0
    assert p < 0.01


def test_perm_z_constant_is_nan():
    x = [5] * 8
    y = list(range(8))
    rho, z, p = p2.perm_z(x, y)
    assert np.isnan(rho) and np.isnan(z)


def test_rank_partial_recovers_direct_when_control_independent():
    a = list(range(1, 11))
    b = list(range(1, 11))                # perfectly correlated
    control = [3, 9, 1, 7, 5, 10, 2, 8, 4, 6]  # an unrelated permutation
    assert p2.rank_partial(a, b, [control]) > 0.8


def test_rank_partial_removes_shared_confound():
    # a and b share a common control but have independent noise -> partial ~ 0,
    # even though their raw correlation is strongly positive.
    rng = np.random.default_rng(0)
    control = rng.normal(size=300)
    a = control + rng.normal(size=300)
    b = control + rng.normal(size=300)
    raw = np.corrcoef(a, b)[0, 1]
    partial = p2.rank_partial(list(a), list(b), [list(control)])
    assert raw > 0.3                     # confounded raw correlation
    assert abs(partial) < 0.15           # removed once control is partialled out


def test_rank_partial_degenerate_is_nan():
    a = list(range(1, 11))
    assert np.isnan(p2.rank_partial(a, a, [a]))
