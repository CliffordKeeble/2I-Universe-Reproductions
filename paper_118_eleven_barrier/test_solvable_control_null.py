"""
Acceptance tests for Paper 118 Part 2 (solvable-control null).

Guards the anchors, the pre-registered verdict (H_unique), the two-route agreement,
and the exploratory characterisation. Uses analyse()/build_population() directly so no
output files are written.

Run:  pytest -v test_solvable_control_null.py
"""

import pytest

import solvable_control_null as s


@pytest.fixture(scope="module")
def res():
    return [s.analyse(*g) for g in s.build_population()]


@pytest.fixture(scope="module")
def by_name(res):
    return {r["name"]: r for r in res}


def test_population_shape(res):
    assert len(res) == 91
    assert sum(1 for r in res if r["solvable"]) == 90
    assert sum(1 for r in res if not r["solvable"]) == 1


def test_machinery_anchors(by_name):
    # classical invariant degrees -- if these fail the machinery is wrong
    assert by_name["2T"]["l_gap"] == 6
    assert by_name["2O"]["l_gap"] == 8
    assert by_name["2I"]["l_gap"] == 12


def test_landed_on_H_unique(res):
    barriers = [r["name"] for r in res if r["barrier"]]
    assert barriers == ["2I"]
    assert not [r["name"] for r in res if r["barrier"] and r["solvable"]]


def test_known_barrier_values(by_name):
    assert by_name["2I"]["m_at_p1"] == 0        # the barrier
    assert by_name["2T"]["m_at_p1"] == 1        # return == gap
    assert by_name["2O"]["m_at_p1"] == 1        # return == gap
    assert by_name["2D5"]["m_at_p1"] == 1       # odd-n dihedral, nearest miss
    assert by_name["2I"]["lam1"] == 168


def test_nearest_solvable_approach_is_one(res):
    # no solvable group reaches m(p1) = 0; the closest is 1.
    assert min(r["m_at_p1"] for r in res if r["solvable"]) == 1


def test_exploratory_characterisation(res):
    # Barrier(Gamma)  <=>  p1 < l_gap, across the whole population.
    assert all(r["barrier"] == (r["p1"] < r["l_gap"]) for r in res)


@pytest.mark.parametrize("name", ["Z12", "Z25", "2D7", "2D30", "2T", "2O", "2I"])
def test_two_routes_agree(by_name, name):
    # analyse() already asserts A == B internally for every l; this re-checks a
    # sample independently against the exact route.
    r = by_name[name]
    if r["kind"] == "poly":
        exact = s.molien(*s.MOLIEN[name], r["L"])
        route_b = lambda l: exact[l]                       # noqa: E731
    elif r["kind"] == "cyclic":
        n = int(name[1:])
        route_b = lambda l: s.m_cyclic_exact(n, l)         # noqa: E731
    else:
        n = int(name[2:])
        route_b = lambda l: s.m_dihedral_exact(n, l)       # noqa: E731
    for l in range(r["L"] + 1):
        ma, resid = s.m_char(r_classes(name), r_order(name), l)
        assert ma == route_b(l) == r["m"][l]
        assert resid < 1e-30


def r_classes(name):
    for g in s.build_population():
        if g[0] == name:
            return g[3]
    raise KeyError(name)


def r_order(name):
    for g in s.build_population():
        if g[0] == name:
            return g[4]
    raise KeyError(name)


def test_m_zero_is_one_everywhere(res):
    assert all(r["m"][0] == 1 for r in res)
