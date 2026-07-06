"""
pytest suite for the Paper 116 pipeline.

Covers AT1 (multiplicity identity, two independent routes and two independent
modules), the A1/F1 reference moments, the rotation-angle vs half-angle
reconciliation (the assertion the brief asks for, not a comment), the GUE/GOE
switch direction, discriminant self-consistency, and pipeline determinism.

Run:  python -m pytest -q     (from this folder)
"""
import numpy as np
import pytest
from mpmath import mpf, pi, sin

import mult
import spectral_stats as ss
from mult import eb

SEED = 20260706

# groups spanning all three families for the route cross-check
AT1_GROUPS = [("2I", "poly", None), ("2O", "poly", None), ("2T", "poly", None),
              ("Z12", "cyclic", 12), ("Z120", "cyclic", 120),
              ("2D5", "dihedral", 5), ("2D15", "dihedral", 15)]


# ----------------------------------------------------------------------
# AT1 -- the multiplicities feeding Z_Gamma equal the repo's, exactly.
# ----------------------------------------------------------------------
@pytest.mark.parametrize("name,kind,param", AT1_GROUPS)
def test_at1_route_A_equals_route_B(name, kind, param):
    L = 60
    route_a = mult.m_char_series(name, kind, param, L)          # mpmath char sum
    route_b = mult.m_exact_series(name, kind, param, None, None, L)  # exact int
    assert route_a == route_b


def test_at1_2I_two_modules_agree():
    L = 120
    a = mult.m_2I_eleven_barrier(L)                       # eleven_barrier Molien
    b = mult.m_exact_series("2I", "poly", None, None, None, L)  # scn Molien
    assert a == b
    assert b[0] == 1
    assert b[10] == 0 and b[12] == 1     # the 11 barrier: gap at l=12, not l=10


# ----------------------------------------------------------------------
# A1/F1 -- reference surmises have unit mean and the corrected variances.
# ----------------------------------------------------------------------
@pytest.mark.parametrize("name,var_ref", [("Poisson", 1.0),
                                          ("GOE", 4 / np.pi - 1),
                                          ("GUE", 3 * np.pi / 8 - 1)])
def test_reference_moments(name, var_ref):
    surv = 1.0 - ss._CDF[name]
    grid = ss._SGRID
    mean = np.sum(0.5 * (surv[1:] + surv[:-1]) * np.diff(grid))      # int(1-F)ds
    m2 = np.sum(0.5 * ((2 * grid * surv)[1:] + (2 * grid * surv)[:-1])
                * np.diff(grid))                                     # int 2s(1-F)ds
    assert abs(mean - 1.0) < 2e-3
    assert abs((m2 - mean**2) - var_ref) < 3e-3


# ----------------------------------------------------------------------
# Reconciliation -- repo half-angle theta = (rotation angle)/2 (assertion).
# 116 s2.1 writes chi_l at the rotation angle alpha: sin((l+1)alpha/2)/sin(alpha/2).
# The repo stores theta = p*pi/q with alpha = 2*p*pi/q, so theta = alpha/2.
# ----------------------------------------------------------------------
def test_rotation_angle_half_angle_reconciliation():
    # non-endpoint 2I classes (skip Gen: theta=0, pi)
    classes = [(p, q) for _lab, p, q, _sz, fib in eb.CLASSES_2I if fib != "Gen"]
    for (p, q) in classes:
        alpha = 2 * mpf(p) / q * pi          # rotation angle
        for l in range(0, 25):
            paper = sin((l + 1) * alpha / 2) / sin(alpha / 2)   # 116 s2.1 form
            repo = eb.chi_mp(l, p, q)                            # repo half-angle
            assert abs(paper - repo) < mpf("1e-40")


# ----------------------------------------------------------------------
# The GUE/GOE switch direction must reproduce (the core claim).
# ----------------------------------------------------------------------
def test_switch_direction():
    m2I = mult.m_exact_series("2I", "poly", None, None, None, 200)
    res_2I, _z, _s = ss.run_series(m2I, target_N=250, l_max=150)
    assert res_2I["class_cvm"] == "GUE"
    assert res_2I["ratio_cvm"] < 1.0

    mZ = mult.m_exact_series("Z120", "cyclic", 120, None, None, 600)
    res_Z, _z2, _s2 = ss.run_series(mZ, target_N=250, l_max=400)
    assert res_Z["class_cvm"] == "GOE"
    assert res_Z["ratio_cvm"] > 1.0


def test_2I_sub_GUE_variance():
    """2I spacing variance sits below the GUE surmise (0.178) -- the sub-GUE regime."""
    m2I = mult.m_exact_series("2I", "poly", None, None, None, 200)
    res, _z, _s = ss.run_series(m2I, target_N=400, l_max=200)
    assert res["var"] < ss.REF_VAR["GUE"]


# ----------------------------------------------------------------------
# Pipeline determinism -- no hidden RNG; identical inputs, identical spacings.
# ----------------------------------------------------------------------
def test_determinism():
    m2I = mult.m_exact_series("2I", "poly", None, None, None, 150)
    _r1, _z1, s1 = ss.run_series(m2I, l_max=150, t_max=120.0)
    _r2, _z2, s2 = ss.run_series(m2I, l_max=150, t_max=120.0)
    assert np.array_equal(s1, s2)


# ----------------------------------------------------------------------
# Discriminant self-consistency -- surmise-sampled spacings classify to self.
# ----------------------------------------------------------------------
@pytest.mark.parametrize("name", ["Poisson", "GOE", "GUE"])
def test_classifier_recovers_reference(name):
    rng = np.random.default_rng(SEED)
    u = rng.random(4000)
    s = np.interp(u, ss._CDF[name], ss._SGRID)      # inverse-CDF sampling
    s = s / np.mean(s)
    cls, _ratio, _msd = ss.classify_cvm(s)
    assert cls == name


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
