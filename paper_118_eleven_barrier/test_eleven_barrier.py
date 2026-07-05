"""
Acceptance tests for Paper 118 Part 1 (Brief 118-CODE-1).

AT1-AT9 reproduce the DERIVED tier of "The 11 Barrier"; each load-bearing number
is checked across all three routes (character sum / Molien / exact Q(sqrt5)).

Run:  pytest -v test_eleven_barrier.py
"""

from fractions import Fraction as F

from mpmath import mpf

import eleven_barrier as eb
from eleven_barrier import Q5, PHI, PSI, F as _F  # noqa: F401


# --------------------------------------------------------------------------
# Substrate: class sizes and the three routes agree everywhere
# --------------------------------------------------------------------------
def test_class_sizes_sum_to_order():
    assert eb.assert_class_sizes() == 120


def test_three_routes_agree_and_integral():
    """(a) mpmath, (b) Molien, (c) exact Q(sqrt5) agree for 0 <= l <= 120,
    and route (a)'s integrality residual is below 1e-30."""
    mol = eb.molien_coeffs(120)
    for l in range(121):
        mc, residual = eb.m_char(l)
        assert residual < mpf("1e-30"), f"l={l}: integrality residual {residual}"
        assert mc == mol[l] == eb.m_exact(l), f"route disagreement at l={l}"


# --------------------------------------------------------------------------
# AT1 -- the barrier and the gap
# --------------------------------------------------------------------------
def test_AT1_barrier_and_gap():
    assert eb.m_exact(0) == 1
    for l in range(1, 12):
        assert eb.m_exact(l) == 0, f"expected m({l}) = 0"
    assert eb.m_exact(12) == 1
    gap_l, gap_lambda = eb.spectral_gap()
    assert gap_l == 12
    assert gap_lambda == 12 * 14 == 168


# --------------------------------------------------------------------------
# AT2 -- the exact cancellation at l = 10
# --------------------------------------------------------------------------
def test_AT2_cancellation_at_10():
    fib = eb.fibre_decomposition(10)
    assert fib["Gen"] == F(11, 60)
    assert fib["Vtx"] == F(2, 5)
    assert fib["Face"] == F(-1, 3)
    assert fib["Edge"] == F(-1, 4)
    # Frobenius sum in 120ths: 22 + 48 - 40 - 30 = 0.
    assert (F(22, 120) + F(48, 120) - F(40, 120) - F(30, 120)) == 0
    assert sum(fib.values()) == 0
    assert eb.m_exact(10) == 0


# --------------------------------------------------------------------------
# AT3 -- every vertex character equals 1 at l = 10
# --------------------------------------------------------------------------
def test_AT3_vertex_characters_at_10():
    chars = eb.vertex_characters(10)
    assert len(chars) == 4
    assert all(c == Q5(1, 0) for c in chars)


# --------------------------------------------------------------------------
# AT4 -- the golden consensus at l = 12
# --------------------------------------------------------------------------
def test_AT4_decomposition_at_12():
    fib = eb.fibre_decomposition(12)
    assert fib["Gen"] == F(13, 60)
    assert fib["Vtx"] == F(1, 5)
    assert fib["Face"] == F(1, 3)
    assert fib["Edge"] == F(1, 4)
    assert sum(fib.values()) == 1


def test_AT4_golden_vertex_values_at_12():
    chars = eb.vertex_characters(12)
    # {psi, phi, phi, psi} in some class order.
    assert sorted((c.a, c.b) for c in chars) == \
        sorted((c.a, c.b) for c in [PHI, PHI, PSI, PSI])
    # sum = 2(phi + psi) = 2, exact in Q(sqrt5).
    total = Q5(0, 0)
    for c in chars:
        total = total + c
    assert total == Q5(2, 0)


# --------------------------------------------------------------------------
# AT5 -- Table 5: vertex returns within one full period (60)
# --------------------------------------------------------------------------
def test_AT5_table5():
    # l -> (m, face_slot, face_value, edge_slot, edge_value)
    expected = {
        0:  (1, 0, F(1, 3),  0, F(1, 4)),
        10: (0, 4, F(-1, 3), 2, F(-1, 4)),
        20: (1, 2, F(0),     0, F(1, 4)),
        30: (1, 0, F(1, 3),  2, F(-1, 4)),
        40: (1, 4, F(-1, 3), 0, F(1, 4)),
        50: (1, 2, F(0),     2, F(-1, 4)),
    }
    for l, (m, fslot, fval, eslot, eval_) in expected.items():
        fib = eb.fibre_decomposition(l)
        assert eb.m_exact(l) == m, f"m({l})"
        assert l % 6 == fslot, f"face slot at l={l}"
        assert l % 4 == eslot, f"edge slot at l={l}"
        assert fib["Face"] == fval, f"Face({l})"
        assert fib["Edge"] == eval_, f"Edge({l})"


# --------------------------------------------------------------------------
# AT6 -- Molien identity, coefficient by coefficient to l = 120
# --------------------------------------------------------------------------
def test_AT6_molien_identity():
    mol = eb.molien_coeffs(120)
    for l in range(121):
        assert eb.m_exact(l) == mol[l], f"Molien mismatch at l={l}"
    # spot-check the closed form's structure: only 12,20,30,32,... carry weight
    assert mol[0] == 1 and mol[12] == 1 and mol[20] == 1 and mol[30] == 1
    assert mol[10] == 0 and mol[11] == 0


# --------------------------------------------------------------------------
# AT7 -- fibre cycles match Table 1 exactly for l = 0..59
# --------------------------------------------------------------------------
VTX_CYCLE = [F(2, 5), F(0), F(1, 5), F(0), F(0),
             F(0), F(-1, 5), F(0), F(-2, 5), F(0)]          # period 10
FACE_CYCLE = [F(1, 3), F(0), F(0), F(0), F(-1, 3), F(0)]    # period 6
EDGE_CYCLE = [F(1, 4), F(0), F(-1, 4), F(0)]                # period 4


def test_AT7_fibre_cycles():
    for l in range(60):
        fib = eb.fibre_decomposition(l)
        assert fib["Vtx"] == VTX_CYCLE[l % 10], f"Vtx({l})"
        assert fib["Face"] == FACE_CYCLE[l % 6], f"Face({l})"
        assert fib["Edge"] == EDGE_CYCLE[l % 4], f"Edge({l})"
        expected_gen = F(l + 1, 60) if l % 2 == 0 else F(0)
        assert fib["Gen"] == expected_gen, f"Gen({l})"


# --------------------------------------------------------------------------
# AT8 -- global uniqueness of the barrier among vertex returns to l = 120
# --------------------------------------------------------------------------
def test_AT8_unique_barrier():
    zeros = [l for l in range(0, 121, 10) if eb.m_exact(l) == 0]
    assert zeros == [10], f"expected only l=10 to vanish, got {zeros}"
    # l = 70: same destructive face/edge, but the generic baseline compensates.
    fib = eb.fibre_decomposition(70)
    assert fib["Gen"] == F(71, 60)
    assert fib["Vtx"] == F(2, 5)
    assert fib["Face"] == F(-1, 3)
    assert fib["Edge"] == F(-1, 4)
    assert eb.m_exact(70) == 1
    assert (F(71, 60) + F(24, 60) - F(20, 60) - F(15, 60)) == 1


# --------------------------------------------------------------------------
# AT9 -- characterisation of l = 12 (Paper 118 section 11.2)
# --------------------------------------------------------------------------
def test_AT9_twelve_characterisation():
    def qualifies(l):
        return l % 4 == 0 and l % 6 == 0 and eb.fibre_decomposition(l)["Vtx"] > 0

    smallest = next(l for l in range(1, 200) if qualifies(l))
    assert smallest == 12
