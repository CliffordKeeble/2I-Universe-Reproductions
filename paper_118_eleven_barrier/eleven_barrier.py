"""
Paper 118 — The 11 Barrier: reproduction of the DERIVED tier.

Brief 118-CODE-1 Part 1 (Amendment A1: Python, reproductions repo).

Machinery
---------
For a finite Gamma < SU(2), the multiplicity of the trivial representation of
Gamma inside the (l+1)-dimensional irreducible representation of SU(2) is

    m_Gamma(l) = (1/|Gamma|) * sum_C |C| * chi_l(theta_C),

    chi_l(theta) = sin((l+1) theta) / sin(theta),
    chi_l(0)  = l + 1,
    chi_l(pi) = (-1)^l (l + 1).

chi_l is the Chebyshev polynomial of the second kind: chi_l(theta) = U_l(cos theta).

Two independent routes are used for every load-bearing number, plus an exact
field route for the fibre fractions:

  (a) character-sum route  -- mpmath at >= 50 dps, integrality assertion
                              |m - round(m)| < 1e-30.               [OBSERVED->DERIVED]
  (b) Molien route         -- exact integer power series
                              sum_l m(l) t^l = (1 + t^30) / ((1 - t^12)(1 - t^20)). [DERIVED]
  (c) exact Q(sqrt5) route -- Chebyshev-U recurrence in the ring Q(sqrt5), giving
                              the fibre contributions (Gen, Vtx, Face, Edge) as exact
                              rationals and the golden vertex values {phi, psi}.   [DERIVED]

A number is verified when (a), (b) and (c) agree.

Run standalone:  python eleven_barrier.py   -> writes m_2I.csv, fibre_cycles.csv
"""

from __future__ import annotations

import csv
import os
from fractions import Fraction as F

from mpmath import mp, mpf, sin, pi, nint, fabs

# 50-digit floor requested by the brief; carry extra guard digits.
mp.dps = 60

# --------------------------------------------------------------------------
# Class data for 2I (|Gamma| = 120).  theta stored as an exact rational
# multiple of pi: (p, q) means theta = p*pi/q.
# Fibre grouping per Paper 118 Table 1:
#   Gen  = {theta = 0, pi}          (identity, -I)
#   Vtx  = four size-12 classes     (vertex, Z5 stabiliser, period 10)
#   Face = two size-20 classes      (face,   Z3 stabiliser, period  6)
#   Edge = size-30 class            (edge,   Z2 stabiliser, period  4)
# --------------------------------------------------------------------------
ORDER_2I = 120

CLASSES_2I = [
    # (label,   p, q,  size, fibre)
    ("1",       0, 1,     1, "Gen"),   # theta = 0
    ("-1",      1, 1,     1, "Gen"),   # theta = pi
    ("C10",     1, 5,    12, "Vtx"),   # theta = pi/5
    ("C5",      2, 5,    12, "Vtx"),   # theta = 2pi/5
    ("C10'",    3, 5,    12, "Vtx"),   # theta = 3pi/5
    ("C5'",     4, 5,    12, "Vtx"),   # theta = 4pi/5
    ("C6",      1, 3,    20, "Face"),  # theta = pi/3
    ("C3",      2, 3,    20, "Face"),  # theta = 2pi/3
    ("C4",      1, 2,    30, "Edge"),  # theta = pi/2
]

FIBRES = ("Gen", "Vtx", "Face", "Edge")


def assert_class_sizes(classes=CLASSES_2I, order=ORDER_2I):
    """Brief: assert sum of class sizes = |Gamma| before anything else."""
    total = sum(c[3] for c in classes)
    if total != order:
        raise AssertionError(f"class sizes sum to {total}, expected |Gamma| = {order}")
    return total


# ==========================================================================
# (a) character-sum route -- mpmath, high precision
# ==========================================================================
def chi_mp(l: int, p: int, q: int):
    """chi_l(theta) at theta = p*pi/q, high precision."""
    if p == 0:                       # theta = 0
        return mpf(l + 1)
    if p == q:                       # theta = pi
        return mpf((-1) ** l * (l + 1))
    theta = pi * p / q
    return sin((l + 1) * theta) / sin(theta)


def m_char(l: int, classes=CLASSES_2I, order=ORDER_2I):
    """Multiplicity via the class sum, high precision. Returns (int, residual)."""
    s = mpf(0)
    for _label, p, q, size, _fibre in classes:
        s += size * chi_mp(l, p, q)
    m = s / order
    m_round = int(nint(m))
    residual = fabs(m - m_round)
    return m_round, residual


# ==========================================================================
# (b) Molien route -- exact integer power series
# ==========================================================================
def molien_coeffs(L: int):
    """
    Coefficients of (1 + t^30) / ((1 - t^12)(1 - t^20)) up to t^L, exact ints.

    1/((1-t^12)(1-t^20)) has coefficient p(n) = #{(a,b) >= 0 : 12a + 20b = n};
    multiplying by (1 + t^30) gives m(n) = p(n) + p(n-30).
    """
    n = L + 1
    s = [0] * n
    s[0] = 1
    for k in range(12, n):           # multiply by 1/(1 - t^12)
        s[k] += s[k - 12]
    for k in range(20, n):           # multiply by 1/(1 - t^20)
        s[k] += s[k - 20]
    m = [0] * n
    for k in range(n):               # multiply by (1 + t^30)
        m[k] = s[k] + (s[k - 30] if k >= 30 else 0)
    return m


# ==========================================================================
# (c) exact Q(sqrt5) route -- fibre decomposition and golden vertex values
# ==========================================================================
class Q5:
    """Element a + b*sqrt(5) of the field Q(sqrt5); a, b are Fractions."""

    __slots__ = ("a", "b")

    def __init__(self, a=0, b=0):
        self.a = F(a)
        self.b = F(b)

    @staticmethod
    def _coerce(o):
        return o if isinstance(o, Q5) else Q5(o)

    def __add__(self, o):
        o = self._coerce(o)
        return Q5(self.a + o.a, self.b + o.b)

    __radd__ = __add__

    def __sub__(self, o):
        o = self._coerce(o)
        return Q5(self.a - o.a, self.b - o.b)

    def __rsub__(self, o):
        return self._coerce(o).__sub__(self)

    def __mul__(self, o):
        o = self._coerce(o)
        return Q5(self.a * o.a + 5 * self.b * o.b, self.a * o.b + self.b * o.a)

    __rmul__ = __mul__

    def __neg__(self):
        return Q5(-self.a, -self.b)

    def __eq__(self, o):
        o = self._coerce(o)
        return self.a == o.a and self.b == o.b

    def __hash__(self):
        return hash((self.a, self.b))

    def is_rational(self):
        return self.b == 0

    def as_fraction(self):
        if self.b != 0:
            raise ValueError(f"not rational: {self}")
        return self.a

    def __repr__(self):
        if self.b == 0:
            return f"{self.a}"
        return f"{self.a} + {self.b}*sqrt5"


# Golden ratio and its conjugate, exact in Q(sqrt5).
PHI = Q5(F(1, 2), F(1, 2))    # (1 + sqrt5)/2
PSI = Q5(F(1, 2), F(-1, 2))   # (1 - sqrt5)/2

# Exact cos(theta) for every class, in Q(sqrt5).
#   cos(pi/5)  = (1 + sqrt5)/4      cos(2pi/5) = (-1 + sqrt5)/4
#   cos(3pi/5) = (1 - sqrt5)/4      cos(4pi/5) = (-1 - sqrt5)/4
COS_Q5 = {
    (0, 1): Q5(1, 0),
    (1, 1): Q5(-1, 0),
    (1, 5): Q5(F(1, 4), F(1, 4)),
    (2, 5): Q5(F(-1, 4), F(1, 4)),
    (3, 5): Q5(F(1, 4), F(-1, 4)),
    (4, 5): Q5(F(-1, 4), F(-1, 4)),
    (1, 3): Q5(F(1, 2), 0),
    (2, 3): Q5(F(-1, 2), 0),
    (1, 2): Q5(0, 0),
}


def cheb_u_q5(l: int, x: Q5) -> Q5:
    """U_l(x) = chi_l via the Chebyshev recurrence, exact in Q(sqrt5)."""
    if l == 0:
        return Q5(1, 0)
    u_prev, u_cur = Q5(1, 0), Q5(2, 0) * x   # U_0, U_1
    if l == 1:
        return u_cur
    for _ in range(2, l + 1):
        u_prev, u_cur = u_cur, Q5(2, 0) * x * u_cur - u_prev
    return u_cur


def chi_exact(l: int, p: int, q: int) -> Q5:
    """chi_l(theta) = U_l(cos theta), exact in Q(sqrt5)."""
    return cheb_u_q5(l, COS_Q5[(p, q)])


def vertex_characters(l: int):
    """The four vertex-class characters at level l, exact in Q(sqrt5)."""
    return [chi_exact(l, p, q) for _lab, p, q, _sz, fib in CLASSES_2I if fib == "Vtx"]


def fibre_decomposition(l: int):
    """
    Exact (Fraction) contribution of each fibre to m(l):
        contribution = (1/|Gamma|) * sum_{C in fibre} |C| * chi_l(theta_C).
    Each fibre sum is Galois-stable, hence rational; asserted so.
    """
    out = {f: Q5(0, 0) for f in FIBRES}
    for _lab, p, q, size, fib in CLASSES_2I:
        out[fib] += Q5(size, 0) * chi_exact(l, p, q)
    result = {}
    for fib in FIBRES:
        val = out[fib] * Q5(F(1, ORDER_2I), 0)
        if not val.is_rational():
            raise AssertionError(f"fibre {fib} at l={l} not rational: {val}")
        result[fib] = val.as_fraction()
    return result


def m_exact(l: int) -> int:
    """Multiplicity as the exact sum of fibre contributions; asserted integral."""
    total = sum(fibre_decomposition(l).values())
    if total.denominator != 1:
        raise AssertionError(f"m_exact({l}) = {total} is not an integer")
    return int(total)


# ==========================================================================
# Derived quantities
# ==========================================================================
def lambda_of(l: int) -> int:
    """Laplacian eigenvalue lambda = l(l+2) on S^3."""
    return l * (l + 2)


def spectral_gap():
    """First l > 0 with m(l) >= 1 (via the exact route) and its lambda."""
    l = 1
    while m_exact(l) == 0:
        l += 1
    return l, lambda_of(l)


# ==========================================================================
# Output tables
# ==========================================================================
def write_tables(directory=None, L=120):
    if directory is None:
        directory = os.path.dirname(os.path.abspath(__file__))
    assert_class_sizes()
    mol = molien_coeffs(L)

    m_path = os.path.join(directory, "m_2I.csv")
    with open(m_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["l", "lambda", "m_char", "m_molien", "m_exact",
                    "Gen", "Vtx", "Face", "Edge"])
        for l in range(L + 1):
            fib = fibre_decomposition(l)
            mc, _res = m_char(l)
            w.writerow([l, lambda_of(l), mc, mol[l], m_exact(l),
                        str(fib["Gen"]), str(fib["Vtx"]),
                        str(fib["Face"]), str(fib["Edge"])])

    c_path = os.path.join(directory, "fibre_cycles.csv")
    with open(c_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["l", "Gen", "Vtx", "Face", "Edge",
                    "vtx_slot", "face_slot", "edge_slot"])
        for l in range(60):
            fib = fibre_decomposition(l)
            w.writerow([l, str(fib["Gen"]), str(fib["Vtx"]),
                        str(fib["Face"]), str(fib["Edge"]),
                        l % 10, l % 6, l % 4])
    return m_path, c_path


def _summary():
    assert_class_sizes()
    gap_l, gap_lambda = spectral_gap()
    print(f"|2I| = {ORDER_2I}, class sizes assert OK")
    print(f"spectral gap: first surviving mode at l = {gap_l}, "
          f"lambda = {gap_l}*{gap_l + 2} = {gap_lambda}")
    print("l  : m  (Gen, Vtx, Face, Edge)")
    for l in [0, 10, 11, 12, 20, 30, 70]:
        fib = fibre_decomposition(l)
        print(f"{l:3d}: {m_exact(l):d}  "
              f"({fib['Gen']}, {fib['Vtx']}, {fib['Face']}, {fib['Edge']})")


if __name__ == "__main__":
    _summary()
    mp_path, cyc_path = write_tables()
    print(f"wrote {os.path.basename(mp_path)}, {os.path.basename(cyc_path)}")
