#!/usr/bin/env python3
"""Paper 174 (The Icosahedral Spectrum) -- Appendix A reproduction.

Target quantity
---------------
    m_k = dim( Sym^k(C^2) )^{2I}

the multiplicity of the trivial representation of the binary icosahedral
group 2I (order 120) inside the degree-k SU(2) harmonics, for k = 0 .. 120.

Four independent methods are implemented and cross-checked at every k:

  1. Character sum   (paper eq. 2)  -- average of the SU(2) character
     chi_k over the nine 2I conjugacy classes.  Uses mpmath sin/cos, so it
     is a floating-point route: OBSERVED, rounded to the nearest integer.
  2. Reduced monomial basis (paper eq. 6) -- count of lattice points
     (a,b,c) with 12a + 20b + 30c = k, a,b >= 0, c in {0,1}.  Exact
     integer arithmetic: DERIVED.
  3. Klein complete-intersection Hilbert series (paper eq. 5) -- power
     series coefficients of (1 - t^60) / ((1-t^12)(1-t^20)(1-t^30)).
     Exact integer arithmetic: DERIVED.
  4. Kostant form -- power series coefficients of
     (1 + t^30) / ((1-t^12)(1-t^20)).  Exact integer arithmetic: DERIVED.

The three exact methods (2,3,4) agreeing is a DERIVED statement; method 1
matching them after rounding corroborates the character-theoretic route.

Standalone:  python icosahedral_spectrum.py

Writes icosahedral_spectrum_results.md alongside itself and exits non-zero
on ANY acceptance failure (the reproduction is worthless if it disagrees
with the paper and silently passes).
"""

import os
import sys

import mpmath as mp

mp.mp.dps = 50  # 50 decimal digits: character sum is trivially within 1e-10

K_MAX = 120

# ---------------------------------------------------------------------------
# The nine conjugacy classes of the binary icosahedral group 2I.
# Each entry: (class label, class size, half-angle alpha).
# The SU(2) element is a rotation whose eigenvalues are exp(+/- i*alpha);
# the character of Sym^k(C^2) on it is sin((k+1)alpha)/sin(alpha).
# ---------------------------------------------------------------------------
CLASSES = [
    ("{1}",    1, mp.mpf(0)),
    ("{-1}",   1, mp.pi),
    ("C5",    12, mp.pi / 5),
    ("C5^2",  12, 2 * mp.pi / 5),
    ("-C5^2", 12, 3 * mp.pi / 5),
    ("-C5",   12, 4 * mp.pi / 5),
    ("C3",    20, mp.pi / 3),
    ("-C3",   20, 2 * mp.pi / 3),
    ("C2",    30, mp.pi / 2),
]

GROUP_ORDER = 120


def _chi(k, alpha):
    """SU(2) character of Sym^k(C^2): sin((k+1)alpha)/sin(alpha).

    Limits: alpha -> 0 gives k+1; alpha -> pi gives (-1)^k (k+1).
    """
    if alpha == 0:
        return mp.mpf(k + 1)
    if alpha == mp.pi:
        return mp.mpf((-1) ** k * (k + 1))
    return mp.sin((k + 1) * alpha) / mp.sin(alpha)


def method1_character_sum(k_max=K_MAX):
    """m_k via (1/|G|) sum_classes size * chi_k(alpha). Returns (ints, max_dev)."""
    assert sum(size for _, size, _ in CLASSES) == GROUP_ORDER, "class sizes must sum to 120"
    out, max_dev = [], mp.mpf(0)
    for k in range(k_max + 1):
        s = mp.mpf(0)
        for _, size, alpha in CLASSES:
            s += size * _chi(k, alpha)
        val = s / GROUP_ORDER
        nearest = int(mp.nint(val))
        max_dev = max(max_dev, abs(val - nearest))
        out.append(nearest)
    return out, max_dev


def method2_reduced_monomial(k_max=K_MAX):
    """m_k = #{(a,b,c): 12a+20b+30c=k, a,b>=0, c in {0,1}}. Exact."""
    out = []
    for k in range(k_max + 1):
        count = 0
        for c in (0, 1):
            rem_c = k - 30 * c
            if rem_c < 0:
                continue
            b = 0
            while 20 * b <= rem_c:
                rem_b = rem_c - 20 * b
                if rem_b % 12 == 0:
                    count += 1
                b += 1
        out.append(count)
    return out


def _series_geom(d, k_max):
    """Power series of 1/(1 - t^d): 1 at multiples of d, else 0."""
    s = [0] * (k_max + 1)
    for k in range(0, k_max + 1, d):
        s[k] = 1
    return s


def _series_mul(a, b, k_max):
    """Truncated convolution of two coefficient lists (exact integers)."""
    out = [0] * (k_max + 1)
    for i, ai in enumerate(a):
        if ai == 0 or i > k_max:
            continue
        for j, bj in enumerate(b):
            if i + j > k_max:
                break
            if bj:
                out[i + j] += ai * bj
    return out


def method3_klein_hilbert(k_max=K_MAX):
    """Coefficients of (1 - t^60) / ((1-t^12)(1-t^20)(1-t^30)). Exact."""
    denom = _series_mul(
        _series_mul(_series_geom(12, k_max), _series_geom(20, k_max), k_max),
        _series_geom(30, k_max),
        k_max,
    )
    numer = [0] * (k_max + 1)
    numer[0] = 1
    if 60 <= k_max:
        numer[60] = -1
    return _series_mul(numer, denom, k_max)


def method4_kostant(k_max=K_MAX):
    """Coefficients of (1 + t^30) / ((1-t^12)(1-t^20)). Exact."""
    denom = _series_mul(_series_geom(12, k_max), _series_geom(20, k_max), k_max)
    numer = [0] * (k_max + 1)
    numer[0] = 1
    if 30 <= k_max:
        numer[30] = 1
    return _series_mul(numer, denom, k_max)


def frobenius_number(generators, bound=1000):
    """Largest non-negative integer not representable as a non-negative
    integer combination of `generators` (Frobenius number of the numerical
    semigroup). Assumes gcd(generators) == 1 so it is finite."""
    reachable = [False] * (bound + 1)
    reachable[0] = True
    for n in range(1, bound + 1):
        for g in generators:
            if n >= g and reachable[n - g]:
                reachable[n] = True
                break
    last = max(n for n in range(bound + 1) if not reachable[n])
    return last


# ---------------------------------------------------------------------------
# Acceptance checks -- every one must pass or the script exits non-zero.
# ---------------------------------------------------------------------------
def run():
    lines = []

    def emit(msg):
        print(msg)
        lines.append(msg)

    emit("Paper 174 -- The Icosahedral Spectrum: m_k reproduction (k = 0..120)")
    emit("=" * 70)

    m1, max_dev = method1_character_sum()
    m2 = method2_reduced_monomial()
    m3 = method3_klein_hilbert()
    m4 = method4_kostant()

    failures = []

    def check(cond, label):
        status = "PASS" if cond else "FAIL"
        emit(f"[{status}] {label}")
        if not cond:
            failures.append(label)

    # --- four-method cross-check ---------------------------------------
    mismatches = [k for k in range(K_MAX + 1)
                  if not (m1[k] == m2[k] == m3[k] == m4[k])]
    check(len(mismatches) == 0,
          f"Four-method agreement k=0..120 (mismatches: {mismatches})")

    # --- character-sum float accuracy ----------------------------------
    check(max_dev < mp.mpf("1e-10"),
          f"Character sum within 1e-10 of integers (max deviation {mp.nstr(max_dev, 4)})")

    # --- series identity: method 3 == method 4 -------------------------
    check(m3 == m4,
          "Series identity: Klein (1-t^60)/(1-t^30) == Kostant (1+t^30) for all k<=120")

    m = m2  # canonical (exact) sequence for the remaining checks

    # --- lambda_1 = 168 at k=12 ----------------------------------------
    nonzero_k = [k for k in range(1, K_MAX + 1) if m[k] != 0]
    first = nonzero_k[0] if nonzero_k else None
    check(first == 12, f"First non-trivial survivor at k=12 (got k={first})")
    lam1 = first * (first + 2) if first is not None else None
    check(lam1 == 168, f"lambda_1 = k(k+2) = 12*14 = 168 (got {lam1})")

    # --- all odd k killed ----------------------------------------------
    odd_alive = [k for k in range(1, K_MAX + 1, 2) if m[k] != 0]
    check(len(odd_alive) == 0, f"All odd k killed (alive odd k: {odd_alive})")

    # --- exactly 15 killed even modes ----------------------------------
    killed_even = [k for k in range(2, K_MAX + 1, 2) if m[k] == 0]
    expected_killed = [2, 4, 6, 8, 10, 14, 16, 18, 22, 26, 28, 34, 38, 46, 58]
    check(killed_even == expected_killed,
          f"Exactly 15 killed even modes == expected (got {killed_even})")
    check(len(killed_even) == 15 and max(killed_even) == 58,
          f"15 killed even modes, largest = 58 (got count={len(killed_even)}, "
          f"max={max(killed_even) if killed_even else None})")

    # --- Frobenius(<6,10,15>) = 29  =>  58 largest even non-member ------
    frob = frobenius_number([6, 10, 15])
    check(frob == 29, f"Frobenius(<6,10,15>) = 29 (got {frob})")
    check(2 * frob == 58, f"2 * Frobenius = 58 largest even non-member (got {2*frob})")
    above = [k for k in range(60, K_MAX + 1, 2) if m[k] == 0]
    check(len(above) == 0, f"Every even k > 58 survives (dead even k>58: {above})")

    # --- spot table values ---------------------------------------------
    spot = {12: 1, 20: 1, 24: 1, 60: 2}
    for k, want in spot.items():
        check(m[k] == want, f"m_{k} = {want} (got {m[k]})")

    # --- the m_k table -------------------------------------------------
    emit("")
    emit("m_k for k = 0..120 (survivors only shown with their Laplacian eigenvalue k(k+2)):")
    survivors = [(k, m[k], k * (k + 2)) for k in range(K_MAX + 1) if m[k] != 0]
    for k, mk, lam in survivors:
        emit(f"  k={k:3d}  m_k={mk}  lambda=k(k+2)={lam}")
    emit(f"Total survivors in 0..120: {len(survivors)}")

    emit("")
    emit("=" * 70)
    if failures:
        emit(f"RESULT: FAIL -- {len(failures)} acceptance target(s) failed:")
        for f in failures:
            emit(f"    - {f}")
    else:
        emit("RESULT: PASS -- all acceptance targets reproduced.")

    # write results file next to this script
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "icosahedral_spectrum_results.md"), "w",
              encoding="utf-8") as fh:
        fh.write("# Paper 174 -- icosahedral spectrum: reproduction results\n\n")
        fh.write("```\n")
        fh.write("\n".join(lines))
        fh.write("\n```\n")

    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(run())
