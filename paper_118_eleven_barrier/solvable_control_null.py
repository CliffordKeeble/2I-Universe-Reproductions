"""
Paper 118 Part 2 — the solvable-control null (pre-registered).

Pre-registration: PRE-REGISTRATION.md, committed 4d74af8 (Gate G1) before this ran.
Substrate: Part 1 machinery, commit 96a290f (AT1-AT9 green).

Runs the barrier predicate across every family of finite Gamma < SU(2) and reports,
per Gamma, whether the top-fibre perfect return kills the mode. Tests the pre-registered
H_unique (barrier only for 2I) vs H_arith (barrier for some solvable Gamma). Findings only.

Two independent routes per multiplicity (Gate G3):
  (A) mpmath >= 50 dps character sum, integrality assertion |m - round(m)| < 1e-30;
  (B) exact integer -- Molien series {2T,2O,2I}, weight-count {Z_n},
      m_{2D_n}(l) = (m_{Z_2n}(l) + chi_l(pi/2)) / 2 {2D_n}.

Run:  python solvable_control_null.py
"""

from __future__ import annotations

import csv
import os
from fractions import Fraction as F

from mpmath import mp, mpf, nint, fabs

import eleven_barrier as eb   # chi_mp (Part 1, commit 96a290f)

mp.dps = 60

PREREG_COMMIT = "4d74af8"
SUBSTRATE_COMMIT = "96a290f"


# ==========================================================================
# Population -- class data grouped by theta = p*pi/q, as (p, q, size)
# ==========================================================================
def classes_cyclic(n):
    """Z_n < SU(2): generator diag(z, z^-1), z = e^{2 pi i / n}.
    n classes at theta_k = 2 pi k / n = (2k) pi / n, size 1 each."""
    return [(2 * k, n, 1) for k in range(n)], n


def classes_dihedral(n):
    """Binary dihedral 2D_n (dicyclic), order 4n:
       I (theta 0), -I (theta pi), n-1 cyclic classes theta = k pi / n (size 2),
       and two order-4 'reflection' classes at theta = pi/2 (size n each -> 2n)."""
    cls = [(0, 1, 1), (1, 1, 1)]
    cls += [(k, n, 2) for k in range(1, n)]
    cls += [(1, 2, 2 * n)]           # both period-4 fibres, theta = pi/2
    return cls, 4 * n


CLASSES_2T = ([(0, 1, 1), (1, 1, 1), (1, 2, 6), (1, 3, 8), (2, 3, 8)], 24)
CLASSES_2O = ([(0, 1, 1), (1, 1, 1), (1, 4, 6), (3, 4, 6),
               (1, 3, 8), (2, 3, 8), (1, 2, 18)], 48)
CLASSES_2I = ([(0, 1, 1), (1, 1, 1),
               (1, 5, 12), (2, 5, 12), (3, 5, 12), (4, 5, 12),
               (1, 3, 20), (2, 3, 20), (1, 2, 30)], 120)


# ==========================================================================
# Route A -- mpmath character sum (general)
# ==========================================================================
def m_char(classes, order, l):
    s = mpf(0)
    for p, q, size in classes:
        s += size * eb.chi_mp(l, p, q)
    m = s / order
    mr = int(nint(m))
    return mr, fabs(m - mr)


def chi_val(p, q, l):
    """Numeric character value chi_l(theta), theta = p pi / q (for slot reports)."""
    return eb.chi_mp(l, p, q)


# ==========================================================================
# Route B -- exact integer
# ==========================================================================
def m_cyclic_exact(order, l):
    """dim Sym^l(C^2)^{Z_order}: weights w = 2i - l with order | w."""
    return sum(1 for i in range(l + 1) if (2 * i - l) % order == 0)


def chi_half_pi(l):
    """chi_l(pi/2) = sin((l+1) pi/2), exact in {0, +1, -1}."""
    if l % 2 == 1:
        return 0
    return 1 if (l // 2) % 2 == 0 else -1


def m_dihedral_exact(n, l):
    """m_{2D_n}(l) = (m_{Z_2n}(l) + chi_l(pi/2)) / 2."""
    num = m_cyclic_exact(2 * n, l) + chi_half_pi(l)
    assert num % 2 == 0, f"2D_{n}: non-integral at l={l} (num={num})"
    return num // 2


def molien(num_deg, dens, L):
    """Coefficients of (1 + t^num_deg) / prod(1 - t^d) up to t^L, exact ints."""
    s = [0] * (L + 1)
    s[0] = 1
    for d in dens:
        for k in range(d, L + 1):
            s[k] += s[k - d]
    return [s[k] + (s[k - num_deg] if k >= num_deg else 0) for k in range(L + 1)]


MOLIEN = {  # (num_deg, dens)
    "2T": (12, [6, 8]),
    "2O": (18, [8, 12]),
    "2I": (30, [12, 20]),
}


# ==========================================================================
# Per-group analysis
# ==========================================================================
def analyse(name, kind, param, classes, order, p1):
    # G3 substrate checks
    assert sum(sz for _p, _q, sz in classes) == order, f"{name}: sizes != order"
    L = max(120, 2 * p1)

    # route B array
    if kind == "poly":
        mol = molien(*MOLIEN[name], L)
        route_b = lambda l: mol[l]                       # noqa: E731
    elif kind == "cyclic":
        route_b = lambda l: m_cyclic_exact(param, l)     # noqa: E731
    else:  # dihedral
        route_b = lambda l: m_dihedral_exact(param, l)   # noqa: E731

    m = [0] * (L + 1)
    max_res = mpf(0)
    for l in range(L + 1):
        ma, res = m_char(classes, order, l)
        mb = route_b(l)
        assert ma == mb, f"{name}: route disagreement at l={l} (A={ma}, B={mb})"
        assert res < mpf("1e-30"), f"{name}: integrality residual {res} at l={l}"
        m[l] = ma
        max_res = max(max_res, res)

    assert m[0] == 1, f"{name}: m(0) = {m[0]} != 1"

    l_gap = next(l for l in range(1, L + 1) if m[l] >= 1)
    lam1 = l_gap * (l_gap + 2)
    m_at_p1 = m[p1]
    barrier = (m_at_p1 == 0)

    if l_gap < p1:
        gap_vs_p1 = "gap<p1"
    elif l_gap == p1:
        gap_vs_p1 = "gap==p1"
    else:
        gap_vs_p1 = "gap>p1"

    # fibre-slot (per-class character) values at p1
    slots = []
    for p, q, size in classes:
        c = chi_val(p, q, p1)
        slots.append((p, q, size, c))

    return {
        "name": name, "kind": kind, "order": order, "solvable": name != "2I",
        "L": L, "m": m, "l_gap": l_gap, "lam1": lam1, "p1": p1,
        "m_at_p1": m_at_p1, "barrier": barrier, "gap_vs_p1": gap_vs_p1,
        "max_res": max_res, "slots": slots,
    }


def build_population():
    pop = []
    for n in range(2, 61):
        cls, order = classes_cyclic(n)
        pop.append(("Z%d" % n, "cyclic", n, cls, order, 2 * n))
    for n in range(2, 31):
        cls, order = classes_dihedral(n)
        pop.append(("2D%d" % n, "dihedral", n, cls, order, 2 * n))
    pop.append(("2T", "poly", None, CLASSES_2T[0], CLASSES_2T[1], 6))
    pop.append(("2O", "poly", None, CLASSES_2O[0], CLASSES_2O[1], 8))
    pop.append(("2I", "poly", None, CLASSES_2I[0], CLASSES_2I[1], 10))
    return pop


# ==========================================================================
# Outputs
# ==========================================================================
def run(directory=None):
    if directory is None:
        directory = os.path.dirname(os.path.abspath(__file__))

    results = [analyse(*g) for g in build_population()]
    by_name = {r["name"]: r for r in results}

    # Machinery-validation anchors -- classical, not hypotheses.
    anchors = {"2T": 6, "2O": 8, "2I": 12}
    anchor_report = {}
    for nm, expect in anchors.items():
        got = by_name[nm]["l_gap"]
        anchor_report[nm] = (got, expect, got == expect)
        assert got == expect, (
            f"ANCHOR FAIL: first invariant of {nm} at l={got}, expected {expect} "
            f"-- machinery is wrong, stopping.")

    # barrier_summary
    sum_path = os.path.join(directory, "barrier_summary.csv")
    with open(sum_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Gamma", "order", "solvable", "l_gap", "lambda1",
                    "lambda1_vs_order", "p1", "m_at_p1", "barrier", "gap_vs_p1"])
        for r in results:
            lam_cmp = ("=" if r["lam1"] == r["order"]
                       else (">" if r["lam1"] > r["order"] else "<"))
            w.writerow([r["name"], r["order"], r["solvable"], r["l_gap"], r["lam1"],
                        f"{r['lam1']}{lam_cmp}{r['order']}", r["p1"],
                        r["m_at_p1"], r["barrier"], r["gap_vs_p1"]])

    # full m(l) tables, combined
    m_path = os.path.join(directory, "m_tables.csv")
    with open(m_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Gamma", "l", "m"])
        for r in results:
            for l in range(r["L"] + 1):
                w.writerow([r["name"], l, r["m"][l]])

    # fibre slots at p1
    slot_path = os.path.join(directory, "fibre_slots_at_p1.csv")
    with open(slot_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Gamma", "p1", "theta", "size", "chi_at_p1", "returned_to_1"])
        for r in results:
            for p, q, size, c in r["slots"]:
                theta = "0" if p == 0 else ("pi" if p == q else f"{p}pi/{q}")
                w.writerow([r["name"], r["p1"], theta, size,
                            mp.nstr(c, 8), bool(fabs(c - 1) < mpf("1e-20"))])

    barriers = [r["name"] for r in results if r["barrier"]]
    solvable_barriers = [r["name"] for r in results if r["barrier"] and r["solvable"]]
    landed = ("H_arith" if solvable_barriers
              else ("H_unique" if barriers == ["2I"] else "NEITHER"))

    return results, by_name, anchor_report, barriers, solvable_barriers, landed, \
        (sum_path, m_path, slot_path)


def _print_summary(results, anchor_report, barriers, solvable_barriers, landed):
    print(f"pre-reg {PREREG_COMMIT}, substrate {SUBSTRATE_COMMIT}")
    print("anchors:", {k: f"{v[0]} (want {v[1]}) {'OK' if v[2] else 'FAIL'}"
                       for k, v in anchor_report.items()})
    print(f"population: {len(results)} groups "
          f"({sum(1 for r in results if r['solvable'])} solvable, 1 non-solvable)")
    print(f"barrier holds for: {barriers}")
    print(f"solvable barriers: {solvable_barriers or 'none'}")
    print(f"sweep landed on: {landed}")
    print("\nkey rows (Gamma: order  l_gap lam1  p1  m@p1  barrier):")
    for nm in ["Z5", "Z10", "2D5", "2D15", "2T", "2O", "2I"]:
        r = next(x for x in results if x["name"] == nm)
        print(f"  {nm:5s}: {r['order']:4d}  {r['l_gap']:3d}  {r['lam1']:4d}  "
              f"{r['p1']:3d}  {r['m_at_p1']:3d}   {r['barrier']}")


if __name__ == "__main__":
    res, by_name, anchors, barr, solv_barr, landed, paths = run()
    _print_summary(res, anchors, barr, solv_barr, landed)
    print("wrote:", ", ".join(os.path.basename(p) for p in paths))
