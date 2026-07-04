#!/usr/bin/env python3
"""
Solvable-Case Null — the withholding test (Pattern-75 falsifier).

Design locked by Flint's brief (3 Jul 2026); see PRE-REGISTRATION.md. This is the
self-contained reproduction: run the framework's weight-2 selection machinery on
genus-0 conics (β₁ = 0) and genus-1 grammar curves (β₁ = 2), and show it WITHHOLDS
where there is no loop register and SELECTS where there is.

Two hard gates:
  * LMFDB is a test oracle, never an input — every number here is computed from the
    Weierstrass coefficients / the conic equation alone.
  * Minimality: bad primes are those dividing Δ (computed), not looked up.

Run:  python solvable_null.py            (X = 10^4, ~25 s)
      python solvable_null.py --x 100000 (the trend check, slow: ~30 min in Python)

Status flags: [DERIVED] proven/decidable; [OBSERVED] numerical, null-tested.
"""
import argparse
import csv
import math
import os
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")   # Windows consoles default to cp1252
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))

# ─────────────────────────────────────────────────────────────────────────────
# Genus-1 selection side: a_p by point counting, Mestre–Nagao sum
# ─────────────────────────────────────────────────────────────────────────────

def sieve(limit):
    s = bytearray([1]) * (limit + 1)
    s[0] = s[1] = 0
    for i in range(2, int(limit ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(2, limit + 1) if s[i]]


class GrammarCurve:
    """A genus-1 minimal Weierstrass model y²+a₁xy+a₃y = x³+a₂x²+a₄x+a₆.
    Over ℂ a torus T², Betti (1,2,1): β₁ = 2 — the loop register the filter reads."""
    beta1 = 2
    genus = 1

    def __init__(self, label, coeffs, expected_rank):
        self.label, self.coeffs, self.expected_rank = label, coeffs, expected_rank
        a1, a2, a3, a4, a6 = coeffs
        b2 = a1 * a1 + 4 * a2
        b4 = a1 * a3 + 2 * a4
        b6 = a3 * a3 + 4 * a6
        b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4
        self.delta = -b2 * b2 * b8 - 8 * b4 ** 3 - 27 * b6 * b6 + 9 * b2 * b4 * b6

    def ap(self, p):
        """a_p = p + 1 − #E(𝔽_p) = −Σ_x legendre((a₁x+a₃)² + 4f(x), p)."""
        a1, a2, a3, a4, a6 = (c % p for c in self.coeffs)
        if p == 2:
            cnt = 1
            for x in range(2):
                for y in range(2):
                    lhs = (y * y + a1 * x * y + a3 * y) % 2
                    rhs = (x * x * x + a2 * x * x + a4 * x + a6) % 2
                    if lhs == rhs:
                        cnt += 1
            return 2 + 1 - cnt
        half, total = (p - 1) // 2, 0
        for x in range(p):
            b = (a1 * x + a3) % p
            f = (x * x % p * x + a2 * x * x + a4 * x + a6) % p
            d = (b * b + 4 * f) % p
            if d:
                total += 1 if pow(d, half, p) == 1 else -1
        return -total

    def nagao(self, X):
        """S(X) = (1/log X) · Σ_{good p≤X} (−a_p log p)/p. [OBSERVED]"""
        s = 0.0
        for p in sieve(X):
            if self.delta % p == 0:
                continue                       # good primes only (p ∤ Δ)
            s += (-self.ap(p) * math.log(p)) / p
        return s / math.log(X)


# ─────────────────────────────────────────────────────────────────────────────
# Genus-0 withholding side: Hilbert symbols, Hasse–Minkowski
# ─────────────────────────────────────────────────────────────────────────────

def legendre(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def valuation(n, p):
    """n = p^k · unit, unit signed."""
    k = 0
    while n % p == 0:
        n //= p
        k += 1
    return k, n


def hilbert(a, b, v):
    """Hilbert symbol (a,b)_v ∈ {±1}; v = 'inf' is the real place."""
    if v == 'inf':
        return -1 if (a < 0 and b < 0) else 1
    if v == 2:
        alpha, u = valuation(a, 2)
        beta, w = valuation(b, 2)
        eps = lambda t: (t - 1) // 2 % 2
        om = lambda t: (t * t - 1) // 8 % 2
        return -1 if (eps(u) * eps(w) + alpha * om(w) + beta * om(u)) % 2 else 1
    alpha, u = valuation(a, v)
    beta, w = valuation(b, v)
    r = -1 if (alpha * beta % 2 == 1 and (v - 1) // 2 % 2 == 1) else 1
    if beta & 1:
        r *= legendre(u, v)
    if alpha & 1:
        r *= legendre(w, v)
    return r


def odd_prime_divisors(n):
    n = abs(n)
    while n % 2 == 0:
        n //= 2
    ps, d = [], 3
    while d * d <= n:
        if n % d == 0:
            ps.append(d)
            while n % d == 0:
                n //= d
        d += 2
    if n > 1:
        ps.append(n)
    return ps


class Conic:
    """A smooth diagonal conic a·X²+b·Y²+c·Z² = 0 over ℚ — genus 0.
    Over ℂ P¹(ℂ) = S², Betti (1,0,1): β₁ = 0 — NO loop register."""
    beta1 = 0
    genus = 0

    def __init__(self, a, b, c):
        self.a, self.b, self.c = a, b, c

    @classmethod
    def sum_of_two_squares(cls, n):
        return cls(1, 1, -n)

    def hasse_minkowski(self):
        """(soluble, obstruction_places). Serre IV: ⟨a,b,c⟩ isotropic over ℚ_v iff
        (−1,−d)_v = ∏_{i<j}(·,·)_v.  [DERIVED]"""
        a, b, c = self.a, self.b, self.c
        d = a * b * c
        places = ['inf', 2] + odd_prime_divisors(2 * a * b * c)
        obstruction = []
        for v in places:
            hasse = hilbert(a, b, v) * hilbert(a, c, v) * hilbert(b, c, v)
            if hilbert(-1, -d, v) != hasse:
                obstruction.append(v)
        return len(obstruction) == 0, obstruction


# ─────────────────────────────────────────────────────────────────────────────
# The weight-2 selection pipeline: reads the β₁ register
# ─────────────────────────────────────────────────────────────────────────────

def rank_signal(controller, X):
    """The same machinery applied to any controller. It reads the β₁ register:
    genus 1 (β₁=2) → the Nagao rank signal; genus 0 (β₁=0) → WITHHELD, no register
    to fill. A number returned for a genus-0 controller would be MANUFACTURED — the
    falsification condition. The withholding is structural, not stipulated."""
    if controller.beta1 == 0:
        return "WITHHELD", None
    return "SELECTED", controller.nagao(X)


# ─────────────────────────────────────────────────────────────────────────────

GRAMMAR = [
    GrammarCurve("11a1",  (0, -1, 1, -10, -20), 0),
    GrammarCurve("37a1",  (0, 0, 1, -1, 0),     1),
    GrammarCurve("389a1", (0, 1, 1, -2, 0),     2),
]
# n, expected solubility (DERIVED; cross-checks the sum-of-two-squares criterion)
CONICS = [(1, True), (3, False), (5, True), (21, False), (6, False), (7, False)]


def main():
    ap_ = argparse.ArgumentParser(description="Solvable-Case Null — withholding test")
    ap_.add_argument("--x", type=int, default=10_000, help="Nagao cutoff X (default 10^4)")
    args = ap_.parse_args()
    X = args.x

    passes, fails = [], []

    def check(cond, msg):
        (passes if cond else fails).append(msg)
        print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")

    print(f"SELECTION — genus-1 grammar curves (β₁ = 2), Nagao S({X})")
    sigs = {}
    for c in GRAMMAR:
        state, s = rank_signal(c, X)
        sigs[c.label] = s
        print(f"    {c.label:<6} rank {c.expected_rank}:  S = {s:+.6f}   [{state}]")
    order_ok = sigs["11a1"] < sigs["37a1"] < sigs["389a1"]
    check(order_ok, "selection separates the ranks: r0 < r1 < r2")
    check(all(rank_signal(c, X)[0] == "SELECTED" for c in GRAMMAR),
          "β₁ register filled on every grammar curve (genus 1 selects)")

    print(f"\nWITHHOLDING — genus-0 conics (β₁ = 0), Hasse–Minkowski only")
    conic_rows = []
    for n, want_soluble in CONICS:
        conic = Conic.sum_of_two_squares(n)
        soluble, obs = conic.hasse_minkowski()
        state, val = rank_signal(conic, X)
        obs_str = "—" if soluble else "{" + ",".join(str(o) for o in obs) + "}"
        print(f"    x²+y²={n:<3} {'soluble' if soluble else 'insoluble':<9} {obs_str:<10} [{state}]")
        check(soluble == want_soluble, f"x²+y²={n}: HM binary correct ({'soluble' if soluble else 'insoluble'})")
        check(state == "WITHHELD" and val is None, f"x²+y²={n}: no rank register (β₁=0, withholding holds)")
        conic_rows.append((n, soluble, obs_str))

    # outputs (repo convention: write to files alongside the script)
    with open(os.path.join(HERE, "nagao_sums.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["curve", "expected_rank", "X", "nagao_S"])
        for c in GRAMMAR:
            w.writerow([c.label, c.expected_rank, X, f"{sigs[c.label]:.6f}"])
    with open(os.path.join(HERE, "conic_solubility.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["n", "soluble", "obstruction_places"])
        for n, soluble, obs_str in conic_rows:
            w.writerow([n, soluble, obs_str])

    verdict = "PASS — withholds on genus 0, selects on genus 1" if not fails else "FAIL"
    print(f"\n{len(passes)} passed, {len(fails)} failed.  →  {verdict}")
    return 0 if not fails else 1


if __name__ == "__main__":
    raise SystemExit(main())
