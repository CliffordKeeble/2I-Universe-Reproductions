"""
spectrum.py -- 2I-invariant spectrum of the Laplacian on S^3 / 2I
(Poincare homology sphere), reproducing Paper 101 v3.0 section 4.3.

Independent derivation. The binary icosahedral group 2I < SU(2) is built
explicitly as the 120 unit icosian quaternions. A degree-l harmonic on S^3
has Laplace eigenvalue lambda_l = l(l+2) and lives in the (l+1)^2-dim space
Sym^l(C^2) (x) Sym^l(C^2)-bar under SU(2)_L x SU(2)_R. The number of modes
surviving the 2I quotient at level l is the multiplicity of the trivial rep
of 2I in Sym^l(C^2), i.e. the dimension of 2I-invariant degree-l binary forms.

That multiplicity a_l is read off the Molien series of 2I acting on C^2:

    M(t) = (1/|2I|) * sum_{g in 2I} 1 / det(1 - t*g)
         = (1/120)  * sum_{g in 2I} 1 / (1 - 2*Re(g)*t + t^2)
         = sum_{l>=0} a_l t^l ,    with  a_l = (1/120) sum_g U_l(Re g),

where U_l is the Chebyshev polynomial of the 2nd kind (the eigenvalues of a
unit quaternion g=(w,x,y,z) as an SU(2) element are e^{+-i theta}, w = cos theta,
so the t-series of 1/(1-2w t + t^2) has coefficients U_l(w)).

A level l "survives" iff a_l >= 1. This script prints the surviving levels up
to l = 120 and their eigenvalues, and asserts the two anchor modes published
in section 4.3:  l=12 -> lambda=168,  l=42 -> lambda=1848.

Run: python spectrum.py
"""

import numpy as np
from itertools import permutations, product

PHI = (1.0 + np.sqrt(5.0)) / 2.0
LMAX = 120


def _parity_even(perm):
    """True if perm (a tuple of 0..n-1) is an even permutation."""
    inv = 0
    n = len(perm)
    for i in range(n):
        for j in range(i + 1, n):
            if perm[i] > perm[j]:
                inv += 1
    return inv % 2 == 0


def build_2I():
    """Return the 120 unit icosian quaternions as (w,x,y,z) tuples."""
    q = set()

    def add(vec):
        q.add(tuple(round(c, 9) for c in vec))

    # 8 unit quaternions: (+-1,0,0,0) and coordinate permutations
    for i in range(4):
        for s in (1.0, -1.0):
            v = [0.0, 0.0, 0.0, 0.0]
            v[i] = s
            add(v)

    # 16 quaternions: (+-1/2, +-1/2, +-1/2, +-1/2)
    for signs in product((0.5, -0.5), repeat=4):
        add(list(signs))

    # 96 quaternions: even permutations of (0, +-1/2, +-1/(2phi), +-phi/2)
    mags = [0.0, 0.5, 1.0 / (2.0 * PHI), PHI / 2.0]
    for perm in permutations(range(4)):
        if not _parity_even(perm):
            continue
        arranged = [mags[perm[i]] for i in range(4)]
        nz = [i for i in range(4) if arranged[i] != 0.0]
        for signs in product((1.0, -1.0), repeat=len(nz)):
            v = list(arranged)
            for k, pos in enumerate(nz):
                v[pos] = arranged[pos] * signs[k]
            add(v)

    return [np.array(v) for v in q]


def molien_multiplicities(group, lmax=LMAX):
    """a_l = (1/|G|) sum_g U_l(Re g) for l = 0..lmax, rounded to int."""
    w = np.array([g[0] for g in group])          # real parts = cos(theta)
    n = len(group)
    # Chebyshev-U recurrence: U_0=1, U_1=2w, U_k = 2w U_{k-1} - U_{k-2}
    Uprev = np.ones_like(w)                        # U_0
    Ucur = 2.0 * w                                 # U_1
    a = np.zeros(lmax + 1)
    a[0] = Uprev.sum() / n
    if lmax >= 1:
        a[1] = Ucur.sum() / n
    for l in range(2, lmax + 1):
        Unext = 2.0 * w * Ucur - Uprev
        a[l] = Unext.sum() / n
        Uprev, Ucur = Ucur, Unext
    a_int = np.rint(a).astype(int)
    # sanity: the Molien coefficients must be (near-)integers
    max_dev = np.max(np.abs(a - a_int))
    return a_int, max_dev


def surviving_levels(lmax=LMAX):
    group = build_2I()
    assert len(group) == 120, f"expected 120 icosians, got {len(group)}"
    a, max_dev = molien_multiplicities(group, lmax)
    levels = [(l, l * (l + 2), int(a[l])) for l in range(lmax + 1) if a[l] >= 1]
    return levels, max_dev


if __name__ == "__main__":
    levels, max_dev = surviving_levels(LMAX)
    print(f"Molien coefficient max deviation from integer: {max_dev:.2e}")
    print(f"Surviving 2I levels up to l={LMAX} (l, lambda=l(l+2), multiplicity a_l):")
    for l, lam, mult in levels:
        print(f"  l={l:3d}  lambda={lam:6d}  a_l={mult}")
    # exclude l=0 (constant / zero mode) when counting physical modes
    nz = [x for x in levels if x[0] > 0]
    print(f"\nSurviving levels (l>0) up to l={LMAX}: {len(nz)}")

    anchors = {l: lam for l, lam, _ in levels}
    assert anchors.get(12) == 168,  "section 4.3 anchor l=12 -> 168 FAILED"
    assert anchors.get(42) == 1848, "section 4.3 anchor l=42 -> 1848 FAILED"
    print("Anchor check: l=12 -> 168 OK, l=42 -> 1848 OK (reproduces section 4.3).")
