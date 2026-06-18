"""
icosians.py  —  The imported object for Paper 213, ARM 2.

    *** 2I IS IMPORTED. THIS FILE IS THE IMPORT. ***

This module hand-builds the 120 unit icosians (the binary icosahedral group
2I, det +1, golden coordinates) and their 4x4 real left-multiplication
representation in SO(4). Choosing icosians as the holonomy alphabet is the act
of hand-coding 2I into the model. Nothing here "emerges"; the whole point of
Arm 2 is to ask whether a discrete graph carrying THIS imported holonomy can
*reconstruct* the continuum S^3/2I even-sector spectrum. Reconstruction is a
spectral-convergence result, not an ignition.

The build FAILS HARD if the group does not match 2I exactly:
  - 120 distinct unit quaternions,
  - closed under quaternion multiplication,
  - element-order distribution {1:1, 2:1, 3:20, 4:30, 5:24, 6:20, 10:24}
    (= SL(2,5)).

Standalone:  python icosians.py
"""

import itertools
import numpy as np

PHI = (1.0 + np.sqrt(5.0)) / 2.0          # golden ratio
INV_PHI = PHI - 1.0                        # 1/phi = phi - 1
TOL = 1e-9                                  # min separation between icosians >> this

# Target order distribution for 2I = SL(2,5). Build fails if not matched.
TARGET_ORDER_DIST = {1: 1, 2: 1, 3: 20, 4: 30, 5: 24, 6: 20, 10: 24}

IDENTITY = np.array([1.0, 0.0, 0.0, 0.0])


# ---------------------------------------------------------------------------
# Quaternion arithmetic (q = a + b i + c j + d k, stored as [a,b,c,d])
# ---------------------------------------------------------------------------
def qmul(p, q):
    a1, b1, c1, d1 = p
    a2, b2, c2, d2 = q
    return np.array([
        a1*a2 - b1*b2 - c1*c2 - d1*d2,
        a1*b2 + b1*a2 + c1*d2 - d1*c2,
        a1*c2 - b1*d2 + c1*a2 + d1*b2,
        a1*d2 + b1*c2 - c1*b2 + d1*a2,
    ])


def qconj(q):
    a, b, c, d = q
    return np.array([a, -b, -c, -d])


# ---------------------------------------------------------------------------
# Build the 120 unit icosians
# ---------------------------------------------------------------------------
def _even_permutations(n=4):
    """Even permutations of (0..n-1) — the alternating group A4 for n=4."""
    out = []
    for perm in itertools.permutations(range(n)):
        # sign of permutation
        inv = sum(1 for i in range(n) for j in range(i+1, n) if perm[i] > perm[j])
        if inv % 2 == 0:
            out.append(perm)
    return out


def build_icosians():
    """Return (120, 4) array of the unit icosians."""
    elems = []

    # Group 1: 8 of (+-1, 0, 0, 0) type  ->  +-1, +-i, +-j, +-k
    for pos in range(4):
        for s in (+1.0, -1.0):
            v = np.zeros(4)
            v[pos] = s
            elems.append(v)

    # Group 2: 16 of 1/2 (+-1, +-1, +-1, +-1)
    for signs in itertools.product((+1.0, -1.0), repeat=4):
        elems.append(0.5 * np.array(signs))

    # Group 3: 96 of 1/2 * even-permutations of (0, +-1, +-phi, +-1/phi)
    magnitudes = (0.0, 1.0, PHI, INV_PHI)
    for perm in _even_permutations(4):          # 12 even permutations
        for nz_signs in itertools.product((+1.0, -1.0), repeat=3):  # 8 sign sets
            coords = np.zeros(4)
            sign_iter = iter(nz_signs)
            for slot, mag_idx in enumerate(perm):
                mag = magnitudes[mag_idx]
                if mag == 0.0:
                    coords[slot] = 0.0
                else:
                    coords[slot] = next(sign_iter) * mag
            elems.append(0.5 * coords)

    arr = np.array(elems)
    assert arr.shape == (120, 4), f"expected 120 elements, got {arr.shape[0]}"
    return arr


# ---------------------------------------------------------------------------
# Group structure helpers (float matching against the canonical 120)
# ---------------------------------------------------------------------------
def _match_index(q, group):
    """Index of the element of `group` equal to q (within TOL), else -1."""
    if len(group) == 0:
        return -1
    d = np.linalg.norm(group - q[None, :], axis=1)
    j = int(np.argmin(d))
    return j if d[j] < TOL else -1


def element_order(q, group, max_order=64):
    """Smallest m>=1 with q^m == identity (matched within TOL)."""
    acc = q.copy()
    for m in range(1, max_order + 1):
        if np.linalg.norm(acc - IDENTITY) < TOL:
            return m
        acc = qmul(acc, q)
    return -1  # not a finite-order element within bound (should never happen)


def is_closed(group):
    """True iff every product of two elements lands back in the group."""
    for i in range(len(group)):
        for j in range(len(group)):
            if _match_index(qmul(group[i], group[j]), group) < 0:
                return False
    return True


def order_distribution(group):
    dist = {}
    for q in group:
        o = element_order(q, group)
        dist[o] = dist.get(o, 0) + 1
    return dict(sorted(dist.items()))


# ---------------------------------------------------------------------------
# 4x4 real left-multiplication representation (edge holonomies, in SO(4))
# ---------------------------------------------------------------------------
def left_matrix(q):
    """4x4 real matrix of left-multiplication by quaternion q. In SO(4) for |q|=1."""
    a, b, c, d = q
    return np.array([
        [a, -b, -c, -d],
        [b,  a, -d,  c],
        [c,  d,  a, -b],
        [d, -c,  b,  a],
    ])


def all_left_matrices(group):
    return np.array([left_matrix(q) for q in group])


# ---------------------------------------------------------------------------
# A verified 2-generator set {s, t} with <s,t> = all 120
# ---------------------------------------------------------------------------
def generated_subgroup(generators, group):
    """BFS closure of the subgroup generated by `generators`. Returns index set."""
    gens = list(generators)
    frontier = [IDENTITY.copy()]
    seen_idx = set()
    start = _match_index(IDENTITY, group)
    seen_idx.add(start)
    seen_q = [IDENTITY.copy()]
    while frontier:
        nxt = []
        for q in frontier:
            for g in gens:
                p = qmul(q, g)
                idx = _match_index(p, group)
                if idx >= 0 and idx not in seen_idx:
                    seen_idx.add(idx)
                    seen_q.append(p)
                    nxt.append(p)
        frontier = nxt
    return seen_idx


def find_generator_pair(group):
    """
    Find an (order-4, order-{6 or 10}) pair of (2,3,5)-type generating all of 2I.
    Returns (i, j, order_i, order_j). Deterministic scan, no randomness.
    """
    orders = np.array([element_order(q, group) for q in group])
    order4 = [i for i in range(120) if orders[i] == 4]
    order_high = [j for j in range(120) if orders[j] in (6, 10)]
    for i in order4:
        for j in order_high:
            if len(generated_subgroup([group[i], group[j]], group)) == 120:
                return i, j, int(orders[i]), int(orders[j])
    raise RuntimeError("no generating (order-4, order-{6,10}) pair found")


# ---------------------------------------------------------------------------
# Self-test / build verification
# ---------------------------------------------------------------------------
def verify(verbose=True):
    G = build_icosians()

    # 1. unit norm
    norms = np.linalg.norm(G, axis=1)
    assert np.allclose(norms, 1.0, atol=TOL), "non-unit icosian found"

    # 2. distinctness
    n_distinct = 0
    for i in range(120):
        if _match_index(G[i], G[:i]) < 0:
            n_distinct += 1
    assert n_distinct == 120, f"only {n_distinct} distinct elements"

    # 3. closure
    assert is_closed(G), "group is NOT closed under quaternion multiplication"

    # 4. order distribution == SL(2,5)
    dist = order_distribution(G)
    assert dist == TARGET_ORDER_DIST, (
        f"order distribution mismatch:\n  got    {dist}\n  target {TARGET_ORDER_DIST}"
    )

    # 5. left-multiplication matrices lie in SO(4)
    M = all_left_matrices(G)
    for k in range(120):
        assert abs(np.linalg.det(M[k]) - 1.0) < 1e-9, "left matrix not in SO(4) (det != 1)"
        assert np.allclose(M[k] @ M[k].T, np.eye(4), atol=1e-9), "left matrix not orthogonal"

    # 6. a generating 2-element set
    i, j, oi, oj = find_generator_pair(G)

    if verbose:
        print("icosians.py — 2I build VERIFIED (2I is IMPORTED)")
        print(f"  elements        : {len(G)} unit quaternions")
        print(f"  order dist      : {dist}")
        print(f"  matches SL(2,5) : {dist == TARGET_ORDER_DIST}")
        print(f"  closed          : True")
        print(f"  left-mult rep   : 120 matrices in SO(4)")
        print(f"  generator pair  : s=idx{i} (order {oi}), t=idx{j} (order {oj}) -> <s,t> = all 120")
    return G, (i, j)


if __name__ == "__main__":
    verify()
