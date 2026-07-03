"""
icosian_2I.py  —  the binary icosahedral group 2I as SU(2) target rotations.

The 120 unit icosians (= the 600-cell vertices = the binary icosahedral group 2I)
are the unit quaternions that S^3/2I quotients by. For Paper 214 we need them not as
a triangulation (that was Paper 213) but as concrete SU(2) elements acting on the
Skyrme TARGET S^3 = SU(2).  S^3/2I identifies field values U ~ q U for q in 2I.

This module:
  * builds the 120 unit quaternions (verified construction, same as Paper 213 Arm 2),
  * verifies they form a group of order 120 under quaternion multiplication
    (closure -> finite group; |2I| = 120),
  * exposes each as a 2x2 SU(2) matrix for left-action on the target.

Self-contained (brief: no cross-directory deps). Standalone:  python icosian_2I.py
"""

import itertools
import numpy as np

PHI = (1.0 + np.sqrt(5.0)) / 2.0
INV_PHI = PHI - 1.0


def _even_perms4():
    out = []
    for p in itertools.permutations(range(4)):
        inv = sum(1 for i in range(4) for j in range(i + 1, 4) if p[i] > p[j])
        if inv % 2 == 0:
            out.append(p)
    return out


def build_icosians():
    """120 unit icosians (2I) as quaternions (w,x,y,z) on S^3 in R^4."""
    elems = []
    for pos in range(4):                                     # 8 of (+-1,0,0,0)
        for s in (+1.0, -1.0):
            v = np.zeros(4); v[pos] = s; elems.append(v)
    for signs in itertools.product((+1.0, -1.0), repeat=4):  # 16 of 1/2(+-1,+-1,+-1,+-1)
        elems.append(0.5 * np.array(signs))
    mags = (0.0, 1.0, PHI, INV_PHI)                          # 96 even perms of 1/2(0,+-1,+-phi,+-1/phi)
    for perm in _even_perms4():
        for nz in itertools.product((+1.0, -1.0), repeat=3):
            c = np.zeros(4); it = iter(nz)
            for slot, mi in enumerate(perm):
                m = mags[mi]
                c[slot] = 0.0 if m == 0.0 else next(it) * m
            elems.append(0.5 * c)
    arr = np.array(elems)
    assert arr.shape == (120, 4)
    assert np.allclose(np.linalg.norm(arr, axis=1), 1.0)
    return arr


def quat_mul(a, b):
    """Hamilton product of quaternions a,b given as (w,x,y,z) arrays."""
    w1, x1, y1, z1 = a
    w2, x2, y2, z2 = b
    return np.array([
        w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
        w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
        w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
        w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
    ])


def verify_group(elems, tol=1e-9):
    """Confirm the 120 quaternions are closed under multiplication (a group of order 120).

    Returns dict with closure flag, identity present, and inverse-closure.
    """
    n = len(elems)
    # Build a lookup by rounding to match products back to the element set.
    keys = {tuple(np.round(e, 6)): i for i, e in enumerate(elems)}

    def find(q):
        k = tuple(np.round(q, 6))
        if k in keys:
            return keys[k]
        d = np.linalg.norm(elems - q, axis=1)   # tolerant fallback to nearest
        j = int(np.argmin(d))
        return j if d[j] < tol else -1

    closed = True
    for i in range(n):
        for j in range(n):
            p = quat_mul(elems[i], elems[j])
            if find(p) < 0:
                closed = False
                break
        if not closed:
            break

    has_identity = find(np.array([1.0, 0, 0, 0])) >= 0
    # inverse of unit quaternion = conjugate; check each conjugate is in the set
    inv_closed = all(find(np.array([e[0], -e[1], -e[2], -e[3]])) >= 0 for e in elems)
    return {"order": n, "closed_under_mul": closed,
            "has_identity": has_identity, "inverse_closed": inv_closed}


def quat_to_su2(q):
    """Unit quaternion (w,x,y,z) -> 2x2 SU(2) matrix.  U = w 1 + x(i sx) + y(i sy) + z(i sz)."""
    w, x, y, z = q
    return np.array([[w + 1j * z, y + 1j * x],
                     [-y + 1j * x, w - 1j * z]], dtype=complex)


def icosians_as_su2():
    """Return (120, 2, 2) array of SU(2) matrices for the 2I elements."""
    Q = build_icosians()
    return np.array([quat_to_su2(q) for q in Q]), Q


if __name__ == "__main__":
    Q = build_icosians()
    info = verify_group(Q)
    print("=== binary icosahedral group 2I ===")
    print(f"  order (num elements)     = {info['order']}        (expected 120)")
    print(f"  closed under quat mult   = {info['closed_under_mul']}")
    print(f"  contains identity        = {info['has_identity']}")
    print(f"  closed under inverse     = {info['inverse_closed']}")
    # SU(2) sanity: each matrix unitary with det 1
    U, _ = icosians_as_su2()
    dets = np.array([np.linalg.det(u) for u in U])
    unit = np.array([np.max(np.abs(u.conj().T @ u - np.eye(2))) for u in U])
    print(f"  max |det - 1|            = {np.max(np.abs(dets - 1.0)):.2e}")
    print(f"  max ||U^dag U - 1||      = {np.max(unit):.2e}")
    ok = (info['order'] == 120 and info['closed_under_mul']
          and info['has_identity'] and info['inverse_closed'])
    print(f"  GROUP VERIFIED           = {ok}")
