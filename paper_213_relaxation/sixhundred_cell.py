"""
sixhundred_cell.py  —  the 2I arm substrate for Paper 213 Stage 2.

The 600-cell {3,3,5}: its 120 vertices ARE the unit icosians (the binary icosahedral
group 2I on S^3), its 600 cells are regular tetrahedra. It is the 2I-symmetric
triangulation of the round S^3 (2I IMPORTED via the symmetry). We build it as the convex
hull of the 120 icosian points in R^4 and carry the ANALYTIC round metric (edge lengths =
chord distances on the unit S^3) — we do NOT relax it (Stage 1 showed generic relaxation
glasses; the chirality question is independent).

Verify: V=120, E=720, F=1200, C=600, chi=0; 2I-symmetric (vertex/edge-transitive,
5 tets/edge). Sanity: the round scalar spectrum reproduces k(k+2).

Standalone:  python sixhundred_cell.py
"""

import itertools
import numpy as np
from geometry import Complex, verify_sphere, symmetry_probe

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
    """120 unit icosians (2I) as points on S^3 in R^4 (verified construction, Arm 2)."""
    elems = []
    for pos in range(4):                                   # 8 of (+-1,0,0,0)
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


def build_600cell():
    """Return (Complex, points). 600-cell = convex hull of the 120 icosians."""
    from scipy.spatial import ConvexHull
    pts = build_icosians()
    hull = ConvexHull(pts)
    tets = [tuple(s) for s in hull.simplices]
    cx = Complex(120, tets)
    cx.hull_points = pts
    return cx, pts


def round_metric(cx):
    P = cx.hull_points
    return np.array([np.linalg.norm(P[i] - P[j]) for (i, j) in cx.edges])


if __name__ == "__main__":
    cx, pts = build_600cell()
    print(f"600-cell: V={cx.N} E={len(cx.edges)} F={len(cx.triangles)} C={len(cx.tets)}")
    print(f"  expected: V=120 E=720 F=1200 C=600")
    counts_ok = (cx.N == 120 and len(cx.edges) == 720 and
                 len(cx.triangles) == 1200 and len(cx.tets) == 600)
    print(f"  counts match: {counts_ok}")
    chk = verify_sphere(cx, max_link_checks=120)
    print(f"  chi={chk['euler_char']} closed={chk['closed_no_boundary']} links_ok={chk['links_ok']}")
    sym = symmetry_probe(cx)
    print(f"  symmetry: vertex_deg_unique={sym['vertex_deg_unique']} "
          f"edge_tets_unique={sym['edge_tets_unique']} "
          f"mean_tets/edge={sym['mean_tets_per_edge']:.2f} "
          f"looks_2I_regular={sym['looks_2I_regular']}")
    ell = round_metric(cx)
    print(f"  round metric: edge length min={ell.min():.4f} max={ell.max():.4f} "
          f"(regular 600-cell: all equal)")