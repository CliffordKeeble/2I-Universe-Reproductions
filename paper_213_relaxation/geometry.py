"""
geometry.py  —  substrate + Regge machinery for Paper 213 Approach B (relaxation).

We GRANT the closed-3-sphere topology (the "first closure" axiom) and test whether
the GEOMETRY relaxes to round S3. This module builds the substrate and the discrete
geometry it carries:

  generic_s3(N, seed) : a GENERIC combinatorial 3-sphere with N vertices, built as
                        the boundary of the 4D convex hull of N points sampled on S^3
                        (then coordinates are DISCARDED — only the abstract simplicial
                        complex + free edge lengths survive). Generic => no 2I symmetry
                        => the anti-circularity test on the main arm is meaningful.

  Regge calculus (3D): hinges are EDGES; the deficit around edge e is
                        delta_e = 2*pi - sum_{tets t ⊃ e} theta(t,e), theta = dihedral
                        angle. Regge action S = sum_e ell_e * delta_e, and Schlafli's
                        identity gives dS/d ell_e = delta_e (the gradient that drives
                        normalized discrete Ricci flow in Stage 1).

Dihedral angles are intrinsic (functions of the 6 edge lengths of a tet), so each tet
is embedded independently in R^3 from its edge lengths; an invalid (non-realizable)
tet -> NaN, which flags a non-embeddable metric.

Standalone:  python geometry.py        # build + verify a generic S^3
"""

import numpy as np
from scipy.spatial import ConvexHull
from itertools import combinations

# canonical within-tet edge order (pairs of local vertex indices 0..3)
TET_EDGES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
# for each within-tet edge, the two "other" local vertices (for dihedral faces)
TET_EDGE_OPP = {(0, 1): (2, 3), (0, 2): (1, 3), (0, 3): (1, 2),
                (1, 2): (0, 3), (1, 3): (0, 2), (2, 3): (0, 1)}


# ---------------------------------------------------------------------------
# Substrate: generic combinatorial 3-sphere
# ---------------------------------------------------------------------------
class Complex:
    """Abstract simplicial 3-sphere: vertices, tets (4-tuples), triangles, edges."""
    def __init__(self, n_vertices, tets):
        self.N = n_vertices
        self.tets = [tuple(sorted(t)) for t in tets]
        tri, edg = set(), set()
        for t in self.tets:
            for f in combinations(t, 3):
                tri.add(f)
            for e in combinations(t, 2):
                edg.add(e)
        self.triangles = sorted(tri)
        self.edges = sorted(edg)
        self.edge_index = {e: i for i, e in enumerate(self.edges)}
        # incidence: edge -> list of (tet_index, local edge pair)
        self._edge_tets = {e: [] for e in self.edges}
        for ti, t in enumerate(self.tets):
            loc = {v: k for k, v in enumerate(t)}
            for (a, b) in combinations(t, 2):
                self._edge_tets[(a, b)].append((ti, (loc[a], loc[b])))
        # triangle -> tets (for boundary/closedness check)
        self._tri_tets = {f: [] for f in self.triangles}
        for ti, t in enumerate(self.tets):
            for f in combinations(t, 3):
                self._tri_tets[f].append(ti)

    @property
    def euler(self):
        return self.N - len(self.edges) + len(self.triangles) - len(self.tets)


def _hull_complex(pts, N):
    hull = ConvexHull(pts)
    tets = [tuple(s) for s in hull.simplices]
    cx = Complex(N, tets)
    cx.hull_points = pts  # kept ONLY for the round-reference metric / verification
    return cx


def generic_s3(N, seed=20260618):
    """Generic combinatorial 3-sphere = boundary of the 4D convex hull of N pts on S^3."""
    rng = np.random.default_rng(seed)
    pts = rng.standard_normal((N, 4))
    pts /= np.linalg.norm(pts, axis=1, keepdims=True)
    return _hull_complex(pts, N)


def defective_s3(N, seed=20260618, concentration=6.0):
    """
    DEFECTIVE-mesh null (brief v1.1 control #2): points heavily clustered toward one
    pole, so the triangulation has dense and starved regions -> sliver tets that cannot
    reach round S^3 at fixed connectivity. A valid combinatorial S^3 (convex hull always
    is), but a known-can't-reach-round negative control: gives "the variance got small"
    its discriminative teeth.
    """
    rng = np.random.default_rng(seed)
    pts = rng.standard_normal((N, 4))
    pts[:, 0] += concentration            # bias toward the (1,0,0,0) pole
    pts /= np.linalg.norm(pts, axis=1, keepdims=True)
    return _hull_complex(pts, N)


# ---------------------------------------------------------------------------
# Verification: is it a closed combinatorial 3-sphere, and is it 2I-free?
# ---------------------------------------------------------------------------
def verify_sphere(cx, check_links=True, max_link_checks=200):
    """Return a dict of checks; a closed 3-sphere has chi=0, no boundary, S^2 links."""
    out = {}
    out["euler_char"] = cx.euler
    out["euler_ok"] = (cx.euler == 0)
    # closed: every triangle in exactly 2 tets
    counts = [len(v) for v in cx._tri_tets.values()]
    out["max_tri_tets"] = max(counts)
    out["min_tri_tets"] = min(counts)
    out["closed_no_boundary"] = (max(counts) == 2 and min(counts) == 2)
    # vertex links are 2-spheres (chi = 2)
    if check_links:
        bad = 0
        verts = list(range(cx.N))
        for v in verts[:max_link_checks]:
            link_tris = set()
            for t in cx.tets:
                if v in t:
                    others = tuple(sorted(x for x in t if x != v))
                    link_tris.add(others)
            # link is a 2-complex of triangles 'others'; compute its chi
            lv, le = set(), set()
            for f in link_tris:
                for x in f:
                    lv.add(x)
                for e in combinations(f, 2):
                    le.add(e)
            chi = len(lv) - len(le) + len(link_tris)
            if chi != 2:
                bad += 1
        out["links_checked"] = min(cx.N, max_link_checks)
        out["links_bad"] = bad
        out["links_ok"] = (bad == 0)
    return out


def symmetry_probe(cx):
    """
    Cheap 2I-freeness probe: the 600-cell (2I-symmetric) is vertex-transitive with a
    SINGLE vertex-degree and a single edge-degree (5 tets/edge). A generic complex has
    a spread of degrees. We report the spreads; large spread => not 2I-symmetric.
    """
    vert_deg = np.zeros(cx.N, dtype=int)
    for t in cx.tets:
        for v in t:
            vert_deg[v] += 1
    edge_tet_count = np.array([len(cx._edge_tets[e]) for e in cx.edges])
    return {
        "vertex_deg_unique": int(len(np.unique(vert_deg))),
        "vertex_deg_min": int(vert_deg.min()), "vertex_deg_max": int(vert_deg.max()),
        "edge_tets_unique": int(len(np.unique(edge_tet_count))),
        "edge_tets_min": int(edge_tet_count.min()),
        "edge_tets_max": int(edge_tet_count.max()),
        "mean_tets_per_edge": float(edge_tet_count.mean()),
        # 2I/600-cell fingerprint would be: all degrees equal AND edge_tets all == 5
        "looks_2I_regular": bool(len(np.unique(vert_deg)) == 1 and
                                 len(np.unique(edge_tet_count)) == 1 and
                                 edge_tet_count.max() == 5),
    }


# ---------------------------------------------------------------------------
# Regge geometry from edge lengths
# ---------------------------------------------------------------------------
def tet_embed(L):
    """
    Embed a tetrahedron in R^3 from its 6 edge lengths L = [l01,l02,l03,l12,l13,l23].
    Returns 4x3 coords, or None if non-realizable (negative square in the construction).
    """
    l01, l02, l03, l12, l13, l23 = L
    p0 = np.array([0.0, 0.0, 0.0])
    p1 = np.array([l01, 0.0, 0.0])
    if l01 <= 0:
        return None
    x2 = (l01**2 + l02**2 - l12**2) / (2 * l01)
    y2sq = l02**2 - x2**2
    if y2sq <= 1e-12:
        return None
    y2 = np.sqrt(y2sq)
    p2 = np.array([x2, y2, 0.0])
    x3 = (l01**2 + l03**2 - l13**2) / (2 * l01)
    y3 = (l02**2 + l03**2 - l23**2 - 2 * x2 * x3) / (2 * y2)
    z3sq = l03**2 - x3**2 - y3**2
    if z3sq <= 1e-12:
        return None
    z3 = np.sqrt(z3sq)
    p3 = np.array([x3, y3, z3])
    return np.array([p0, p1, p2, p3])


def tet_volume(L):
    """Signed-free volume of a tet from its embedding (0 if non-realizable)."""
    P = tet_embed(L)
    if P is None:
        return np.nan
    return abs(np.linalg.det(P[1:] - P[0])) / 6.0


def _dihedral(P, a, b, c, d):
    """Dihedral angle along edge (a,b) between faces (a,b,c) and (a,b,d)."""
    u = P[b] - P[a]; u /= np.linalg.norm(u)
    v = P[c] - P[a]; w = P[d] - P[a]
    vp = v - (v @ u) * u
    wp = w - (w @ u) * u
    nv, nw = np.linalg.norm(vp), np.linalg.norm(wp)
    if nv < 1e-12 or nw < 1e-12:
        return np.nan
    cos = np.clip((vp @ wp) / (nv * nw), -1.0, 1.0)
    return np.arccos(cos)


def tet_dihedrals(L):
    """6 dihedral angles in TET_EDGES order; NaN if non-realizable."""
    P = tet_embed(L)
    if P is None:
        return np.full(6, np.nan)
    ang = np.empty(6)
    for i, (a, b) in enumerate(TET_EDGES):
        c, d = TET_EDGE_OPP[(a, b)]
        ang[i] = _dihedral(P, a, b, c, d)
    return ang


def edge_deficits(cx, ell):
    """
    delta_e = 2*pi - sum of dihedral angles of incident tets at e.
    `ell` is a vector over cx.edges (in cx.edge_index order).
    Returns (deficits[over edges], n_invalid_tets).
    """
    two_pi = 2 * np.pi
    deficits = np.full(len(cx.edges), two_pi)
    n_invalid = 0
    # precompute each tet's 6 dihedrals once
    for ti, t in enumerate(cx.tets):
        loc = {v: k for k, v in enumerate(t)}
        L = np.array([ell[cx.edge_index[tuple(sorted((t[a], t[b])))]]
                      for (a, b) in TET_EDGES])
        ang = tet_dihedrals(L)
        if np.any(np.isnan(ang)):
            n_invalid += 1
            continue
        for i, (a, b) in enumerate(TET_EDGES):
            e = tuple(sorted((t[a], t[b])))
            deficits[cx.edge_index[e]] -= ang[i]
    return deficits, n_invalid


def total_volume(cx, ell):
    vol = 0.0
    for t in cx.tets:
        L = np.array([ell[cx.edge_index[tuple(sorted((t[a], t[b])))]]
                      for (a, b) in TET_EDGES])
        v = tet_volume(L)
        if not np.isnan(v):
            vol += v
    return vol


def regge_action(cx, ell):
    deficits, _ = edge_deficits(cx, ell)
    return float(np.sum(ell * deficits))


if __name__ == "__main__":
    cx = generic_s3(200, seed=20260618)
    print(f"generic S^3: N={cx.N}  edges={len(cx.edges)}  "
          f"triangles={len(cx.triangles)}  tets={len(cx.tets)}")
    chk = verify_sphere(cx)
    print("verify:", {k: chk[k] for k in
          ("euler_char", "euler_ok", "closed_no_boundary", "links_ok", "links_bad")})
    sym = symmetry_probe(cx)
    print("symmetry probe:", sym)
    ell = np.ones(len(cx.edges))
    deficits, ninv = edge_deficits(cx, ell)
    print(f"frustrated metric (all ell=1): invalid_tets={ninv}  "
          f"deficit mean={deficits.mean():.4f} std={deficits.std():.4f} "
          f"[min={deficits.min():.4f} max={deficits.max():.4f}]")
    print(f"total volume={total_volume(cx, ell):.4f}  Regge action={regge_action(cx, ell):.4f}")
    # regular-tet sanity: equilateral dihedral = arccos(1/3) = 1.23096 rad
    print(f"regular-tet dihedral check: {tet_dihedrals(np.ones(6))[0]:.5f} "
          f"(expect {np.arccos(1/3):.5f})")
