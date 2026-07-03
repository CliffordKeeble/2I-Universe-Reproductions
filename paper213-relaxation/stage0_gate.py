"""
stage0_gate.py  —  rigorous symmetry gate (brief v1.1 §1, Mr A #4).

Before Stage 1 we must show the generic complex is NOT secretly symmetric — otherwise
the equilateral ell=1 init smuggles in an isometry group (the metric's isometry group =
the combinatorial automorphism group when all edges are equal). The degree-spread screen
from the earlier run is NOT the gate; the automorphism computation is.

Method: Weisfeiler-Leman (1-WL) colour refinement on the 1-skeleton, seeded with each
vertex's simplicial signature (graph degree, #tets, sorted incident edge-tet counts).
Refine to a stable colouring. Any automorphism must preserve WL colours, so:

    all colour classes singletons  =>  Aut = {id} (RIGOROUS upper bound |Aut| = 1)
                                    =>  no order-120 / 2I subgroup => 2I-free. PASS.

If some class has size > 1, automorphism orbits are unions of those classes; we report
the largest class as an upper bound on any orbit and FLAG (the gate is not cleanly passed
and the equilateral init is not licensed).

Standalone:  python stage0_gate.py [N]
"""

import sys
import numpy as np
from collections import defaultdict
from geometry import generic_s3


def wl_refine(cx, max_iter=50):
    """Return stable WL colour per vertex (as integers 0..C-1)."""
    N = cx.N
    # adjacency from edges
    adj = defaultdict(list)
    for (a, b) in cx.edges:
        adj[a].append(b); adj[b].append(a)
    # simplicial seed signature per vertex
    tets_per_vertex = np.zeros(N, dtype=int)
    for t in cx.tets:
        for v in t:
            tets_per_vertex[v] += 1
    edge_tet_count = {e: len(cx._edge_tets[e]) for e in cx.edges}
    seed = []
    for v in range(N):
        inc = sorted(edge_tet_count[tuple(sorted((v, u)))] for u in adj[v])
        seed.append((len(adj[v]), int(tets_per_vertex[v]), tuple(inc)))
    # initial colours by distinct seed
    colour = _relabel(seed)
    for _ in range(max_iter):
        new = [(colour[v], tuple(sorted(colour[u] for u in adj[v]))) for v in range(N)]
        new_colour = _relabel(new)
        if max(new_colour) == max(colour) and new_colour == colour:
            break
        if max(new_colour) == max(colour):
            colour = new_colour
            break
        colour = new_colour
    return colour


def _relabel(signatures):
    order = {}
    out = []
    for s in signatures:
        if s not in order:
            order[s] = len(order)
        out.append(order[s])
    return out


def gate(cx, verbose=True):
    colour = wl_refine(cx)
    classes = defaultdict(int)
    for c in colour:
        classes[c] += 1
    n_classes = len(classes)
    max_class = max(classes.values())
    trivial_aut = (max_class == 1)
    res = dict(N=cx.N, n_colour_classes=n_classes, max_class_size=max_class,
               aut_trivial=bool(trivial_aut),
               aut_upper_bound="1 (identity only)" if trivial_aut
               else f"orbits <= largest class = {max_class}",
               twoI_free=bool(trivial_aut),  # trivial Aut => no order-120 subgroup
               PASS=bool(trivial_aut))
    if verbose:
        print(f"Stage 0 gate (N={cx.N}): colour classes={n_classes}/{cx.N}, "
              f"max class={max_class}")
        print(f"  Aut upper bound: {res['aut_upper_bound']}")
        print(f"  2I-free: {res['twoI_free']}   GATE PASS: {res['PASS']}")
    return res


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    cx = generic_s3(N, seed=20260618)
    gate(cx)
