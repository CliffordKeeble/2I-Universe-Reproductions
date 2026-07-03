"""
growth.py -- Paper 213 Arm 1: self-reference rewrite growth rule GR-1.

Builds a coupling graph G with NO coordinates anywhere. Nodes are reference
events; edges are couplings carrying a holonomy word in {I, R}. Time = causal
order of rewrites (node 'birth' index).

GR-1 (brief sec.3) is sketched, not fully specified. The one design decision the
brief leaves open is HOW CYCLES FORM -- and cycles are what make the closure
predicate C bite and what set the spectral dimension d_s. This module makes that
explicit and LOGS every (rho, C, d_max) variant, exactly as the brief demands.

  *** This concrete cycle rule is Mr Code's instantiation of the brief's sketch.
      It moves d_s (which the brief pre-registers as OPEN) but NOT the signature
      verdict (which rests on <I,R> = Z/2, independent of growth topology). ***

Cycle mechanism: when frontier node n spawns child c, c is connected to n AND
(optionally) to one existing neighbour 'sib' of n, closing a triangle
(n, c, sib). The closure predicate C decides whether that closing edge is kept:

  C = 'none'           : always close (max cycle density)
  C = 'flat-triangles' : close iff triangle holonomy is flat (product = I);
                         first try the other element on edge n->c to repair flatness
  C = 'flat-shortest'  : sib chosen 2 hops from n; close iff the *shortest* cycle
                         the closing edge creates is flat (heavier; debug/medium N)

Edge holonomy element is set by rho (the independent variable, brief sec.3):
  GR-1 default:  rho(n) = R if deg(n) is odd else I.
For Arm 1 every rho output is in {I, R}, so the holonomy stays abelian (Z/2).
"""

from collections import deque
import json
import numpy as np
import networkx as nx

from seed import ELEMENTS, is_flat, reflection_from_offdiag, I2

# ---------------------------------------------------------------------------
# rho variants (per-edge element rule)
# ---------------------------------------------------------------------------

def rho_GR1(deg):
    """GR-1 default: R on odd degree, I on even."""
    return "R" if (deg % 2 == 1) else "I"

def rho_allR(deg):
    """Diagnostic: every edge carries R."""
    return "R"


# ---------------------------------------------------------------------------
# Growth
# ---------------------------------------------------------------------------

def grow(N_target, d_max=4, k=2, C="flat-triangles", rho=rho_GR1,
         seed_a=5.0, rng_seed=20260618, scramble=False):
    """
    Grow a GR-1 graph.

      N_target : target node count (growth stops at >= this, or frontier empties)
      d_max    : branching cap (free choice; SWEPT, never tuned -- brief sec.3/5.4)
      k        : max children proposed per frontier node per visit
      C        : closure predicate variant in {'none','flat-triangles','flat-shortest'}
      rho      : per-edge element rule (the independent variable)
      seed_a   : off-diagonal of the order-2 involution; 5 = golden seed R,
                 2/3 = non-golden controls (sec.5.1)
      scramble : if True, assign edge elements as random {I,R} (control, sec.5.2)

    Returns (G, meta). G is an undirected networkx Graph; each edge carries
      'elem' (label in {'I','R'}) and 'src' (the endpoint whose rho generated it,
      needed to orient the non-symmetric R for the connection Laplacian).
    Node 0 carries the seed self-reference; v[0] = (1,0). (The seed self-loop is
    recorded in meta, NOT added as a literal graph self-loop, so that L_s = D - A
    stays the standard combinatorial Laplacian. Documented in RESULTS.md.)
    """
    rng = np.random.default_rng(rng_seed)
    # The element matrices for this run (golden seed unless a control overrides a).
    elem_mat = {"I": I2, "R": reflection_from_offdiag(seed_a)}

    G = nx.Graph()
    G.add_node(0, v=np.array([1.0, 0.0]), birth=0)
    frontier = deque([0])
    next_id = 1
    tick = 0
    saturated = set()

    def pick_elem(parent_deg):
        if scramble:
            return "R" if rng.random() < 0.5 else "I"
        return rho(parent_deg)

    def add_edge(a, b, label, src):
        G.add_edge(a, b, elem=label, src=src)

    while G.number_of_nodes() < N_target and frontier:
        n = frontier.popleft()
        if n in saturated or G.degree(n) >= d_max:
            continue
        tick += 1
        for _ in range(k):
            if G.number_of_nodes() >= N_target or G.degree(n) >= d_max:
                break
            c = next_id
            g_nc = pick_elem(G.degree(n))

            # candidate sibling for triangle closure
            sib = None
            if C != "none-leaf":
                neigh = [m for m in G.neighbors(n) if G.degree(m) < d_max]
                if C == "flat-shortest":
                    # 2-hop: a neighbour of a neighbour, not n itself
                    pool = []
                    for m in G.neighbors(n):
                        pool += [w for w in G.neighbors(m)
                                 if w != n and G.degree(w) < d_max]
                    pool = list(dict.fromkeys(pool))
                    sib = rng.choice(pool) if pool else None
                else:
                    sib = rng.choice(neigh) if neigh else None
                if sib is not None:
                    sib = int(sib)

            # create the child + tree edge n->c
            G.add_node(c, v=np.array([1.0, 0.0]), birth=tick)
            add_edge(n, c, g_nc, src=n)
            next_id += 1
            frontier.append(c)

            # decide on the closing edge c--sib
            if sib is None or G.degree(c) >= d_max or G.degree(sib) >= d_max:
                continue
            g_sc = pick_elem(G.degree(sib))

            if C == "none":
                add_edge(sib, c, g_sc, src=sib)
            elif C == "flat-triangles":
                g_sn = G[sib][n]["elem"]
                if is_flat_mat([g_nc, g_sc, g_sn], elem_mat):
                    add_edge(sib, c, g_sc, src=sib)
                else:
                    # "try another g" on the n->c edge to repair flatness
                    other = "I" if g_nc == "R" else "R"
                    if is_flat_mat([other, g_sc, g_sn], elem_mat):
                        G[n][c]["elem"] = other
                        add_edge(sib, c, g_sc, src=sib)
                    # else: leave c as a leaf (no cycle); n may later saturate
            elif C == "flat-shortest":
                try:
                    path = nx.shortest_path(G, sib, n)
                except nx.NetworkXNoPath:
                    path = None
                if path is not None:
                    cyc_labels = [g_nc, g_sc]
                    for u, w in zip(path[:-1], path[1:]):
                        cyc_labels.append(G[u][w]["elem"])
                    if is_flat_mat(cyc_labels, elem_mat):
                        add_edge(sib, c, g_sc, src=sib)

        if G.degree(n) < d_max:
            frontier.append(n)   # revisit later
        else:
            saturated.add(n)

    meta = {
        "N": G.number_of_nodes(),
        "E": G.number_of_edges(),
        "d_max": d_max,
        "k": k,
        "C": C,
        "rho": rho.__name__,
        "seed_a": seed_a,
        "scramble": scramble,
        "rng_seed": rng_seed,
        "ticks": tick,
        "seed_self_loop_elem": "R",
        "mean_degree": 2.0 * G.number_of_edges() / max(1, G.number_of_nodes()),
        "n_triangles": int(sum(nx.triangles(G).values()) // 3),
    }
    return G, meta, elem_mat


def is_flat_mat(labels, elem_mat, tol=1e-9):
    """Flatness using this run's element matrices (handles control seed_a)."""
    M = I2.copy()
    for lab in labels:
        M = M @ elem_mat[lab]
    return bool(np.allclose(M, I2, atol=tol))


if __name__ == "__main__":
    G, meta, _ = grow(1000, d_max=4, k=2, C="flat-triangles")
    print("GR-1 debug grow (N=1000, d_max=4, C=flat-triangles):")
    print(json.dumps(meta, indent=2, default=str))
    print("connected:", nx.is_connected(G))
