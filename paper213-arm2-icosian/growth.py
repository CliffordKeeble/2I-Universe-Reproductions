"""
growth.py  —  GR-1 self-reference growth engine for Paper 213 (Arms 1 & 2).

Coupling-agnostic by design. The TOPOLOGY (which nodes/edges exist) is generated
from the RNG alone in the 'no_constraint' closure variant, so it is identical for
any coupling under the same seed — this is what lets compare_arms.py assert that
L_s (hence the spectral dimension d_s) matches across arms, isolating the COUPLING
as the only difference (Arm 1 §4).

A `Coupling` supplies:
  - fibre_dim  : 2 (Arm 1, golden spinor) or 4 (Arm 2, icosian fibre)
  - identity   : I_f
  - rho(G, n)  : per-edge element rule -> a holonomy matrix (the independent variable)
  - name / kind

ARM 2 NOTE: when the coupling is the icosian coupling, 2I IS IMPORTED — the holonomy
alphabet {s,t} is hand-coded 2I. Nothing emerges; we test reconstruction only.

The flatness-gated closure variants (flat_triangles / flat_shortest) DO couple
topology to holonomy and so are used only inside the per-arm rule sweep, never for
the cross-arm comparison.

Standalone smoke test:  python growth.py
"""

import numpy as np
import networkx as nx

TOL = 1e-7  # flatness tolerance for closure predicates (finite group, well separated)


# ---------------------------------------------------------------------------
# Couplings
# ---------------------------------------------------------------------------
class Coupling:
    """Base coupling: identity holonomy on a 1-D fibre (degenerate; for testing)."""
    name = "trivial"
    kind = "trivial"
    fibre_dim = 1

    def __init__(self):
        self.identity = np.eye(self.fibre_dim)

    def rho(self, G, n):
        return np.eye(self.fibre_dim)


class SeedCoupling(Coupling):
    """
    ARM 1 coupling. Import nothing about 2I.
    Seed Gamma = [[0,1],[sqrt5,0]]; connection-admissible normalisation
    R = Gamma / 5^(1/4) is the order-2 reflection (det -1, R^2 = I).
    Edge holonomies are words in {I, R} -> abelian (Z/2). GR-1 default:
    rho(n) = R if deg(n) is odd else I.
    """
    name = "arm1_seed"
    kind = "seed"
    fibre_dim = 2

    def __init__(self, root=5.0):
        s = root ** 0.25
        self.R = np.array([[0.0, 1.0], [root, 0.0]]) / s   # det -1, R^2 = I
        self.I = np.eye(2)
        self.identity = self.I

    def rho(self, G, n):
        return self.R if (G.degree(n) % 2 == 1) else self.I


class IcosianCoupling(Coupling):
    """
    ARM 2 coupling.  *** 2I IS IMPORTED — this coupling is the import. ***
    Holonomy alphabet = a verified generating pair {s,t} of 2I, as 4x4 SO(4)
    left-multiplication matrices. GR-1 element rule: rho(n) = s if deg(n) odd
    else t  (both nontrivial, non-abelian, 5-torsion present -> can carry 2I).
    A deterministic context-hash into the full 120 is available via rho_hash.
    """
    name = "arm2_icosian"
    kind = "icosian"
    fibre_dim = 4

    def __init__(self):
        from icosians import build_icosians, find_generator_pair, left_matrix
        G2I = build_icosians()
        i, j, oi, oj = find_generator_pair(G2I)
        self.s = left_matrix(G2I[i])
        self.t = left_matrix(G2I[j])
        self.s_order, self.t_order = oi, oj
        self.I = np.eye(4)
        self.identity = self.I
        self._mats = np.array([left_matrix(q) for q in G2I])  # for context-hash

    def rho(self, G, n):
        return self.s if (G.degree(n) % 2 == 1) else self.t

    def rho_hash(self, G, n):
        """Deterministic context-hash of node n into the full 120-element 2I."""
        idx = (hash(("icosian", int(n), int(G.degree(n)))) % 120 + 120) % 120
        return self._mats[idx]


# ---------------------------------------------------------------------------
# Holonomy bookkeeping
# ---------------------------------------------------------------------------
def _cycle_holonomy(hol, cycle_nodes):
    """Product of oriented edge holonomies around a closed node cycle."""
    f = hol[(cycle_nodes[0], cycle_nodes[0])].shape[0] if False else None
    M = None
    for a, b in zip(cycle_nodes, cycle_nodes[1:] + cycle_nodes[:1]):
        U = hol[(a, b)]
        M = U if M is None else U @ M
    return M


def _is_flat(M):
    return np.linalg.norm(M - np.eye(M.shape[0])) < TOL


# ---------------------------------------------------------------------------
# GR-1 growth
# ---------------------------------------------------------------------------
def grow(N, coupling, d_max=4, k_children=2, p_close=0.5,
         closure="no_constraint", seed=20260618):
    """
    Grow a GR-1 graph to N nodes.

    Returns (G, hol) where:
      G   : networkx.Graph (no coordinates anywhere; nodes are reference events)
      hol : dict {(i,j): U_ij} oriented edge holonomies, with hol[(j,i)] = U_ij^T

    closure variants:
      'no_constraint'  : topology from RNG only -> coupling-independent (use for compare).
      'flat_triangles' : a closure edge is kept only if the triangle holonomy is flat.
      'flat_shortest'  : closure edge kept only if the shortest induced cycle is flat.

    Determinism: fixed `seed`; no unseeded randomness.
    """
    rng = np.random.default_rng(seed)
    G = nx.Graph()
    hol = {}

    def set_edge(a, b, U):
        G.add_edge(a, b)
        hol[(a, b)] = U
        hol[(b, a)] = U.T

    # Seed: node 0 with a self-loop holonomy (recorded in hol, not as a graph self-edge
    # to keep L_s simple; the self-loop's holonomy seeds rho via degree parity).
    G.add_node(0)
    self_U = coupling.rho(G, 0)
    hol[(0, 0)] = self_U  # bookkeeping only

    frontier = [0]
    next_id = 1

    while G.number_of_nodes() < N and frontier:
        n = frontier.pop(0)
        if G.degree(n) >= d_max:
            continue
        for _ in range(k_children):
            if G.number_of_nodes() >= N:
                break
            c = next_id
            next_id += 1
            U = coupling.rho(G, n)
            set_edge(n, c, U)
            frontier.append(c)

            # closure: optionally add a cross-edge making a cycle
            if rng.random() < p_close:
                cands = [m for m in frontier
                         if m != c and m != n and G.degree(m) < d_max
                         and not G.has_edge(c, m)]
                if cands:
                    m = int(rng.choice(cands))
                    Uc = coupling.rho(G, c)
                    keep = True
                    if closure == "flat_triangles":
                        # triangle c-m-n-c if n-m edge exists, else the cycle c-m-...-n-c
                        if G.has_edge(m, n):
                            tri = [c, m, n]
                            # tentatively place edge to test flatness
                            hol[(c, m)] = Uc
                            hol[(m, c)] = Uc.T
                            keep = _is_flat(_cycle_holonomy(hol, tri))
                            if not keep:
                                del hol[(c, m)]
                                del hol[(m, c)]
                        else:
                            keep = True
                    elif closure == "flat_shortest":
                        try:
                            path = nx.shortest_path(G, c, m)
                            cyc = path + [c]
                            hol[(c, m)] = Uc
                            hol[(m, c)] = Uc.T
                            keep = _is_flat(_cycle_holonomy(hol, cyc[:-1]))
                            if not keep:
                                del hol[(c, m)]
                                del hol[(m, c)]
                        except nx.NetworkXNoPath:
                            keep = True
                    if keep and not G.has_edge(c, m):
                        set_edge(c, m, Uc)

            if G.degree(n) >= d_max:
                break
        # n may still have capacity; if so, keep it available
        if G.degree(n) < d_max and n not in frontier:
            frontier.append(n)

    return G, hol


def topology_fingerprint(G):
    """A cheap invariant to confirm two graphs are the same topology."""
    return (G.number_of_nodes(), G.number_of_edges(),
            tuple(sorted(d for _, d in G.degree())))


if __name__ == "__main__":
    for cpl in (SeedCoupling(), IcosianCoupling()):
        G, hol = grow(2000, cpl, seed=20260618)
        print(f"{cpl.name:14s} fibre={cpl.fibre_dim}  "
              f"N={G.number_of_nodes()} E={G.number_of_edges()} "
              f"<deg>={2*G.number_of_edges()/G.number_of_nodes():.2f}  "
              f"connected={nx.is_connected(G)}")
    # identical-topology check under no_constraint (same seed)
    Ga, _ = grow(2000, SeedCoupling(), seed=20260618)
    Gb, _ = grow(2000, IcosianCoupling(), seed=20260618)
    print("topology identical across arms (no_constraint):",
          topology_fingerprint(Ga) == topology_fingerprint(Gb))
