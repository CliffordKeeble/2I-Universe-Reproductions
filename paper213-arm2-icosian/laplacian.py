"""
laplacian.py  —  scalar and connection Laplacians for Paper 213.

  L_s  (scalar, combinatorial) = D - A of the underlying graph.
        Carries metric/topology -> spectral dimension. IDENTICAL across arms
        under the same growth rule (d_s is a property of topology, not coupling).

  L_c  (connection / bundle) over the f-dim fibres:
        Hermitian block matrix, diagonal block deg(i)*I_f,
        off-diagonal block (i,j) = -U[i,j].
        Carries the holonomy content (Z/2 for Arm 1, 2I for Arm 2).

Both returned as scipy.sparse CSR (real symmetric).
"""

import numpy as np
import scipy.sparse as sp
import networkx as nx


def scalar_laplacian(G):
    """L_s = D - A, sparse, nodes ordered by sorted(G.nodes())."""
    nodes = sorted(G.nodes())
    return nx.laplacian_matrix(G, nodelist=nodes).astype(float).tocsr()


def connection_laplacian(G, hol, fibre_dim):
    """
    Block connection Laplacian. Size (N*f) x (N*f).
    Diagonal block i  : deg(i) * I_f
    Off-diagonal (i,j): -U[i,j]   (U[j,i] = U[i,j]^T -> overall symmetric)
    """
    nodes = sorted(G.nodes())
    index = {n: k for k, n in enumerate(nodes)}
    f = fibre_dim
    N = len(nodes)

    rows, cols, data = [], [], []

    # diagonal blocks
    for n in nodes:
        i = index[n]
        deg = G.degree(n)
        for a in range(f):
            rows.append(i * f + a)
            cols.append(i * f + a)
            data.append(float(deg))

    # off-diagonal blocks
    for (u, v) in G.edges():
        i, j = index[u], index[v]
        U = hol[(u, v)]  # block at (i,j) is -U ; (j,i) is -U^T
        for a in range(f):
            for b in range(f):
                val = -U[a, b]
                if val != 0.0:
                    rows.append(i * f + a); cols.append(j * f + b); data.append(val)
                    rows.append(j * f + b); cols.append(i * f + a); data.append(val)

    M = sp.coo_matrix((data, (rows, cols)), shape=(N * f, N * f)).tocsr()
    # symmetrise defensively against float asymmetry
    M = 0.5 * (M + M.T)
    return M


if __name__ == "__main__":
    from growth import SeedCoupling, IcosianCoupling, grow
    for cpl in (SeedCoupling(), IcosianCoupling()):
        G, hol = grow(500, cpl, seed=20260618)
        Ls = scalar_laplacian(G)
        Lc = connection_laplacian(G, hol, cpl.fibre_dim)
        symm = abs(Lc - Lc.T).max()
        print(f"{cpl.name:14s}  L_s {Ls.shape}  L_c {Lc.shape}  "
              f"|Lc-Lc^T|max={symm:.2e}")
