"""
laplacian.py -- Paper 213 Arm 1: the two read-out operators (brief sec.4).

  L_s : scalar (combinatorial) Laplacian  D - A  of the underlying simple graph.
        Carries metric/topology -> spectral dimension d_s. IDENTICAL for Arm 1
        and Arm 2 under the same growth rule (d_s is topology, not coupling).

  L_c : connection (bundle) Laplacian over the R^2 fibres. Real-symmetric block
        matrix: diagonal block deg(i)*I2, off-diagonal block (i,j) = -U[i,j],
        with U[j,i] = U[i,j]^T. Carries holonomy content.

Note on L_c well-definedness (the seed-orthogonality wrinkle, see seed.py /
RESULTS.md): R is NOT orthogonal, but L_c is still REAL-SYMMETRIC because we set
the (j,i) block to U[i,j]^T = (-U[i,j])^T. Symmetry -> real spectrum, regardless
of orthogonality. We build exactly the operator the brief defines.
"""

import numpy as np
import scipy.sparse as sp


def scalar_laplacian(G):
    """L_s = D - A as a CSR matrix, nodes indexed 0..N-1 by sorted order."""
    nodes = sorted(G.nodes())
    idx = {n: i for i, n in enumerate(nodes)}
    N = len(nodes)
    rows, cols, vals = [], [], []
    for n in nodes:
        i = idx[n]
        d = G.degree(n)
        rows.append(i); cols.append(i); vals.append(float(d))
        for m in G.neighbors(n):
            j = idx[m]
            rows.append(i); cols.append(j); vals.append(-1.0)
    L = sp.csr_matrix((vals, (rows, cols)), shape=(N, N))
    return L, idx


def connection_laplacian(G, elem_mat):
    """L_c over R^2 fibres: 2N x 2N real-symmetric CSR matrix.

    Diagonal block i  : deg(i) * I2
    Off-diagonal (s,d): -U  where U = elem_mat[elem] for the stored orientation
                        src=s -> dst=d ; and block (d,s) = -U^T.
    """
    nodes = sorted(G.nodes())
    idx = {n: i for i, n in enumerate(nodes)}
    N = len(nodes)
    rows, cols, vals = [], [], []

    def put_block(bi, bj, M):
        for a in range(2):
            for b in range(2):
                v = M[a, b]
                if v != 0.0:
                    rows.append(2 * bi + a)
                    cols.append(2 * bj + b)
                    vals.append(float(v))

    # diagonal deg*I2
    for n in nodes:
        i = idx[n]
        put_block(i, i, G.degree(n) * np.eye(2))

    # off-diagonal -U / -U^T
    for u, w, data in G.edges(data=True):
        s = data["src"]
        d = w if s == u else u           # destination is the other endpoint
        U = elem_mat[data["elem"]]
        si, di = idx[s], idx[d]
        put_block(si, di, -U)
        put_block(di, si, -U.T)

    L = sp.csr_matrix((vals, (rows, cols)), shape=(2 * N, 2 * N))
    return L, idx


def check_symmetry(L, tol=1e-9):
    """Max |L - L^T| -- L_c must be symmetric for a real spectrum."""
    return abs((L - L.T)).max() if L.nnz else 0.0


if __name__ == "__main__":
    from growth import grow
    G, meta, elem_mat = grow(500, d_max=4, C="flat-triangles")
    Ls, _ = scalar_laplacian(G)
    Lc, _ = connection_laplacian(G, elem_mat)
    print(f"N={meta['N']}  L_s {Ls.shape} nnz={Ls.nnz}  "
          f"L_c {Lc.shape} nnz={Lc.nnz}")
    print(f"  L_s symmetric? max|L-L^T| = {check_symmetry(Ls):.2e}")
    print(f"  L_c symmetric? max|L-L^T| = {check_symmetry(Lc):.2e}  "
          f"(must be ~0 despite R non-orthogonal)")
    # row sums of L_s should be 0 (Laplacian)
    rs = np.abs(np.asarray(Ls.sum(axis=1)).ravel()).max()
    print(f"  L_s max|row sum| = {rs:.2e}  (0 => valid graph Laplacian)")
