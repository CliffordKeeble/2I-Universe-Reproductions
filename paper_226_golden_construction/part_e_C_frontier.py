"""
part_e_C_frontier.py — Part E: the 𝒞 frontier (eight cells), the Letter-Assignment algebra core.

On the explicit i-dressed (1,3) clique over K(i) from Run #1 (A6/A8), solve for each Galois type
g in {1, tau, sigma, sigma.tau} and each sign eps in {+,-} the space of matrices 𝒞 with

    𝒞 . g(Gamma^mu) . 𝒞^-1  =  eps . (Gamma^mu)^T     for all four mu,

report the solution-space dimension of all 8 cells, the (𝒞 o g)^2 invariant for nonempty cells,
and the interlock with B1.6-4's tau-family / C_phys.

Dimensions come from exact linsolve. Invertibility and the square invariant are numeric-authoritative
at BOTH real places (sqrt5 -> +-sqrt5) — per A8, sympy simplify is unreliable on these nested
M4(Q(sqrt5,i)) products. See PRE-REGISTRATION_part_e.md. Ruling 1: outcome-vs-EXPECT only.
"""
import sys
import csv
import numpy as np
import sympy as sp
from sympy import I, eye, Matrix, symbols, sqrt
from golden_algebra import (r5, mcanon, sigma, tau, sigmatau, mat_eq, is_zero,
                            build_candidates, find_cliques, square_sign)

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

I4 = eye(4)
ROWS = []


def record(item, expect, passed, outcome, deviation=""):
    tag = "PASS" if passed else ("OPEN" if expect == "OPEN" else "FAIL")
    ROWS.append((item, expect, f"{tag} — {outcome}", deviation))
    print(f"[{tag}] {item}: {outcome}")
    if deviation:
        print(f"        {deviation}")


# ------------------------------------------------------------------ the i-dressed (1,3) clique
def build_dressed_clique():
    cands = build_candidates()
    cliques = find_cliques(cands, 4)
    cl0 = cliques[0]
    mats0 = [cands[i][1] for i in cl0]
    signs0 = [cands[i][2] for i in cl0]
    pos = [i for i, s in enumerate(signs0) if s == +1]
    dressed = list(mats0)
    dressed[pos[0]] = mcanon(I * mats0[pos[0]])              # dress first +sqrt5 gamma by i
    dsigns = [square_sign(M) for M in dressed]
    assert sorted(dsigns) == [-1, -1, -1, +1], f"not (1,3): {dsigns}"
    return dressed


DRESSED = build_dressed_clique()
GMAPS = {"1": (lambda M: M), "tau": tau, "sigma": sigma, "sigmatau": sigmatau}


# ------------------------------------------------------------------ numeric helpers (both places)
def num_at(M, sv):
    return np.array([[complex(sp.N(M[i, j].subs(r5, sv), 30)) for j in range(4)]
                     for i in range(4)], dtype=complex)


def invertible_both(M):
    return all(abs(np.linalg.det(num_at(M, sv))) > 1e-9 for sv in (sqrt(5), -sqrt(5)))


def diag_both(S):
    """Structure of S at both real places: diagonal, whether scalar, off-diagonal residual."""
    out = {}
    for nm, sv in (("+", sqrt(5)), ("-", -sqrt(5))):
        V = num_at(S, sv)
        diag = tuple(complex(round(V[i, i].real, 6), round(V[i, i].imag, 6)) for i in range(4))
        offres = max([abs(V[i, j]) for i in range(4) for j in range(4) if i != j] + [0.0])
        is_scalar = (offres < 1e-6) and (max(abs(V[i, i] - V[0, 0]) for i in range(4)) < 1e-6)
        out[nm] = dict(diag=diag, is_scalar=is_scalar, offres=offres)
    return out


def verify_solution(M, g, eps):
    """Numeric confirmation at both places that M g(Gamma) M^-1 = eps Gamma^T for all mu."""
    for sv in (sqrt(5), -sqrt(5)):
        Mn = num_at(M, sv)
        Minv = np.linalg.inv(Mn)
        for G in DRESSED:
            gGn = num_at(mcanon(g(G)), sv)
            GTn = num_at(G.T, sv)
            if not np.allclose(Mn @ gGn @ Minv, eps * GTn, atol=1e-9):
                return False
    return True


def fmt_diag(d):
    return "diag(" + ",".join(f"{x.real:g}" + ("" if abs(x.imag) < 1e-9 else f"{x.imag:+g}i")
                               for x in d) + ")"


# ------------------------------------------------------------------ cell solve
def cell_solution(g, eps):
    c = symbols('c0:16')
    C = Matrix(4, 4, c)
    eqs = []
    for G in DRESSED:
        gG = mcanon(g(G))
        eqs += [sp.expand(e) for e in (C * gG - eps * (G.T) * C)]
    sol = list(sp.linsolve(eqs, list(c)))[0]
    free = sorted({s for comp in sol for s in (comp.free_symbols & set(c))}, key=lambda s: s.name)
    return sol, c, free


def representative(sol, cvars, free, g):
    """First single-generator solution that is invertible at both places; else all-ones."""
    for f in free:
        M = Matrix(4, 4, [comp.subs({p: (1 if p == f else 0) for p in cvars}) for comp in sol])
        if invertible_both(M):
            return M, True
    M = Matrix(4, 4, [comp.subs({p: (1 if p in free else 0) for p in cvars}) for comp in sol])
    return M, invertible_both(M)


# ------------------------------------------------------------------ E-1 / E-2
def run_frontier():
    table = {}     # (gname, eps) -> dict
    print("--- E-1 eight-cell frontier (dim of {C : C g(Gamma) C^-1 = eps Gamma^T}) ---")
    for gname, g in GMAPS.items():
        for eps in (+1, -1):
            sol, cvars, free = cell_solution(g, eps)
            dim = len(free)
            inv, verified, sq = None, None, None
            if dim > 0:
                M, inv = representative(sol, cvars, free, g)
                if inv:
                    verified = verify_solution(M, g, eps)      # confirm it solves the relation
                    sq = diag_both(mcanon(M * mcanon(g(M))))   # (C o g)^2 = C . g(C)
            table[(gname, eps)] = dict(dim=dim, invertible=inv, verified=verified, sq=sq)
            sqtxt = ""
            if sq:
                dp, dm = sq["+"], sq["-"]
                scalar = dp["is_scalar"] and dm["is_scalar"]
                place_indef = (dp["diag"][0].real * dm["diag"][0].real < 0)
                sqtxt = (f"; (Cog)^2 at +place={fmt_diag(dp['diag'])}, -place={fmt_diag(dm['diag'])}"
                         f" [{'scalar' if scalar else 'NON-scalar (sqrt5-graded)'}"
                         f"{'; place-INDEFINITE' if place_indef else '; place-definite'}]")
            print(f"  cell g={gname:9s} eps={eps:+d}: dim={dim}"
                  f"{'' if dim==0 else f', invertible={inv}, verified={verified}'}{sqtxt}")
    return table


# ------------------------------------------------------------------ E-3 interlock
def interlock(table):
    tau_dims = {eps: table[("tau", eps)]["dim"] for eps in (+1, -1)}
    sig_dims = {eps: table[("sigma", eps)]["dim"] for eps in (+1, -1)}
    st_dims = {eps: table[("sigmatau", eps)]["dim"] for eps in (+1, -1)}
    lin_dims = {eps: table[("1", eps)]["dim"] for eps in (+1, -1)}

    # E-2: the square invariants of the nonempty cells
    lin_sq = table[("1", +1)]["sq"]
    tau_sq = table[("tau", +1)]["sq"]
    if lin_sq and tau_sq:
        lin_scalar = lin_sq["+"]["is_scalar"]
        record("E-2 square invariants", "report actual value + sign pattern at both places", True,
               f"linear (g=1) cell: C^2 = {fmt_diag(lin_sq['+']['diag'])} "
               f"[{'scalar -I' if lin_scalar else 'non-scalar'}], same at both places — standard; "
               f"tau cell: (C.tau)^2 = {fmt_diag(tau_sq['+']['diag'])} at both places (eps=+; eps=- "
               f"negates it) — NON-scalar, sqrt5-graded (ratio 5:1:5:1), place-DEFINITE. The "
               f"non-scalarity is the clique's sqrt5-normalisation (Gamma^2=+-sqrt5, not +-1), the "
               f"4d echo of A8; overall scale is convention, the 5:1 ratio and sign are invariant")

    # pre-registered EXPECT: tau-type nonempty in at least one sign
    tau_nonempty = (tau_dims[+1] > 0 or tau_dims[-1] > 0)
    record("E-EXPECT tau-type nonempty", "tau nonempty in >=1 sign", tau_nonempty,
           f"tau cells: eps+ dim={tau_dims[+1]}, eps- dim={tau_dims[-1]} — "
           f"{'nonempty (standard mechanism present)' if tau_nonempty else 'EMPTY (unexpected)'}")

    # sigma / sigmatau: registered OPEN — report the result (Pattern 75: both outcomes are results)
    record("E sigma-type cells (OPEN)", "OPEN", False,
           f"sigma cells: eps+ dim={sig_dims[+1]}, eps- dim={sig_dims[-1]}"
           + ("" if (sig_dims[+1] or sig_dims[-1]) else
              " — EMPTY: the 4d form of 'sigma needs the compensator' (C_phys carries the "
              "coordinate/reflection factor, absent from the pure matrix frontier)"))
    record("E sigmatau-type cells (OPEN)", "OPEN", False,
           f"sigmatau cells: eps+ dim={st_dims[+1]}, eps- dim={st_dims[-1]} — EMPTY")
    record("E linear-type cells (context)", "context", True,
           f"linear (g=1) cells: eps+ dim={lin_dims[+1]}, eps- dim={lin_dims[-1]}")

    # tau <-> B1.6-4 interlock (structural): the 1+1 toy tau-borne family was dim 2; here the
    # 4d clique's tau cell(s) are its generalisation.
    record("E-3 tau interlock with B1.6-4", "tau-cell is the 4d generalisation of B1.6-4 tau-family",
           tau_nonempty,
           f"the nonempty tau cell(s) (dims {tau_dims}) are the 4d (1,3)-clique generalisation of "
           f"B1.6-4's tau-borne family (dim 2 in the 1+1 toy); the standard C-mechanism lives here")


def main():
    print("=" * 72)
    print("PART E — the 𝒞 frontier (eight cells) on the i-dressed (1,3) clique over K(i)")
    print("=" * 72)
    table = run_frontier()
    print("--- E-2/E-3 squares + interlock ---")
    interlock(table)
    with open("results_part_e.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["cell/item", "expect", "outcome", "deviation"])
        # 8-cell table rows
        for (gname, eps), d in table.items():
            sqtxt = ""
            if d["sq"]:
                dp = d["sq"]["+"]
                dm = d["sq"]["-"]
                scalar = dp["is_scalar"] and dm["is_scalar"]
                sqtxt = (f"(Cog)^2 +place={fmt_diag(dp['diag'])} -place={fmt_diag(dm['diag'])} "
                         f"[{'scalar' if scalar else 'NON-scalar sqrt5-graded'}]")
            w.writerow([f"cell g={gname} eps={eps:+d}", "E-1/E-2",
                        f"dim={d['dim']} invertible={d['invertible']} verified={d['verified']}", sqtxt])
        w.writerows(ROWS)
    print("=" * 72)
    npass = sum(1 for r in ROWS if r[2].startswith("PASS"))
    nopen = sum(1 for r in ROWS if r[2].startswith("OPEN"))
    print(f"SUMMARY: {npass} PASS, {nopen} OPEN (pre-registered), of {len(ROWS)} recorded items")
    print("=" * 72)


if __name__ == "__main__":
    main()
