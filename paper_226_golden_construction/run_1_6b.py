"""
run_1_6b.py — Part B (B1.6-4b): the factorisation frontier re-run with the evaluation
point DERIVED from antilinearity (see B1_6_4b_derivation.md), not chosen.

The momentum sign of each borne type is obtained by applying that type's field automorphism
to the plane-wave factor e^{i k x} and reading off the resulting momentum:
    linear  -> e^{i k x}      (momentum +k)
    sigma   -> e^{i k x}      (sigma fixes i; momentum +k)
    tau     -> e^{-i k x}     (tau conjugates i; momentum -k)
The target operator L_{-q} is then evaluated at (sign)*(k0,k1). No -k is hard-coded.

Reproduces Run #1.6's verdicts:  dims 0/0/0 over K (sigma-coupling);  tau-borne dim 2
over K(i) (i-coupling).  Exact over Q(sqrt5,i). One claim per EXPECT.
"""
import sys
import csv
import sympy as sp
from sympy import I, symbols, exp, Matrix
from golden_algebra import r5, I2, Z2, G_seed, G_adj, sigma, tau, mcanon

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

k0, k1, q, a0, a1 = symbols('k0 k1 q a0 a1', real=True)
ROWS = []


def record(item, expect, passed, outcome, deviation=""):
    tag = "PASS" if passed else "FAIL"
    ROWS.append((item, expect, f"{tag} — {outcome}", deviation))
    print(f"[{tag}] {item}: {outcome}")
    if deviation:
        print(f"        {deviation}")


# ------------------------------------------------------------------ couplings
def Msig(K0, K1, ch):
    """golden sigma-coupling over K:  D = d - q a Z  (Z-weight, carries sqrt5)."""
    return G_seed * ((-I * K0) * I2 - ch * a0 * Z2) + G_adj * ((-I * K1) * I2 - ch * a1 * Z2)


def Mi(K0, K1, ch):
    """honest U(1) i-coupling over K(i):  D = d - i q a  (T_tau, the extension)."""
    return G_seed * (-I * (K0 + ch * a0)) + G_adj * (-I * (K1 + ch * a1))


# ------------------------------------------------------------------ derived evaluation point
def momentum_sign(borne):
    """Apply the borne field-map to e^{i k x} and read the resulting momentum sign.
       Antilinearity (i -> -i) is the ONLY thing that flips it. Returns +1 or -1."""
    x, k = symbols('x k', real=True)
    pw = exp(I * k * x)
    if borne == 'lin':
        mapped = pw
    elif borne == 'tau':
        mapped = sp.conjugate(pw)          # tau: i -> -i  =>  e^{-i k x}
    elif borne == 'sigma':
        mapped = pw.subs(r5, -r5)          # sigma: sqrt5 -> -sqrt5, fixes i  =>  unchanged
    else:
        raise ValueError(borne)
    kprime = sp.simplify(sp.diff(mapped, x) / (I * mapped))   # = +k or -k
    return sp.simplify(kprime / k)


# ------------------------------------------------------------------ frontier dimension
def dim_C(Mfunc, borne):
    """dim of constant 2x2 B giving a `borne` charge conjugation L_q -> L_{-q},
       with the target momentum DERIVED via momentum_sign(borne)."""
    s = momentum_sign(borne)
    Mq = mcanon(Mfunc(k0, k1, q))
    fieldmap = {'lin': (lambda X: X), 'tau': tau, 'sigma': sigma}[borne]
    src = mcanon(fieldmap(Mq))
    tgt = mcanon(Mfunc(s * k0, s * k1, -q))      # evaluated at derived momentum
    bb = symbols('b0:4')
    B = Matrix(2, 2, bb)
    eqs = list(mcanon(B * src - tgt * B))
    sol = list(sp.linsolve(eqs, list(bb)))[0]
    fr = set()
    for c in sol:
        fr |= (c.free_symbols & set(bb))
    return len(fr)


def main():
    print("=" * 72)
    print("RUN #1.6b — B1.6-4b factorisation frontier, evaluation point DERIVED (exact over Q(sqrt5,i))")
    print("=" * 72)

    # B1.6-4b-i : the evaluation points, derived from the borne maps
    s_lin, s_sig, s_tau = int(momentum_sign('lin')), int(momentum_sign('sigma')), int(momentum_sign('tau'))
    derived_ok = (s_lin == 1 and s_sig == 1 and s_tau == -1)
    record("B1.6-4b-i evaluation point derived",
           "tau -> -k (antilinear); lin, sigma -> +k",
           derived_ok,
           f"momentum signs from applying each map to e^(ikx): lin={s_lin:+d}k, "
           f"sigma={s_sig:+d}k, tau={s_tau:+d}k — tau(e^(ikx))=e^(-ikx) gives -k by antilinearity",
           "" if derived_ok else "derived signs unexpected — antilinearity derivation failed")

    # B1.6-4b-ii : the frontier, with those derived points
    sig_row = {b: dim_C(Msig, b) for b in ("lin", "tau", "sigma")}
    icp_row = {b: dim_C(Mi, b) for b in ("lin", "tau", "sigma")}
    K_zero = (sig_row["lin"] == 0 and sig_row["tau"] == 0 and sig_row["sigma"] == 0)
    Ki_tau2 = (icp_row["tau"] == 2)
    reproduces = K_zero and Ki_tau2
    record("B1.6-4b-ii frontier reproduces run #1.6",
           "dims 0/0/0 over K; tau-borne dim 2 over K(i)",
           reproduces,
           f"sigma-coupling (K): lin={sig_row['lin']}, tau={sig_row['tau']}, sigma={sig_row['sigma']} "
           f"=> 0/0/0 [{K_zero}]; i-coupling (K(i)): lin={icp_row['lin']}, tau={icp_row['tau']}, "
           f"sigma={icp_row['sigma']} => tau-borne dim 2 [{Ki_tau2}]. Verdicts match run #1.6; "
           f"the -k was produced by antilinearity, not chosen",
           "" if reproduces else "DEVIATION from run #1.6 — the in-flight correction may have been "
           "fitted; asterisk becomes a strike")

    with open("results_run1_6b.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["item", "expect", "outcome", "deviation"])
        w.writerows(ROWS)
    npass = sum(1 for row in ROWS if row[2].startswith("PASS"))
    print("=" * 72)
    print(f"SUMMARY: {npass}/{len(ROWS)} match EXPECT")
    print("=" * 72)


if __name__ == "__main__":
    main()
