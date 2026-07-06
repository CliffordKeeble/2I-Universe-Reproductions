"""
w_bench_w1_w3.py — Part D, W1-W3: the torus, its two real forms, and the compact reverser.

W1: J = iZ has J^2=-I, e^{th J}=diag(e^{i th}, e^{-i th}); tau(J)=-J, sigma(J)=+J.
W2: exp(R Z) is the split real form (moduli e^{+-th}); exp(R J) the compact one (moduli 1);
    both sit inside exp(C Z).
W3: under the compact coupling D = d + q A J, entrywise tau maps q-solutions to (-q)-solutions
    with NO matrix part (B = I) — B1.6-4's tau-family reverser, exhibited.

Exact over Q(sqrt5,i); both real places checked where sqrt5 enters. See PRE-REGISTRATION_w1_w3.md.
Ruling 1: outcome-vs-EXPECT only.
"""
import sys
import csv
import sympy as sp
from sympy import I, symbols, Matrix, cos, sin, exp, Abs
from golden_algebra import r5, I2, Z2, G_seed, G_adj, sigma, tau, mcanon, mat_eq

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

ROWS = []
J = I * Z2                              # compact generator iZ = diag(i,-i)


def record(item, expect, passed, outcome, deviation=""):
    tag = "PASS" if passed else "FAIL"
    ROWS.append((item, expect, f"{tag} — {outcome}", deviation))
    print(f"[{tag}] {item}: {outcome}")
    if deviation:
        print(f"        {deviation}")


def csimp(e):
    return sp.simplify(sp.expand(e))


def meq(A, B):
    return all(csimp((A - B)[i]) == 0 for i in range(len(A)))


# ---------------------------------------------------------------- W1
def w1():
    th = symbols('theta', real=True)
    jsq = mat_eq(mcanon(J * J), -I2)
    # e^{th J} via J^2=-I: cos th I + sin th J ; compare to diag(e^{i th}, e^{-i th})
    expJ = cos(th) * I2 + sin(th) * J
    target = Matrix([[exp(I * th), 0], [0, exp(-I * th)]])
    e_ok = meq(expJ.applyfunc(lambda e: e.rewrite(exp)), target)
    tau_ok = mat_eq(tau(J), -J)             # tau(J) = -J
    sig_ok = mat_eq(sigma(J), J)            # sigma(J) = +J (no sqrt5 in J)
    passed = jsq and e_ok and tau_ok and sig_ok
    record("W1 J = iZ properties", "J^2=-I; e^{thJ}=diag(e^{ith},e^{-ith}); tau(J)=-J; sigma(J)=+J",
           passed,
           f"J=iZ: J^2=-I [{jsq}]; e^(thJ)=diag(e^(i th),e^(-i th)) [{e_ok}]; tau(J)=-J "
           f"[{tau_ok}]; sigma(J)=+J [{sig_ok}]")


# ---------------------------------------------------------------- W2
def w2():
    th = symbols('theta', real=True)
    # exp(R Z): split real form
    expZ = Matrix([[exp(th), 0], [0, exp(-th)]])            # e^{th Z}
    split_moduli = (csimp(Abs(expZ[0, 0]) - exp(th)) == 0 and csimp(Abs(expZ[1, 1]) - exp(-th)) == 0)
    split_noncompact = (csimp(Abs(expZ[0, 0]) - 1) != 0)   # modulus not identically 1
    # exp(R J): compact real form
    expJ = Matrix([[exp(I * th), 0], [0, exp(-I * th)]])    # e^{th J} = e^{(i th) Z}
    compact_moduli = (csimp(Abs(expJ[0, 0]) - 1) == 0 and csimp(Abs(expJ[1, 1]) - 1) == 0)
    # both inside exp(C Z): e^{th Z} is z=th; e^{th J} is z=i th  (J = iZ => th J = (i th) Z)
    inside = mat_eq(mcanon(th * J), mcanon((I * th) * Z2))  # th J = (i th) Z
    passed = split_moduli and split_noncompact and compact_moduli and inside
    record("W2 two real forms of exp(C Z)", "exp(RZ) split (moduli e^{+-th}); exp(RJ) compact (1)",
           passed,
           f"exp(R Z) split: eigenvalue moduli e^(+-th), not unit [{split_moduli and split_noncompact}]; "
           f"exp(R J) compact: moduli 1 [{compact_moduli}]; both = exp(z Z) with z real (split) "
           f"or z imaginary (compact) since th J = (i th) Z [{inside}]")


# ---------------------------------------------------------------- W3
def w3():
    k0, k1, q, a0, a1 = symbols('k0 k1 q a0 a1', real=True)

    def Mi(K0, K1, ch):
        """compact / i-coupling  D = d - i q a  (the honest U(1), J=iZ acting)."""
        return G_seed * (-I * (K0 + ch * a0)) + G_adj * (-I * (K1 + ch * a1))

    Mq = mcanon(Mi(k0, k1, q))
    # tau-borne, evaluation at -k (antilinearity, per B1.6-4b); B = I (no matrix part)
    tgt = mcanon(Mi(-k0, -k1, -q))
    reverser = mat_eq(mcanon(tau(Mq)), tgt)     # tau(Mq) at -k equals M_{-q}, with B=I
    # confirm it is within B1.6-4's tau-borne family (dim 2): B=I is a member
    bb = symbols('b0:4')
    B = Matrix(2, 2, bb)
    eqs = list(mcanon(B * mcanon(tau(Mq)) - tgt * B))
    sol = list(sp.linsolve(eqs, list(bb)))[0]
    free = set()
    for c in sol:
        free |= (c.free_symbols & set(bb))
    dim2 = (len(free) == 2)
    # is identity B=I in the solution space? (b0=b3=1, b1=b2=0 consistent)
    idsol = all(csimp(sol[i].subs({bb[0]: 1, bb[3]: 1, bb[1]: 0, bb[2]: 0})
                      - [1, 0, 0, 1][i]) == 0 for i in range(4))
    passed = reverser and dim2 and idsol
    record("W3 compact-coupling reverser", "tau maps q->-q with no matrix part (B=I)", passed,
           f"tau(Mq) at -k = M_{{-q}} exactly [{reverser}] — entrywise tau (B=I, no matrix part) "
           f"reverses the compact charge; and B=I lies in B1.6-4's tau-borne family (dim "
           f"{len(free)}) [{idsol}]")


def main():
    print("=" * 72)
    print("PART D — W1-W3 (torus / two real forms / compact reverser), exact over Q(sqrt5,i)")
    print("=" * 72)
    for fn in (w1, w2, w3):
        fn()
    with open("results_w1_w3.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["item", "expect", "outcome", "deviation"])
        w.writerows(ROWS)
    npass = sum(1 for row in ROWS if row[2].startswith("PASS"))
    print("=" * 72)
    print(f"SUMMARY: {npass}/{len(ROWS)} match EXPECT")
    print("=" * 72)


if __name__ == "__main__":
    main()
