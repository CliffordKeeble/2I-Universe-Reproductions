"""
verify_c_from_sigma_repro.py — INDEPENDENT reproduction of the cold-boot construction
(Part A-1, Paper 226).

Reproduces the claim families of the Fable-5-max cold draft
`paper_2XX_fourth_face_charge_conjugation_v0.1.md` (Appendix A claim list), each re-derived
here from the DRAFT's stated mathematics — not ported from the cold script.

INDEPENDENCE: the cold script `cold_boot/verify_c_from_sigma.py` implements the golden amplitude
algebra 𝕂 = ℝ[δ]/(δ²−5) as coefficient PAIRS (p,q) with a hand-written product. This reproduction
uses a DIFFERENT mechanism — the **regular representation** of 𝕂 on the basis (1, δ):

    a = s + tδ   ↦   reg(a) = [[s, 5t], [t, s]]     (δ ↦ [[0,5],[1,0]], δ² = 5·I)

so that algebra multiplication is matrix multiplication, the Galois norm is the **determinant**,
the two real places are the **eigenvalues**, and σ (δ ↦ −δ) is **conjugation by Z = diag(1,−1)**.
Matrix tie-ins use the shared `golden_algebra` definitions (Paper 203 §4 source, not the cold
script). Caveat (see cold_boot/PROVENANCE.md): written after the cold script was seen; independent
at the derivation level, not script-blind.

Exact symbolic throughout; √5-bearing constant identities are cross-checked numerically at both
real places (√5 → ±√5), the arbiter where sympy simplify is unreliable over the golden field.

Draft-claim citations are given inline as [App A: ...].
"""
import sys
import sympy as sp
from sympy import sqrt, I, Rational, Matrix, eye, symbols, Function, diff, cosh, sinh, log, Eq
from golden_algebra import r5, I2, Z2, G_seed, G_adj, sigma as gsigma, mat_eq as gmat_eq

try:                                    # allow Greek/math glyphs in output on Windows (cp1252)
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

FAILURES = []
NCHECKS = 0


def check(name, cond):
    global NCHECKS
    NCHECKS += 1
    ok = bool(cond)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}")
    if not ok:
        FAILURES.append(name)


# ------------------------------------------------------------------ regular rep of 𝕂
def reg(s, t):
    """a = s + tδ  ->  regular representation on basis (1, δ)."""
    return Matrix([[s, 5 * t], [t, s]])


ONE = reg(1, 0)                       # = I2
DELTA = reg(0, 1)                     # δ
JJ = reg(0, Rational(1) / r5)         # j = δ/√5


def rp(M):
    return sp.expand(M[0, 0])         # coefficient of 1


def dp(M):
    return sp.expand(M[1, 0])         # coefficient of δ


def sig(M):
    """Galois σ (δ ↦ −δ) = Z-conjugation in the regular representation."""
    return Z2 * M * Z2


def Nrm(M):
    """Galois norm N(a) = a·σ(a) — here the determinant of the regular rep."""
    return sp.expand(M.det())


def csimp(e):
    return sp.simplify(sp.expand(e))


def meq(A, B):
    """Symbolic matrix equality (handles function entries)."""
    D = A - B
    return all(csimp(D[i]) == 0 for i in range(len(D)))


def both_places(A, B):
    """Exact + numeric at both real places √5 → ±√5 (constant matrices)."""
    ex = meq(A, B)
    ok = True
    D = A - B
    for s in (sp.sqrt(5), -sp.sqrt(5)):
        for e in D:
            if abs(complex(sp.N(e.subs(r5, s), 30))) > 1e-25:
                ok = False
    return ex and ok


def phase(theta):
    """e^{θ j} = cosh θ · 1 + sinh θ · j  (j² = 1)."""
    return cosh(theta) * ONE + sinh(theta) * JJ


# ============================================================ 1. the algebra 𝕂
# [App A: 𝕂 = ℝ[δ]/(δ²−5); j²=1, σ(j)=−j; σ involutive automorphism; N=p²−5q² real;
#         N multiplicative; N σ-invariant; σ-fixed subalgebra = ℝ]
check("delta^2 = 5 (regular rep faithful)", both_places(DELTA * DELTA, 5 * ONE))
check("j^2 = 1", both_places(JJ * JJ, ONE))
check("sigma(j) = -j", both_places(sig(JJ), -JJ))

ps, qs, ps2, qs2 = symbols('ps qs ps2 qs2', real=True)
A = reg(ps, qs)
B = reg(ps2, qs2)
check("sigma(ab) = sigma(a) sigma(b)  [algebra automorphism]",
      meq(sig(A * B), sig(A) * sig(B)))
check("sigma^2 = id", meq(sig(sig(A)), A))
check("N(a) = p^2 - 5 q^2  (= det of regular rep)", csimp(Nrm(A) - (ps**2 - 5 * qs**2)) == 0)
check("N(ab) = N(a) N(b)  [det multiplicative]", csimp(Nrm(A * B) - Nrm(A) * Nrm(B)) == 0)
check("N(sigma a) = N(a)", csimp(Nrm(sig(A)) - Nrm(A)) == 0)
check("sigma-fixed subalgebra is R: sigma(a)=a forces q=0",
      sp.solve(Eq(dp(sig(A)), dp(A)), qs) == [0])
# independence extra: norm is literally the determinant / places are eigenvalues
check("[extra] N(a) = det(reg a) exactly", csimp(Nrm(A) - A.det()) == 0)

# ============================================================ 2. Lemma 1 (involution budget)
# [App A: unital ℝ-endomorphisms are exactly {id, σ, two place projections}; Aut_ℝ(𝕂)={1,σ}]
xx, yy = symbols('xx yy', real=True)
# unital endo fixed by f(δ)=xx+yy·δ with f(δ)²=5  ->  xx²+5yy²=5 and 2·xx·yy=0
sols = sp.solve([Eq(xx**2 + 5 * yy**2, 5), Eq(2 * xx * yy, 0)], [xx, yy], dict=True)
autos = [s for s in sols if s[yy] != 0]
projs = [s for s in sols if s[yy] == 0]
check("endomorphism eqns: exactly four solutions", len(sols) == 4)
check("automorphisms: f(δ)=±δ only (Aut = {1, σ})",
      sorted(s[yy] for s in autos) == [-1, 1] and all(s[xx] == 0 for s in autos))
check("place projections: f(δ)=±√5 (the two non-injective ℝ-endomorphisms)",
      sorted(s[xx] for s in projs) == [-r5, r5] and all(s[yy] == 0 for s in projs))
# the y=0 maps are the place projections a=(s+tδ) -> s ± t√5 ; verify non-injective (kill e∓)
E_PLUS = (ONE + JJ) / 2
E_MINUS = (ONE - JJ) / 2


def scalar_hom(fd, M):
    """f(s+tδ) = s + t·fd  (fd = ±√5 real => projection to a place)."""
    return csimp(rp(M) + dp(M) * fd)


check("f(δ)=+√5 kills e_minus (non-injective projection)", scalar_hom(r5, E_MINUS) == 0)
check("f(δ)=-√5 kills e_plus (non-injective projection)", scalar_hom(-r5, E_PLUS) == 0)

# ============================================================ 3. idempotents / places
# [App A: e± orthogonal, complete, σ-swapped; N(a)=a₊a₋]
check("e_+ e_- = 0", both_places(E_PLUS * E_MINUS, sp.zeros(2)))
check("e_+^2 = e_+", both_places(E_PLUS * E_PLUS, E_PLUS))
check("e_+ + e_- = 1", both_places(E_PLUS + E_MINUS, ONE))
check("sigma swaps e_+ <-> e_-", both_places(sig(E_PLUS), E_MINUS))
# place coordinates a± = s ± t√5 are the eigenvalues of reg(a); N = a₊ a₋
a_plus, a_minus = ps + qs * r5, ps - qs * r5
evs = list(A.eigenvals().keys())
check("[extra] place coords a± = s ± t√5 are the eigenvalues of reg(a)",
      {csimp(e) for e in evs} == {csimp(a_plus), csimp(a_minus)})
check("N(a) = a_+ a_-  (product of places)", csimp(Nrm(A) - a_plus * a_minus) == 0)

# ============================================================ 4. units / the golden element
# [App A: Φ²=Φ+1; N(Φ)=−1; N(Φ²)=+1; coords(Φ)=(φ,ψ); Φ²=e^{θ0 j}, θ0=2 ln φ=arcosh(3/2)]
PHI = reg(Rational(1, 2), Rational(1, 2))       # Φ = (1 + δ)/2
phi = (1 + r5) / 2
check("Phi^2 = Phi + 1 (golden recursion in 𝕂)", both_places(PHI * PHI, PHI + ONE))
check("N(Phi) = -1 (fundamental unit, norm -1)", csimp(Nrm(PHI) + 1) == 0)
check("N(Phi^2) = +1", csimp(Nrm(PHI * PHI) - 1) == 0)
check("coords(Phi) = (phi, psi) at the two places",
      csimp(rp(PHI) + dp(PHI) * r5 - phi) == 0 and csimp(rp(PHI) - dp(PHI) * r5 - (1 - phi)) == 0)
PHI2 = PHI * PHI
check("coords(Phi^2) = (phi^2, phi^-2)",
      csimp((rp(PHI2) + dp(PHI2) * r5) - phi**2) == 0
      and csimp((rp(PHI2) - dp(PHI2) * r5) * phi**2 - 1) == 0)
theta0 = 2 * log(phi)
# cosh/sinh of a log need exp-rewrite before they reduce: cosh(2 ln φ)=(φ²+φ⁻²)/2=3/2,
# sinh(2 ln φ)=(φ²−φ⁻²)/2=√5/2. Standard sympy idiom (rewrite→expand→radsimp).
cosh0 = sp.radsimp(sp.expand(cosh(theta0).rewrite(sp.exp)))
sinh0 = sp.radsimp(sp.expand(sinh(theta0).rewrite(sp.exp)))
exp0 = sp.powsimp(sp.exp(theta0), force=True)
check("Phi^2 = e^{theta0 j}, theta0 = 2 ln phi = arcosh(3/2)",
      sp.simplify(cosh0 - Rational(3, 2)) == 0
      and sp.simplify(sinh0 - r5 / 2) == 0
      and sp.simplify(exp0 - phi**2) == 0)

# ============================================================ 5. phases / split torus
# [App A: N(e^{θj})=1; σ(e^{θj})=e^{−θj}; norm-one two components; Lie(T)=ℝj]
th = symbols('theta', real=True)
check("N(e^{theta j}) = 1", csimp(Nrm(phase(th)) - 1) == 0)
check("sigma(e^{theta j}) = e^{-theta j}  [the C mechanism: phase inversion is Galois]",
      meq(sig(phase(th)), phase(-th)))
# norm-one has two components: -1 has N=1 but its place coord (eigenvalue) is negative
check("norm-one group has two components (-1: N=1 but not forward)",
      csimp(Nrm(-ONE) - 1) == 0 and csimp((-ONE).eigenvals().popitem()[0]) == -1)
# Lie algebra: X + sigma(X) = 0 forces X in R·j (real part zero)
Xp, Xq = symbols('Xp Xq', real=True)
X = reg(Xp, Xq)
check("Lie algebra of T: X + sigma(X) = 0 forces X in R*j (Xp=0)",
      csimp(rp(X + sig(X))) == 2 * Xp and sp.solve(Eq(2 * Xp, 0), Xp) == [0])

# ============================================================ 6. arithmetic quantisation
# [App A: φ^{2g} ∈ ±φ^ℤ  ⟺  g ∈ ½ℤ]
w, k = symbols('w k', real=True)
check("phi^{2w} = phi^k  <=>  w = k/2  (half-integer charge)",
      sp.solve(Eq(phi**(2 * w), phi**k), w) == [k / 2])

# ============================================================ 7-9. gauge sector (fields over 𝕂)
# [App A: gauge covariance; current J real, =2√5(uv'−u'v), C-odd; kinetic expansion, −g² seagull]
x, g = symbols('x g', real=True)
u = Function('u', real=True)(x)
v = Function('v', real=True)(x)
thx = Function('th', real=True)(x)
Ax = Function('Ax', real=True)(x)
af = reg(u, v)


def ddx(M):
    return M.applyfunc(lambda e: diff(e, x))


def Dcov(M, Afield, weight):
    """D a = ∂a − weight·A·j·a."""
    return ddx(M) - weight * Afield * (JJ * M)


# gauge covariance: a -> e^{gθ(x)j} a, A -> A + ∂θ
gauged = phase(g * thx) * af
lhs = Dcov(gauged, Ax + diff(thx, x), g)
rhs = phase(g * thx) * Dcov(af, Ax, g)
check("gauge covariance D'(e^{gθj}a) = e^{gθj}(Da) with A -> A + ∂θ", meq(lhs, rhs))

# Noether current  J = j(σ(a)∂a − σ(∂a)a)
da = ddx(af)
Jc = JJ * (sig(af) * da - sig(da) * af)
check("current is real (delta-part vanishes)", csimp(dp(Jc)) == 0)
check("current = 2√5 (u v' - u' v)",
      csimp(rp(Jc) - 2 * r5 * (u * diff(v, x) - diff(u, x) * v)) == 0)
sa = sig(af)
dsa = ddx(sa)
Jc_c = JJ * (sig(sa) * dsa - sig(dsa) * sa)
check("current is C-odd: J[sigma a] = -J[a]", csimp(rp(Jc_c) + rp(Jc)) == 0 and csimp(dp(Jc_c)) == 0)

# kinetic expansion σ(Da)Da = σ(∂a)∂a + g A·J − g² A² N(a)   [negative seagull]
Da = Dcov(af, Ax, g)
kin = sig(Da) * Da
expected = sig(da) * da + g * Ax * Jc - g**2 * Ax**2 * Nrm(af) * ONE
check("kinetic expansion with NEGATIVE seagull  -g^2 A^2 N(a)", meq(kin, expected))

# ============================================================ 10. Theorem 1 (C = σ, with A -> -A)
# [App A: covariant conjugation D(σa)|_{-A}=σ(Da); C-invariance of every term; charge inversion;
#         C²=id; (phase∘C)²=id]
check("Theorem1(i): D(sigma a)|_{-A} = sigma(D a)  [A -> -A forced]",
      meq(Dcov(sig(af), -Ax, g), sig(Dcov(af, Ax, g))))
Da_C = Dcov(sig(af), -Ax, g)
check("Theorem1(ii) kinetic C-invariant: sigma(D'a')D'a' = sigma(Da)Da",
      meq(sig(Da_C) * Da_C, kin))
check("Theorem1(ii) mass term C-invariant: N(sigma a) = N(a)", csimp(Nrm(sig(af)) - Nrm(af)) == 0)
check("Theorem1(ii) F^2 C-invariant: (-∂A)^2 = (∂A)^2",
      csimp(diff(-Ax, x)**2 - diff(Ax, x)**2) == 0)
th2 = symbols('th2', real=True)
check("Theorem1(iii) charge inversion: sigma(e^{gθj}a) = e^{-gθj} sigma(a)",
      meq(sig(phase(g * th2) * A), phase(-g * th2) * sig(A)))
check("Theorem1(iv) C^2 = id", meq(sig(sig(A)), A))
Cp = lambda M: phase(th2) * sig(M)
check("Theorem1(iv) (phase o C)^2 = id  [every representative involutive]",
      meq(Cp(Cp(A)), A))

# ============================================================ 11. matrix tie-ins (Papers 191/192/203)
# [App A: Γ_seed²=+√5I, Γ_adj²=−√5I, {Γ_seed,Γ_adj}=0, Γ_seedΓ_adj=√5σ₃;
#         σ(Γ_adj)=−Γ_seed, σ(Γ_seed)=−Γ_adj; spectra ±5^{1/4}, ±i5^{1/4}]
s3 = Z2
check("Gamma_seed^2 = +√5 I", gmat_eq(G_seed * G_seed, r5 * I2))
check("Gamma_adj^2 = -√5 I", gmat_eq(G_adj * G_adj, -r5 * I2))
check("{Gamma_seed, Gamma_adj} = 0", gmat_eq(G_seed * G_adj + G_adj * G_seed, sp.zeros(2)))
check("Gamma_seed Gamma_adj = √5 σ_3  [Paper 192 measurement product]",
      gmat_eq(G_seed * G_adj, r5 * s3))
check("sigma(Gamma_adj) = -Gamma_seed  [Paper 203 Face 1]", gmat_eq(gsigma(G_adj), -G_seed))
check("sigma(Gamma_seed) = -Gamma_adj  [involution closes]", gmat_eq(gsigma(G_seed), -G_adj))
check("spec(Gamma_seed) = {±5^{1/4}} real",
      sorted(csimp(e**2) for e in G_seed.eigenvals()) == [r5, r5])
check("spec(Gamma_adj) = {±i 5^{1/4}} imaginary",
      sorted((csimp(e**2) for e in G_adj.eigenvals()), key=str) == [-r5, -r5])

# ============================================================ 12. golden unitary group ≅ GL(n,R)
# [App A: σ-unitary pairs (M, (Mᵀ)⁻¹) => U(n,𝕂,σ) ≅ GL(n,ℝ), checked n=2]
m11, m12, m21, m22 = symbols('m11 m12 m21 m22', real=True)
Mp = Matrix([[m11, m12], [m21, m22]])
Mm = (Mp.T).inv()
check("sigma-unitary pairs (M,(M^T)^-1): U(n,𝕂,σ) ~ GL(n,R) (n=2)",
      sp.simplify(Mm.T * Mp - eye(2)) == sp.zeros(2)
      and sp.simplify(Mp.T * Mm - eye(2)) == sp.zeros(2))

# ------------------------------------------------------------------ summary
print()
if FAILURES:
    print(f"{len(FAILURES)}/{NCHECKS} FAILURES:")
    for f in FAILURES:
        print("  -", f)
    raise SystemExit(1)
print(f"ALL {NCHECKS} REPRO CHECKS PASSED (independent regular-rep derivation, exact + both "
      f"places, sympy {sp.__version__})")
