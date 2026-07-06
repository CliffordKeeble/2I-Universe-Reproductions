#!/usr/bin/env python3
"""
Cold-boot verification: charge conjugation from the minimal interacting
extension of Paper 192's golden QM.

Amplitude algebra: KK = R[d]/(d^2 - 5), elements a = p + q*d represented
as pairs (p, q). sigma: d -> -d. Norm N(a) = a*sigma(a). j := d/sqrt(5),
j^2 = 1. Idempotent coordinates a_plus = p + q*sqrt5, a_minus = p - q*sqrt5.

Exact symbolic throughout (sympy). Each check prints PASS/FAIL.
"""

import sympy as sp

sqrt5 = sp.sqrt(5)
x, theta, g, m, A = sp.symbols('x theta g m A', real=True)

FAILURES = []

def check(name, cond):
    ok = bool(cond)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}")
    if not ok:
        FAILURES.append(name)

# ---------- the split golden algebra KK ----------

def mul(a, b):
    """(p1 + q1 d)(p2 + q2 d) with d^2 = 5."""
    p1, q1 = a; p2, q2 = b
    return (sp.expand(p1*p2 + 5*q1*q2), sp.expand(p1*q2 + q1*p2))

def add(a, b):
    return (sp.expand(a[0] + b[0]), sp.expand(a[1] + b[1]))

def sub(a, b):
    return (sp.expand(a[0] - b[0]), sp.expand(a[1] - b[1]))

def smul(c, a):
    return (sp.expand(c*a[0]), sp.expand(c*a[1]))

def sigma(a):
    return (a[0], -a[1])

def N(a):
    """Galois norm N(a) = a*sigma(a); returns KK element (real part, d part)."""
    return mul(a, sigma(a))

def coords(a):
    """Idempotent coordinates (a_plus, a_minus) at the two real places."""
    p, q = a
    return (sp.expand(p + q*sqrt5), sp.expand(p - q*sqrt5))

def simp(a):
    return (sp.simplify(a[0]), sp.simplify(a[1]))

ONE = (sp.Integer(1), sp.Integer(0))
D   = (sp.Integer(0), sp.Integer(1))            # d, the formal sqrt(5)
J   = (sp.Integer(0), sp.Rational(1, 1)/sqrt5)  # j = d/sqrt5

# j^2 = 1, sigma(j) = -j
check("j^2 = 1", simp(mul(J, J)) == ONE)
check("sigma(j) = -j", sigma(J) == smul(-1, J))

# sigma is an involutive algebra automorphism; N is multiplicative, sigma-invariant
pa, qa, pb, qb = sp.symbols('pa qa pb qb', real=True)
a_gen = (pa, qa); b_gen = (pb, qb)
check("sigma(ab) = sigma(a)sigma(b)",
      simp(sub(sigma(mul(a_gen, b_gen)), mul(sigma(a_gen), sigma(b_gen)))) == (0, 0))
check("sigma^2 = id", sigma(sigma(a_gen)) == a_gen)
check("N(a) = p^2 - 5 q^2 (real, i.e. d-part zero)",
      simp(N(a_gen)) == (sp.expand(pa**2 - 5*qa**2), 0))
check("N(ab) = N(a) N(b)",
      sp.simplify(N(mul(a_gen, b_gen))[0] - N(a_gen)[0]*N(b_gen)[0]) == 0)
check("N(sigma(a)) = N(a)", simp(sub(N(sigma(a_gen)), N(a_gen))) == (0, 0))
check("sigma-fixed subalgebra is the reals: sigma(a) = a forces q = 0",
      sp.solve(sp.Eq(-qa, qa), qa) == [0])

# ---------- Lemma: Aut_R(KK) = {id, sigma} ----------
# unital R-algebra endomorphism: f(d) = xx + yy*d with f(d)^2 = 5
xx, yy = sp.symbols('xx yy', real=True)
fd = (xx, yy)
sq = mul(fd, fd)
sols = sp.solve([sp.Eq(sq[0], 5), sp.Eq(sq[1], 0)], [xx, yy], dict=True)
auto_sols = [s for s in sols if s[yy] != 0]      # invertible (automorphisms)
proj_sols = [s for s in sols if s[yy] == 0]      # non-injective projections
check("endomorphism eqns: exactly the four expected solutions", len(sols) == 4)
check("automorphisms: f(d) = +d or f(d) = -d only",
      sorted([s[yy] for s in auto_sols]) == [-1, 1]
      and all(s[xx] == 0 for s in auto_sols))
# the yy = 0 solutions are the two place-projections, with nontrivial kernel:
# f(d) = +sqrt5 kills e_minus = (1 - j)/2 ; f(d) = -sqrt5 kills e_plus
E_PLUS  = smul(sp.Rational(1, 2), add(ONE, J))
E_MINUS = smul(sp.Rational(1, 2), sub(ONE, J))
def apply_scalar_hom(fd_val, a):
    """f(p + q d) = p + q * fd_val  (fd_val real => projection to R)."""
    p, q = a
    return sp.simplify(p + q*fd_val)
check("f(d)=+sqrt5 is the projection killing e_minus (non-injective)",
      apply_scalar_hom(sqrt5, E_MINUS) == 0)
check("f(d)=-sqrt5 is the projection killing e_plus (non-injective)",
      apply_scalar_hom(-sqrt5, E_PLUS) == 0)

# idempotents and coordinates
check("e_+ e_- = 0", simp(mul(E_PLUS, E_MINUS)) == (0, 0))
check("e_+^2 = e_+", simp(sub(mul(E_PLUS, E_PLUS), E_PLUS)) == (0, 0))
check("e_+ + e_- = 1", simp(add(E_PLUS, E_MINUS)) == ONE)
check("sigma swaps e_+ and e_-", simp(sub(sigma(E_PLUS), E_MINUS)) == (0, 0))
check("N in coordinates is a_plus * a_minus",
      sp.simplify(N(a_gen)[0] - coords(a_gen)[0]*coords(a_gen)[1]) == 0)

# ---------- units: the golden element inside KK ----------
PHI = (sp.Rational(1, 2), sp.Rational(1, 2))      # (1 + d)/2
phi = (1 + sqrt5)/2
check("Phi^2 = Phi + 1 (golden recursion inside KK)",
      simp(sub(mul(PHI, PHI), add(PHI, ONE))) == (0, 0))
check("N(Phi) = -1 (fundamental unit has norm -1)", sp.simplify(N(PHI)[0] + 1) == 0)
check("N(Phi^2) = +1", sp.simplify(N(mul(PHI, PHI))[0] - 1) == 0)
cp, cm = coords(PHI)
check("coords(Phi) = (phi, psi): the two golden roots at the two places",
      sp.simplify(cp - phi) == 0 and sp.simplify(cm - (1 - phi)) == 0)
# Phi^2 as hyperbolic phase e^{theta0 j}: coords (e^t0, e^-t0), t0 = 2 ln phi
PHI2 = mul(PHI, PHI)
c2p, c2m = coords(PHI2)
check("coords(Phi^2) = (phi^2, phi^-2)",
      sp.simplify(c2p - phi**2) == 0 and sp.simplify(c2m - phi**(-2)) == 0)
t0 = 2*sp.log(phi)
cosh_t0 = sp.radsimp(sp.expand(sp.cosh(t0).rewrite(sp.exp)))
sinh_t0 = sp.radsimp(sp.expand(sp.sinh(t0).rewrite(sp.exp)))
check("Phi^2 = e^{theta0 j} with theta0 = 2 ln phi = arcosh(3/2)",
      sp.simplify(cosh_t0 - sp.Rational(3, 2)) == 0
      and sp.simplify(sinh_t0 - sqrt5/2) == 0
      and sp.simplify(sp.expand(sp.exp(t0)) - phi**2) == 0)

# hyperbolic phase e^{theta j} = cosh theta + j sinh theta ; norm-one, forward
def phase(t):
    return add(smul(sp.cosh(t), ONE), smul(sp.sinh(t), J))
check("N(e^{theta j}) = 1", sp.simplify(N(phase(theta))[0] - 1) == 0)
check("sigma(e^{theta j}) = e^{-theta j}  [phase inversion: the C mechanism]",
      simp(sub(sigma(phase(theta)), phase(-theta))) == (0, 0))
check("norm-one group has two components: N(-1)=1 and -1 not a forward phase",
      sp.simplify(N(smul(-1, ONE))[0] - 1) == 0
      and sp.simplify(coords(smul(-1, ONE))[0]) == -1)  # e_+ coord negative

# Lie algebra of the norm-one group is R*j  (sigma-odd elements)
Xp, Xq = sp.symbols('Xp Xq', real=True)
X = (Xp, Xq)
lin = simp(add(X, sigma(X)))       # N(1+eps X) = 1 + eps(X + sigma X) + O(eps^2)
check("Lie algebra condition X + sigma(X) = 0 forces X in R*j (Xp = 0)",
      lin == (sp.expand(2*Xp), 0) and sp.solve(sp.Eq(2*Xp, 0), Xp) == [0])

# arithmetic quantisation: phi^{2w} in +-phi^Z  <=>  2w in Z
w, k = sp.symbols('w k', real=True)
sol_w = sp.solve(sp.Eq(phi**(2*w), phi**k), w)
check("phi^{2w} = phi^k  <=>  w = k/2  (charge quantised in half-integers)",
      sol_w == [k/2])

# ---------- gauge structure: covariance, current, C-invariance ----------
u = sp.Function('u', real=True)(x)
v = sp.Function('v', real=True)(x)
th = sp.Function('th', real=True)(x)
Af = sp.Function('Af', real=True)(x)
af = (u, v)

def d_dx(a):
    return (sp.diff(a[0], x), sp.diff(a[1], x))

def Dcov(a, Afield, weight):
    """D a = da/dx - weight * A * j * a"""
    return simp(sub(d_dx(a), smul(weight*Afield, mul(J, a))))

# local phase with weight g: a -> e^{g th(x) j} a ;  A -> A + dth/dx
gauged_a = simp(mul(phase(g*th), af))
lhs = Dcov(gauged_a, Af + sp.diff(th, x), g)
rhs = simp(mul(phase(g*th), Dcov(af, Af, g)))
check("gauge covariance: D'(e^{g th j} a) = e^{g th j} (D a) with A -> A + dth",
      simp(sub(lhs, rhs)) == (0, 0))

# Noether current  Jc = j (sigma(a) da - sigma(da) a) : real and sigma-odd
da = d_dx(af)
Jc = simp(mul(J, sub(mul(sigma(af), da), mul(sigma(da), af))))
check("current is real (d-part zero)", sp.simplify(Jc[1]) == 0)
Jc_real = sp.simplify(Jc[0])
check("current = 2 sqrt5 (u v' - u' v)",
      sp.simplify(Jc_real - 2*sqrt5*(u*sp.diff(v, x) - sp.diff(u, x)*v)) == 0)
sa = sigma(af); dsa = d_dx(sa)
Jc_conj = simp(mul(J, sub(mul(sigma(sa), dsa), mul(sigma(dsa), sa))))
check("current is C-odd: J[sigma(a)] = -J[a]",
      sp.simplify(Jc_conj[0] + Jc_real) == 0 and sp.simplify(Jc_conj[1]) == 0)

# kinetic expansion:  sigma(Da) Da = sigma(da) da + g A Jc - g^2 A^2 N(a)
Da = Dcov(af, Af, g)
kin = simp(mul(sigma(Da), Da))
free_kin = mul(sigma(da), da)
expected = add(add(free_kin, smul(g*Af, Jc)), smul(-g**2*Af**2, N(af)))
check("expansion sigma(Da)Da = sigma(da)da + gA*J - g^2 A^2 N(a) "
      "[seagull sign is NEGATIVE: structural difference from complex scalar QED]",
      simp(sub(kin, expected)) == (0, 0))
# term-by-term C-invariance: C: a -> sigma(a), A -> -A
Da_C = Dcov(sigma(af), -Af, g)
kin_C = simp(mul(sigma(Da_C), Da_C))
check("C-invariance of full gauged kinetic term: sigma(D'a')D'a' = sigma(Da)Da",
      simp(sub(kin_C, kin)) == (0, 0))
check("D(sigma a) with A -> -A equals sigma(D a)  [covariant conjugation]",
      simp(sub(Da_C, sigma(Da))) == (0, 0))
check("mass term C-invariant: N(sigma a) = N(a)",
      simp(sub(N(sigma(af)), N(af))) == (0, 0))
# field strength: F = dA in 1+0 toy reduces to derivative; sign flip squares away
check("F^2 C-invariant: (-dA/dx)^2 = (dA/dx)^2",
      sp.simplify((sp.diff(-Af, x))**2 - (sp.diff(Af, x))**2) == 0)

# weight inversion: sigma of a weight-g transform is a weight-(-g) transform
lhsw = sigma(simp(mul(phase(g*theta), a_gen)))
rhsw = simp(mul(phase(-g*theta), sigma(a_gen)))
check("sigma(e^{g theta j} a) = e^{-g theta j} sigma(a)  [charge g -> -g]",
      simp(sub(lhsw, rhsw)) == (0, 0))

# C is involutive, and remains involutive after composing with any phase
check("C^2 = id on amplitudes", sigma(sigma(a_gen)) == a_gen)
Cp = lambda a: simp(mul(phase(theta), sigma(a)))     # C' = phase o C
check("(phase o C)^2 = id  [every representative of the sigma-class is involutive]",
      simp(sub(Cp(Cp(a_gen)), a_gen)) == (0, 0))

# ---------- matrix tie-ins: Papers 191/192/203 ----------
G_seed = sp.Matrix([[0, 1], [sqrt5, 0]])
G_adj  = sp.Matrix([[0, -1], [sqrt5, 0]])
I2 = sp.eye(2)
s3 = sp.Matrix([[1, 0], [0, -1]])
check("Gamma_seed^2 = +sqrt5 I", sp.simplify(G_seed**2 - sqrt5*I2) == sp.zeros(2))
check("Gamma_adj^2 = -sqrt5 I", sp.simplify(G_adj**2 + sqrt5*I2) == sp.zeros(2))
check("{Gamma_seed, Gamma_adj} = 0",
      sp.simplify(G_seed*G_adj + G_adj*G_seed) == sp.zeros(2))
check("Gamma_seed * Gamma_adj = +sqrt5 sigma_3  [Paper 192 measurement product]",
      sp.simplify(G_seed*G_adj - sqrt5*s3) == sp.zeros(2))
def sigma_matrix(M):
    return M.subs(sqrt5, -sqrt5)
check("sigma(Gamma_adj) = -Gamma_seed  [Paper 203 Face 1, our convention]",
      sp.simplify(sigma_matrix(G_adj) + G_seed) == sp.zeros(2))
check("sigma(Gamma_seed) = -Gamma_adj  [involution closes]",
      sp.simplify(sigma_matrix(G_seed) + G_adj) == sp.zeros(2))
ev_seed = list(G_seed.eigenvals().keys())
ev_adj = list(G_adj.eigenvals().keys())
check("spec(Gamma_seed) = {+-5^(1/4)} real",
      sorted([sp.simplify(e**2) for e in ev_seed], key=str) == [sqrt5, sqrt5])
check("spec(Gamma_adj) = {+-i 5^(1/4)} imaginary",
      sorted([sp.simplify(e**2) for e in ev_adj], key=str) == [-sqrt5, -sqrt5])

# ---------- the golden analogue of U(n): GL(n, R) ----------
# KK = R x R with swap involution; h-unitary M = (M+, M-) with sigma(M)^T M = I
M11, M12, M21, M22 = sp.symbols('M11 M12 M21 M22', real=True)
Mp = sp.Matrix([[M11, M12], [M21, M22]])
Mm = (Mp.T)**(-1)
check("pairs (M, (M^T)^-1) satisfy the sigma-unitarity condition componentwise "
      "[U(n, KK, sigma) ~ GL(n, R)]",
      sp.simplify(Mm.T*Mp - sp.eye(2)) == sp.zeros(2)
      and sp.simplify(Mp.T*Mm - sp.eye(2)) == sp.zeros(2))

print()
if FAILURES:
    print(f"{len(FAILURES)} FAILURES:")
    for f in FAILURES:
        print("  -", f)
    raise SystemExit(1)
print("ALL CHECKS PASSED (exact symbolic, sympy", sp.__version__ + ")")
