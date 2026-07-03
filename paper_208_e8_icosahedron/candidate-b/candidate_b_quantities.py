import sympy as sp
from itertools import product, combinations, permutations
from collections import Counter

tau = sp.Rational(1,2) + sp.sqrt(5)/2
half = sp.Rational(1,2)
sqrt5 = sp.sqrt(5)

# ---- rebuild the 120 icosians (2I = H4 roots) ----
def add(store, elems, q):
    k = tuple(sp.nsimplify(sp.simplify(t)) for t in q)
    if k not in elems:
        elems.add(k); store.append(tuple(sp.expand(t) for t in q))

store=[]; elems=set()
for i in range(4):
    for s in (1,-1):
        v=[0,0,0,0]; v[i]=s; add(store,elems,tuple(v))
for signs in product((1,-1),repeat=4):
    add(store,elems,tuple(half*s for s in signs))
inv=tau-1
base=[sp.Integer(0),half,tau*half,inv*half]
A4=[p for p in permutations(range(4)) if sum(1 for i in range(4) for j in range(i+1,4) if p[i]>p[j])%2==0]
for perm in A4:
    mags=[base[perm.index(i)] for i in range(4)]
    nz=[i for i in range(4) if mags[i]!=0]
    for signs in product((1,-1),repeat=3):
        v=list(mags)
        for idx,s in zip(nz,signs): v[idx]=v[idx]*s
        add(store,elems,tuple(v))

H4 = store  # 120 roots, the even set = 2I
assert len(H4)==120

# pseudoscalar / Hodge-dual copy: I*(even). In the quaternion (even-subalgebra) model
# the Hodge dual of the even set is realised as the SECOND H4 copy. Dechant: E8 roots
# = H4  U  scaling*(I.H4). The two conventions differ in whether that second copy is
# tau-scaled (A, Dechant's stated op) or bare (B).
#
# The "odd"/dual copy as a point set is an isometric image of H4. The inner-product
# structure that matters for Q1/Q1x is: within-copy products, and CROSS products
# between a root in copy 1 and a (scaled) root in copy 2.
#
# Represent every root as a 4-vector. Euclidean inner product. Entries live in Z[tau].

def dot(x,y):
    return sp.expand(sum(a*b for a,b in zip(x,y)))

# Decompose an exact value v in Q(sqrt5) into a + tau*b with a,b rational.
# v = p + q*sqrt5  ->  tau=(1+sqrt5)/2 => sqrt5 = 2*tau-1 => v = p + q*(2tau-1) = (p-q) + (2q)*tau
def to_a_b(v):
    v = sp.expand(v)
    p = v.subs(sqrt5,0)            # rational part
    # coefficient of sqrt5:
    q = sp.expand((v - p)/sqrt5)
    p = sp.nsimplify(p); q = sp.nsimplify(q)
    a = sp.nsimplify(p - q)
    b = sp.nsimplify(2*q)
    return a, b

def build_roots(scale_second):
    copy1 = [(r,0) for r in H4]                       # tag sector 0
    copy2 = [(tuple(scale_second*c for c in r),1) for r in H4]  # tag sector 1
    return copy1+copy2

def quantities(scale_second, label):
    roots = build_roots(scale_second)
    n=len(roots)
    Q1=sp.Integer(0); Q1x=sp.Integer(0); within=sp.Integer(0)
    bhist=Counter()
    for (x,sx),(y,sy) in combinations(roots,2):
        a,b = to_a_b(dot(x,y))
        ab = abs(b)
        Q1 += ab
        bhist[b]+=1
        if sx!=sy: Q1x += ab
        else: within += ab
    Q2 = Q1x/within if within!=0 else sp.oo
    print(f"--- Convention {label} (second copy scale = {scale_second}) ---")
    print(f"  pairs: {n*(n-1)//2}")
    print(f"  Q1  (total |b|)            = {sp.nsimplify(Q1)}  ~ {float(Q1):.6f}")
    print(f"  Q1x (cross-sector |b|)     = {sp.nsimplify(Q1x)}  ~ {float(Q1x):.6f}")
    print(f"  within (= Q1 - Q1x)        = {sp.nsimplify(within)}  ~ {float(within):.6f}")
    print(f"  Q2  = Q1x/within           = {sp.nsimplify(Q2)}  ~ {float(Q2):.6f}")
    # normalisations
    C=n*(n-1)//2
    for name,div in [("per-pair (/C(240,2))",C),("per-root (/240)",n),("raw",1)]:
        print(f"    Q1   {name:24s}: {float(Q1/div):.6f}")
        print(f"    Q1x  {name:24s}: {float(Q1x/div):.6f}")
    print(f"  Q3 histogram of b (top 12 by |b|):")
    for bval,cnt in sorted(bhist.items(), key=lambda kv:-abs(float(kv[0])))[:12]:
        print(f"      b={sp.nsimplify(bval)} ~{float(bval):+.4f}  count={cnt}")
    return dict(Q1=Q1,Q1x=Q1x,within=within,Q2=Q2,bhist=bhist)

print("="*70)
A = quantities(tau, "A")    # Dechant stated: second copy tau-scaled
print("="*70)
B = quantities(sp.Integer(1), "B")  # bare pinor split, unscaled
print("="*70)

# Targets
mu = 6*sp.pi**5
ainv = sp.pi + sp.pi**2 + 4*sp.pi**3
print(f"TARGETS:  mu=6pi^5 = {float(mu):.6f}   alpha^-1 = pi+pi^2+4pi^3 = {float(ainv):.6f}")
