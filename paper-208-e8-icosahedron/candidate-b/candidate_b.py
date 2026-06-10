"""
Candidate B — E8 = 2I + star(2I), the discarded golden content.
Construction per Dechant 2016, Sec 6. Exact arithmetic over Z[tau].

Plan:
  - Work in the even subalgebra of Cl(3) = quaternions H, spanned by {1, e23, e31, e12}.
    The 120 even spinors = 2I = binary icosahedral group = H4 root system (600-cell vertices).
  - The 120 odd pinors = I * (even set), I = e1e2e3 the pseudoscalar.
  - Pseudoscalar multiplication = Hodge duality (the star).
  - tau = (1+sqrt5)/2.

We generate 2I directly as the 120 unit icosians (standard, matches Dechant's H4 roots),
in the quaternion basis (a + b i + c j + d k) with a,b,c,d in Z[tau]/2.
Then verify: 120 elements, closed under multiplication (it's a group), all unit norm.
"""

import sympy as sp
from itertools import product
from collections import defaultdict

tau = sp.Rational(1,2) + sp.sqrt(5)/2     # golden ratio, exact
half = sp.Rational(1,2)

# A quaternion is (a,b,c,d) with entries in Z[tau]. We'll keep them as sympy exprs.
def qmul(x, y):
    a1,b1,c1,d1 = x; a2,b2,c2,d2 = y
    return (
        sp.expand(a1*a2 - b1*b2 - c1*c2 - d1*d2),
        sp.expand(a1*b2 + b1*a2 + c1*d2 - d1*c2),
        sp.expand(a1*c2 - b1*d2 + c1*a2 + d1*b2),
        sp.expand(a1*d2 + b1*c2 - c1*b2 + d1*a2),
    )

def qnorm2(x):
    return sp.expand(sum(t**2 for t in x))

def key(x):
    # canonical hashable key after simplifying in Z[tau]
    return tuple(sp.nsimplify(sp.simplify(t)) for t in x)

# --- The 120 icosians (binary icosahedral group) as unit quaternions ---
# Standard construction: 
#   8 from (+-1,0,0,0) and perms  -> the 8 "quaternion units" *2... actually:
#   24 Hurwitz-style: (+-1,0,0,0)&perms (8) and (+-1,+-1,+-1,+-1)/2 (16)
#   96 even permutations of (0, +-1/2, +-tau/2, +-(tau-1)/2 ) i.e. (0,+-1,+-tau,+-1/tau)/2 even perms
# This is the standard 600-cell vertex set = 2I.

elems = set()
store = []

def add(q):
    k = key(q)
    if k not in elems:
        elems.add(k)
        store.append(q)

# Group 1: (+-1,0,0,0) and all coordinate permutations -> 8
for i in range(4):
    for s in (1,-1):
        v = [0,0,0,0]; v[i] = s
        add(tuple(v))

# Group 2: (+-1/2,+-1/2,+-1/2,+-1/2) -> 16
for signs in product((1,-1), repeat=4):
    add(tuple(half*s for s in signs))

# Group 3: even permutations of (0, +-1/2, +-tau/2, +-(1/tau)/2)
#   with 1/tau = tau - 1.  Pattern values: 0, +-1/2, +-tau/2, +-(tau-1)/2
inv = tau - 1   # 1/tau
base_vals = [sp.Integer(0), half, tau*half, inv*half]
# even permutations of positions (the alternating group A4 on 4 slots)
A4 = [p for p in __import__('itertools').permutations(range(4))
      if sum(1 for i in range(4) for j in range(i+1,4) if p[i]>p[j]) % 2 == 0]

for perm in A4:
    # assign base magnitudes to slots per this even permutation, with sign choices
    # slot getting 0 has no sign; the other three each +-
    mags = [base_vals[perm.index(i)] for i in range(4)]  # place base_vals by perm
    # find which slot is zero
    nz = [i for i in range(4) if mags[i] != 0]
    for signs in product((1,-1), repeat=3):
        v = list(mags)
        for idx, s in zip(nz, signs):
            v[idx] = v[idx]*s
        add(tuple(v))

print("count after construction:", len(store))
print("distinct norms:", set(sp.simplify(qnorm2(q)) for q in store))

# --- Verify closure: it must be a group (this is what makes it 2I, not a random set) ---
keyset = set(elems)
closed = True
checks = 0
import random
random.seed(0)
sample = store  # full check is 120*120=14400, fine
for x in sample:
    for y in sample:
        if key(qmul(x,y)) not in keyset:
            closed = False
            print("NOT CLOSED at", x, y)
            break
        checks += 1
    if not closed: break
print(f"closure verified over {checks} products:", closed)
