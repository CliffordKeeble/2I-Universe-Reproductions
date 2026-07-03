#!/usr/bin/env python3
"""paper66-ladder-null addendum: the 19 through the null (Mr A, 2nd cold pass, w2).
PRE-REGISTERED before running:
Q1 (definitional anchor test): the count N(r) = r/u, u = lambdabar_p/4, is computed
for every PHYSICALLY DEFINED candidate for "the strong-force/confinement range":
proton rms charge radius 0.8409 fm; charged-pion Yukawa range hbar/(m_pi+ c) ;
neutral-pion Yukawa range; and separately for the DECIMAL CONVENTION 1.000 fm.
DECISION RULE R4: a count is bankable only if anchored to a physically defined
length; a match produced only by the round decimal figure is a fact about the
metric system and is carried at NO CLAIM (definition-dominated, the cosmic-row rule).
Q2 (in-class density): spellings of every integer in [12, 29] under class C of the
main null (<=3 terms, |k|<=7, 21 monomials), to answer "how special is 19 as an
integer" — same promiscuity standard as the ladder.
"""
from itertools import combinations, product
hbarc=197.3269804; mp=938.27208816
u=hbarc/mp/4
print(f"u = {u:.6f} fm")
cands=[("proton rms charge radius",0.8409),("charged-pion Yukawa hbar/m_pi c",hbarc/139.57039),
       ("neutral-pion Yukawa",hbarc/134.9768),("DECIMAL CONVENTION 1.000 fm",1.0)]
print("\nQ1 — counts for physically defined ranges vs the decimal convention:")
for tag,r in cands:
    n=r/u; d=abs(n-round(n))
    print(f"  {tag:34s} r={r:.4f} fm  N={n:6.2f}  nearest int {round(n)} (dist {d:.2f})")
inv={"V":12,"E":30,"F":20,"D":3,"chi":2}
names=list(inv); mono={"1":1}; mono.update(inv)
for i,a in enumerate(names):
    for b in names[i:]: mono[f"{a}*{b}"]=inv[a]*inv[b]
mi=sorted(mono.items(),key=lambda kv:kv[1]); C=[k for k in range(-7,8) if k]
cnt={t:0 for t in range(12,30)}
def rec(v):
    if 12<=v<=29: cnt[v]+=1
for (mn,mv) in mi:
    for k in C: rec(k*mv)
for (an,av),(bn,bv) in combinations(mi,2):
    for ka,kb in product(C,C): rec(ka*av+kb*bv)
for (an,av),(bn,bv),(cn,cv) in combinations(mi,3):
    for ka,kb,kc in product(C,C,C): rec(ka*av+kb*bv+kc*cv)
print("\nQ2 — class-C spellings, integers 12..29:")
print("  " + "  ".join(f"{t}:{cnt[t]:,}" for t in sorted(cnt)))
print(f"\n  19 carries {cnt[19]:,} in-class spellings — same promiscuity as the ladder.")
