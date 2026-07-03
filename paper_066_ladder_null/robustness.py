#!/usr/bin/env python3
"""paper66-ladder-null: cold-pass robustness check (Mr A, 10 Jun 2026).
Settings declared before running: T1 = (max 2 terms, |k|<=3); T2 = (max 3 terms,
|k|<=3); reference T0 = original (max 3 terms, |k|<=7). Same monomial pool as
enumerate.py. Reports coverage of [30,200] and ladder-integer spellings, to show the
DECORATIVE verdict is not an artefact of a loose class boundary."""
from itertools import combinations, product
inv = {"V":12,"E":30,"F":20,"D":3,"chi":2}
names=list(inv); mono={"1":1}; mono.update(inv)
for i,a in enumerate(names):
    for b in names[i:]: mono[f"{a}*{b}"]=inv[a]*inv[b]
mono_items=sorted(mono.items(), key=lambda kv: kv[1]); LO,HI=30,200
ladder=[44,84,137,140,142,168,184]
def run(max_terms,kmax,tag):
    C=[k for k in range(-kmax,kmax+1) if k]
    cnt={t:0 for t in range(LO,HI+1)}
    def rec(v):
        if LO<=v<=HI: cnt[v]+=1
    for (mn,mv) in mono_items:
        for k in C: rec(k*mv)
    if max_terms>=2:
        for (an,av),(bn,bv) in combinations(mono_items,2):
            for ka,kb in product(C,C): rec(ka*av+kb*bv)
    if max_terms>=3:
        for (an,av),(bn,bv),(cn,cv) in combinations(mono_items,3):
            for ka,kb,kc in product(C,C,C): rec(ka*av+kb*bv+kc*cv)
    cov=sum(1 for t in cnt if cnt[t]); n=HI-LO+1
    print(f"{tag}: coverage {cov}/{n} = {100*cov/n:.1f}%   "
          f"mean spellings {sum(cnt.values())/n:,.0f}")
    print("   ladder: " + ", ".join(f"{t}:{cnt[t]:,}" for t in ladder))
    miss=[t for t in cnt if not cnt[t]]
    if miss: print(f"   unexpressible: {miss}")
run(2,3,"T1 (<=2 terms, |k|<=3)")
run(3,3,"T2 (<=3 terms, |k|<=3)")
run(3,7,"T0 (original: <=3 terms, |k|<=7)")
