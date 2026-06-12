#!/usr/bin/env python3
"""paper66-ladder-null addendum 3: the Lambda_QCD rival (briefed-pass steelman, J2).
PRE-REGISTERED candidates: six standard scheme/value choices for Lambda_QCD spanning
PDG MS-bar nf=3..5 and legacy values: 210, 250, 263, 290, 332, 340 MeV.
Question: does the rival physical length hbar*c/Lambda, in units u = lambdabar_p/4,
land on integers stably (signal) or sweep them (definitional mush)?"""
hbarc=197.3269804; mp=938.27208816; u=hbarc/mp/4
print(f"u = {u:.6f} fm")
print("Lambda (MeV)   length (fm)   N = len/u   nearest int (dist)")
for L in (210,250,263,290,332,340):
    l=hbarc/L; n=l/u; d=abs(n-round(n))
    print(f"  {L:>5}        {l:.4f}       {n:7.3f}    {round(n)} ({d:.3f})")
print("\nVerdict: N sweeps ~11-19 across schemes, including a spuriously sharp")
print("15.012 at Lambda=250 MeV -- definitional promiscuity demonstrated on the")
print("nearest physical rival. Strengthens the Sec 3.2 no-claim demotion (rule R4).")
