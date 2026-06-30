"""Minimal N=1000 FSS point: frustrated relaxation + d_frust only (fast path, teardown-
resistant). Control #3 already validated at N=100/200/500; convergence re-run owed.
Combines with the saved 3-point ladder and does the 4-point d_inf fit."""
import json, time
import numpy as np
from geometry import generic_s3
from relax_fast import precompute_fss, relax_fss
from spectral_s3 import low_spectrum, normalized_spectrum, spectrum_match_error
from fast_geometry import precompute as precompute_curv, fast_energy
from scipy.optimize import curve_fit

N = 1000
t0 = time.time()
cx = generic_s3(N, seed=20260618)
pre = precompute_fss(cx); pcv = precompute_curv(cx)
ell1 = np.ones(len(cx.edges))
P = cx.hull_points
ellr = np.array([np.linalg.norm(P[i]-P[j]) for (i,j) in cx.edges])

relaxed, info = relax_fss(cx, ell1, pre=pre, iters=300)
print(f"N=1000 relaxed: E {info['E0']:.0f}->{info['Ef']:.1f} iters={info['iters']} {info['msg'][:30]} ({time.time()-t0:.0f}s)", flush=True)

def sd(a,b,M=12):
    a=normalized_spectrum(a); b=normalized_spectrum(b); m=min(M,len(a),len(b))
    return float(np.sqrt(np.mean(((a[:m]-b[:m])/b[:m])**2)))
spec_i,_=low_spectrum(cx,ell1,k=40); spec_r,_=low_spectrum(cx,relaxed,k=40); spec_o,_=low_spectrum(cx,ellr,k=40)
_,cf,_,_=fast_energy(relaxed,pcv); _,cr,_,_=fast_energy(ellr,pcv)
rec=dict(N=N,edges=len(cx.edges),tets=len(cx.tets),E0=info['E0'],Ef=info['Ef'],
         relax_iters=info['iters'],relax_msg=info['msg'],
         d_frust_init=sd(spec_i,spec_o),d_frust=sd(spec_r,spec_o),
         d_round=spectrum_match_error(spec_o,4),cov_frust=float(cf),cov_round=float(cr),
         cov_gap=float(cf-cr),control3_stays_round="owed (validated at N<=500)")
print(f"  d_frust_init={rec['d_frust_init']:.3f} -> d_frust={rec['d_frust']:.3f}  d_round={rec['d_round']:.3f}", flush=True)

rows = json.load(open("results_fss_3pt.json"))
rows = [r for r in rows if r["N"] != 1000] + [rec]
rows.sort(key=lambda r: r["N"])
Ns=np.array([r["N"] for r in rows],float); ds=np.array([r["d_frust"] for r in rows])
m=lambda N,dinf,c,p: dinf+c*N**(-p)
fit=None
try:
    popt,pcov=curve_fit(m,Ns,ds,p0=[0.02,1,0.5],bounds=([-0.3,-50,0.05],[1,50,4]),maxfev=40000)
    perr=np.sqrt(np.diag(pcov))
    fit=dict(d_inf=float(popt[0]),d_inf_err=float(perr[0]),c=float(popt[1]),p=float(popt[2]),
             rms_resid=float(np.sqrt(np.mean((ds-m(Ns,*popt))**2))))
except Exception as e:
    fit=dict(error=str(e))
json.dump(dict(ladder=rows,fit=fit),open("results_fss.json","w"),indent=2)
print("\n=== 4-POINT LADDER ===")
for r in rows: print(f"  N={r['N']:5d}  d_frust={r['d_frust']:.4f}  d_round={r['d_round']:.3f}")
print("fit:",fit)
if fit and "d_inf" in fit:
    verdict = "ARTEFACT (d_inf~0)" if fit["d_inf"]-2*fit["d_inf_err"]<=0 else "FUNDAMENTAL glass (d_inf>0)"
    print(f"VERDICT: {verdict}  [d_inf={fit['d_inf']:.3f}+/-{fit['d_inf_err']:.3f}]")
