"""
run_arm1.py -- Paper 213 Arm 1 driver (brief sec.6).

Runs, in order:
  (A) finite-size scaling of d_s on the default rule (golden seed, GR-1 rho,
      C=flat-triangles, d_max=4), with the <12,20,30> signature on BOTH L_s and
      L_c at every N;
  (B) controls (non-golden seeds, scrambled coupling) at a fixed N;
  (C) rule-space sweep (rho x C x d_max) logged to logs/rule_trials.jsonl;
  (D) figures + machine-readable results.json.

Usage:
  python run_arm1.py --scale debug     # N up to 8k        (minutes)
  python run_arm1.py --scale prod      # N up to 128k      (longer)

The pre-registered verdict criteria live in RESULTS.md, committed BEFORE this
runs (Pattern 75). This script only PRODUCES numbers; the verdict is read off
the signature scores vs the control band.
"""

import argparse
import json
import time
import numpy as np

from growth import grow, rho_GR1, rho_allR
from laplacian import scalar_laplacian, connection_laplacian
from spectral import low_spectrum, spectral_dimension, signature_score
from controls import CONTROL_SPECS, DMAX_SWEEP, C_VARIANTS

HERE = __file__.rsplit("/", 1)[0] if "/" in __file__ else "."
import os
HERE = os.path.dirname(os.path.abspath(__file__))
LOGDIR = os.path.join(HERE, "logs")
FIGDIR = os.path.join(HERE, "figures")


def kfor(N):
    return int(min(800, max(50, N // 5)))


def measure(G, meta, em, k=None):
    """Compute d_s (from L_s) and signature on both Laplacians for one graph."""
    N = meta["N"]
    k = k or kfor(N)
    Ls, _ = scalar_laplacian(G)
    Lc, _ = connection_laplacian(G, em)
    es = low_spectrum(Ls, k=k)
    ec = low_spectrum(Lc, k=k)
    ds, dse, win = spectral_dimension(es, N)
    return {
        "N": N, "E": meta["E"], "d_max": meta["d_max"], "C": meta["C"],
        "rho": meta["rho"], "seed_a": meta["seed_a"],
        "scramble": meta["scramble"], "mean_degree": meta["mean_degree"],
        "n_triangles": meta["n_triangles"],
        "d_s": ds, "d_s_se": dse, "t_window": list(win),
        "sig_Ls": signature_score(es),
        "sig_Lc": signature_score(ec),
        "lam_Ls_min": float(es.min()), "lam_Lc_min": float(ec.min()),
    }


def checkpoint(results):
    """Write results.json after every stage/N so a session-boundary interruption
    (which has bitten this work before) preserves partial FSS, not nothing."""
    with open(os.path.join(HERE, "results.json"), "w") as f:
        json.dump(results, f, indent=2, default=float)


def run_fss(N_list, results):
    rows = results.setdefault("fss", [])
    done = {r["N"] for r in rows}
    for N in N_list:
        if N in done:        # resume: skip N already checkpointed
            continue
        t0 = time.time()
        G, meta, em = grow(N, d_max=4, k=2, C="flat-triangles", rho=rho_GR1,
                           seed_a=5.0)
        r = measure(G, meta, em)
        r["wall_s"] = round(time.time() - t0, 1)
        rows.append(r)
        checkpoint(results)
        print(f"[FSS] N={N:>7} d_s={r['d_s']:.3f}+/-{r['d_s_se']:.3f}  "
              f"kk2(Ls)={r['sig_Ls']['kk2_score']:.2f} "
              f"kk2(Lc)={r['sig_Lc']['kk2_score']:.2f}  "
              f"({r['wall_s']}s)", flush=True)
    return rows


def run_controls(N, results):
    rows = results.setdefault("controls", [])
    done = {r["control"] for r in rows}
    for spec in CONTROL_SPECS:
        if spec["name"] in done:
            continue
        kw = {k: v for k, v in spec.items() if k != "name"}
        G, meta, em = grow(N, d_max=4, k=2, rho=rho_GR1, **kw)
        r = measure(G, meta, em)
        r["control"] = spec["name"]
        rows.append(r)
        checkpoint(results)
        print(f"[CTRL] {spec['name']:<18} d_s={r['d_s']:.3f}  "
              f"kk2(Ls)={r['sig_Ls']['kk2_score']:.2f} "
              f"kk2(Lc)={r['sig_Lc']['kk2_score']:.2f} "
              f"sg_hit(Lc)={r['sig_Lc']['semigroup_hit']:.2f}", flush=True)
    return rows


def run_rulespace(N, logpath, results):
    """rho x C x d_max sweep -> jsonl. Honest about nulls (maps the rule space)."""
    rows = results.setdefault("rulespace", [])
    done = {(r["rho"], r["C"], r["d_max"]) for r in rows}
    rhos = {"rho_GR1": rho_GR1, "rho_allR": rho_allR}
    f = open(logpath, "a")
    for rho_name, rho in rhos.items():
        for C in C_VARIANTS:
            for dmax in DMAX_SWEEP:
                if (rho_name, C, dmax) in done:
                    continue
                G, meta, em = grow(N, d_max=dmax, k=2, C=C, rho=rho,
                                   seed_a=5.0)
                r = measure(G, meta, em)
                rec = {"rho": rho_name, "C": C, "d_max": dmax,
                       "N": r["N"], "E": r["E"],
                       "mean_degree": round(r["mean_degree"], 3),
                       "n_triangles": r["n_triangles"],
                       "d_s": round(r["d_s"], 4),
                       "d_s_se": round(r["d_s_se"], 4),
                       "kk2_Ls": round(r["sig_Ls"]["kk2_score"], 3),
                       "kk2_Lc": round(r["sig_Lc"]["kk2_score"], 3),
                       "sg_hit_Lc": round(r["sig_Lc"]["semigroup_hit"], 3),
                       "gap_avoid_Lc": round(r["sig_Lc"]["gap_avoidance"], 3),
                       "n_neg_Lc": r["sig_Lc"]["n_negative"],
                       "lambda1_Lc": round(r["sig_Lc"]["lambda1"], 4)}
                f.write(json.dumps(rec) + "\n"); f.flush()
                rows.append(rec)
                checkpoint(results)
                print(f"[RULE] {rho_name}/{C}/d{dmax} d_s={rec['d_s']:.2f} "
                      f"kk2_Lc={rec['kk2_Lc']:.2f}", flush=True)
    f.close()
    print(f"[RULE] {len(rows)} (rho,C,d_max) trials -> {logpath}")
    return rows


def make_figures(fss, rulerows):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"[FIG] matplotlib unavailable ({e}); skipping figures")
        return
    # d_s vs N
    Ns = [r["N"] for r in fss]
    ds = [r["d_s"] for r in fss]
    se = [r["d_s_se"] for r in fss]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.errorbar(Ns, ds, yerr=se, marker="o", capsize=3)
    ax.set_xscale("log"); ax.set_xlabel("N (nodes)"); ax.set_ylabel("d_s")
    ax.axhline(3.0, ls="--", color="grey", lw=1, label="d=3 (reference)")
    ax.set_title("Arm 1: spectral dimension vs N (GR-1, golden seed)")
    ax.legend(); fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "d_s_vs_N.png"), dpi=130)
    plt.close(fig)

    # signature scores across rule space (Lc kk2)
    fig, ax = plt.subplots(figsize=(6, 4))
    kk2 = [r["kk2_Lc"] for r in rulerows]
    sg = [r["sg_hit_Lc"] for r in rulerows]
    ax.scatter(kk2, sg, alpha=0.6)
    ax.set_xlabel("kk2 score (L_c)"); ax.set_ylabel("semigroup hit (L_c)")
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
    ax.set_title("Arm 1: signature scores across rho x C x d_max")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "signature_scatter.png"), dpi=130)
    plt.close(fig)
    print(f"[FIG] wrote d_s_vs_N.png, signature_scatter.png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scale", choices=["debug", "prod"], default="debug")
    args = ap.parse_args()

    if args.scale == "debug":
        N_list = [1000, 2000, 4000, 8000]
        N_ctrl = 8000
        N_rule = 4000
    else:
        N_list = [1000, 2000, 4000, 8000, 16000, 32000, 64000, 128000]
        N_ctrl = 32000
        N_rule = 16000

    print(f"=== Paper 213 Arm 1 :: scale={args.scale} ===", flush=True)
    t0 = time.time()

    # Resume from a prior checkpoint of the SAME scale if present (a long prod run
    # may span a session boundary; debug-scale results.json must not be reused).
    rpath = os.path.join(HERE, "results.json")
    results = {"scale": args.scale}
    if os.path.exists(rpath):
        try:
            prev = json.load(open(rpath))
            if prev.get("scale") == args.scale:
                results = prev
                print(f"[RESUME] continuing from checkpoint "
                      f"({len(prev.get('fss', []))} FSS, "
                      f"{len(prev.get('controls', []))} ctrl, "
                      f"{len(prev.get('rulespace', []))} rule)", flush=True)
        except Exception:
            pass

    rule_log = os.path.join(LOGDIR, "rule_trials.jsonl")
    if not results.get("rulespace") and os.path.exists(rule_log):
        open(rule_log, "w").close()   # fresh sweep if not resuming rulespace

    fss = run_fss(N_list, results)
    ctrl = run_controls(N_ctrl, results)
    rule = run_rulespace(N_rule, rule_log, results)
    make_figures(fss, rule)

    results["wall_s"] = round(time.time() - t0, 1)
    checkpoint(results)
    print(f"=== done in {results['wall_s']}s -> results.json ===", flush=True)


if __name__ == "__main__":
    main()
