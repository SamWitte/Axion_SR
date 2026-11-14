#!/usr/bin/env python3
"""Summary of the extended-SG evolution runs: final spin / mass, levels above a threshold, top of each tower.
    python3 summarize_runs.py [--set full|core] [--thr 1e-20] [--csv out.csv]
Reads runs/<set>/gw_<model>/M_<mass>/<variant>/alpha_<a>/{Time,Spin,MassBH,Modes,States}_*.dat (States: one row per level).
"""
import glob, os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
def opt(name, default):
    return sys.argv[sys.argv.index(name) + 1] if name in sys.argv else default
sset = opt("--set", "full"); thr = float(opt("--thr", "1e-20")); csv = opt("--csv", None)
rows = []
for d in sorted(glob.glob(os.path.join(HERE, "runs", sset, "gw_*", "M_*", "*", "alpha_*"))):
    gw, mass, var, al = d.split(os.sep)[-4:]
    ft = glob.glob(os.path.join(d, "Time_*.dat"))
    if not ft: rows.append((gw, mass, var, al, "no output", "", "", "", "")); continue
    t = np.loadtxt(ft[0]); spin = np.loadtxt(glob.glob(os.path.join(d, "Spin_*.dat"))[0]); mbh = np.loadtxt(glob.glob(os.path.join(d, "MassBH_*.dat"))[0])
    modes = np.loadtxt(glob.glob(os.path.join(d, "Modes_*.dat"))[0])[:, :3].astype(int)
    pk = []
    with open(glob.glob(os.path.join(d, "States_*.dat"))[0]) as f:
        for i, line in enumerate(f):
            r = np.fromstring(line, sep="\t"); pk.append((r.max(), r[-1]))
    pk = np.array(pk); sel = np.where(pk[:, 0] > thr)[0]; sel = sel[np.argsort(-pk[sel, 0])]
    lab = lambda m: f"{m[0]}{m[1]}{m[2]}" if max(m) < 10 else f"{m[0]}-{m[1]}-{m[2]}"
    tower = {}
    for i in sel:
        n, l, m = modes[i]
        if l == m and n > tower.get((l, m), 0): tower[(l, m)] = n
    done = "done" if t[-1] >= 0.99e9 else f"t={t[-1]:.2e}"
    rows.append((gw, mass, var, al, done, f"{spin[-1]:.4f}", f"{mbh[-1] / mbh[0]:.4f}", str(len(sel)),
                 " ".join(f"{lab(modes[i])}:{pk[i, 0]:.1e}" for i in sel[:12]) + "  | tower tops (l=m): " +
                 " ".join(f"({l},{m})->n={n}" for (l, m), n in sorted(tower.items()))))
hdr = ("gw", "mass", "variant", "alpha", "status", "a_final", "M_final/M0", f"n_levels>{thr:g}", "peak occupations (largest first)")
for r in [hdr] + rows: print(" | ".join(r))
if csv:
    with open(csv, "w") as f:
        for r in [hdr] + rows: f.write(",".join('"' + x + '"' for x in r) + "\n")
