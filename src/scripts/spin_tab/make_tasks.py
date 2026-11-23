#!/usr/bin/env python3
"""Task lists of the spin-table campaign (nothing is submitted here).   python3 make_tasks.py
  levels.txt     the 215 levels
  fast_bh.txt    every channel whose fourth quantum can be bound ("S1 S2 S3 S4"): input of spin_fast.jl.  Left out: the
                 horizon channels with |m4| >= 20 (l4 = 24 or 30 at the horizon: the rate underflows, zero in the tables)
  fast_n8.txt    the subset with a level n = 8: recomputed above a = 0.95 with that level hydrogenic (the relativistic
                 n = 8 radial function is not reliable at a = 0.998)
  fast_m4.txt    the subset with m4 != 0 (57 channels; for local tests - they are part of fast_bh)
  full_check.txt full-path points (spin_points.jl: "S1 S2 S3 S4 spin alpha hfac BHlmax") to validate spin_fast.jl under
                 the campaign switches: 270 random m4 = 0 horizon channels (3 spins + 0.95, 2-3 alpha) and 6 m4 != 0 ones"""
import random, collections, os
X = "/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext"
T = X + "/spin_tab/tasks"
st = lambda t: tuple(int(x) for x in t.split("-")) if "-" in t else tuple(int(c) for c in t)
chans = [tuple(l.split()[:4]) for l in open(X + "/sgx_channels.txt") if l.strip()]
sec = {}
for l in open(X + "/spin_dep/plan_scratch/census-cost/pts_all_perchan.tsv"):
    p = l.rstrip("\n").split("\t"); sec[p[0]] = {x.split(":")[0]: float(x.split(":")[1]) for x in p[1:]}
active = {tuple(l.split()[1:5]) for l in open(X + "/spin_dep/plan_scratch/runs-importance/active_all_1e-20.txt")}
def info(c):
    s = [st(x) for x in c[:3]]
    return dict(m4=s[0][2] + s[1][2] - s[2][2], nmax=max(x[0] for x in s), nmin=min(x[0] for x in s), ns=[x[0] for x in s],
                bind=1 / s[0][0] ** 2 + 1 / s[1][0] ** 2 - 1 / s[2][0] ** 2)
I = {c: info(c) for c in chans}
os.makedirs(T, exist_ok=True)
levels = sorted({st(x) for c in chans for x in c[:3]})
open(T + "/levels.txt", "w").write("\n".join("%d-%d-%d" % l for l in levels) + "\n")
fast = [c for c in chans if (c[3] == "BH" and abs(I[c]["m4"]) < 20) or (c[3] != "BH" and I[c]["bind"] >= -1e-9)]
fast.sort(key=lambda c: (c not in active, I[c]["nmax"]))          # populated channels and low n first
open(T + "/fast_bh.txt", "w").write("\n".join(" ".join(c) for c in fast) + "\n")
n8 = [c for c in fast if 8 in I[c]["ns"]]
open(T + "/fast_n8.txt", "w").write("\n".join(" ".join(c) for c in n8) + "\n")
m4 = [c for c in fast if c[3] == "BH" and I[c]["m4"] != 0]
open(T + "/fast_m4.txt", "w").write("\n".join(" ".join(c) for c in sorted(m4, key=lambda c: (I[c]["nmax"], abs(I[c]["m4"])))) + "\n")
print("%d levels; fast_bh %d channels (BH label %d, m4 != 0 %d, Inf label with a possibly bound fourth quantum %d); fast_n8 %d" %
      (len(levels), len(fast), sum(c[3] == "BH" for c in fast), len(m4), sum(c[3] != "BH" for c in fast), len(n8)))
print("left out (|m4| >= 20):", sum(1 for c in chans if c[3] == "BH" and abs(I[c]["m4"]) >= 20))
# full-path validation sample
random.seed(20261120)
bh0 = [c for c in chans if c[3] == "BH" and I[c]["m4"] == 0]
cls = {"i": [c for c in bh0 if I[c]["nmax"] <= 8],
       "ii": [c for c in bh0 if 8 < I[c]["nmax"] <= 20 and I[c]["nmin"] <= 8],
       "iii": [c for c in bh0 if I[c]["nmin"] > 8 and I[c]["nmax"] <= 20]}
L = []; seen = set(); cost = 0.0
def add(c, al, spins, hfac="1"):
    global cost
    for a in list(spins) + ["0.95"]:
        k = (c, a, al)
        if k in seen: continue
        seen.add(k); L.append(" ".join(c) + " %s %s %s 2" % (a, al, hfac)); cost += sec.get("_".join(c), {}).get(al, 60.0)
for name, n, als in (("i", 60, ["0.05", "0.1", "0.2"]), ("ii", 150, ["0.05", "0.1", "0.2"]), ("iii", 60, ["0.05", "0.2"])):
    for c in random.sample(cls[name], min(n, len(cls[name]))):
        for al in als: add(c, al, ["0", "0.5", "0.9"])
for c in (("422", "422", "433", "BH"), ("422", "422", "411", "BH"), ("411", "433", "422", "BH"), ("644", "644", "655", "BH"), ("633", "633", "644", "BH"), ("644", "644", "633", "BH")):
    assert c in I
    for al in ["0.1", "0.2"]: add(c, al, ["0", "0.5", "0.9"])
random.shuffle(L)                                                # mix cheap and expensive points over the workers
open(T + "/full_check.txt", "w").write("\n".join(L) + "\n")
print("full_check %d points, ~%.0f core-h at the 0.95 timings (the m4 != 0 points with the fine step cost more)" % (len(L), cost / 3600))
