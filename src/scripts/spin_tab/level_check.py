#!/usr/bin/env python3
"""Screen the level energies of the campaign (levels_out) for eigenvalue searches that converged to a wrong root.
    python3 level_check.py [threshold = 0.3] [out list]
n_eff = alpha / sqrt(2 (1 - Re w / alpha)) must be close to n (fine structure shifts it by a few alpha^2 / (l + 1/2));
an entry with |n_eff - n| > threshold, or with Re w >= alpha, is listed.  Writes 'n l m spin alpha' of the bad entries."""
import glob, math, sys, collections
X = "/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext/spin_tab/"
thr = float(sys.argv[1]) if len(sys.argv) > 1 else 0.3
E = {}
for f in glob.glob(X + "levels_out/lev_*.txt"):
    for ln in open(f):
        p = ln.split()
        if len(p) == 9: E[(int(p[0]), int(p[1]), int(p[2]), float(p[3]), float(p[4]))] = (float(p[5]), float(p[6]), int(p[7]))
bad = []; dev = []
for (n, l, m, a, al), (w, wi, ok) in E.items():
    x = 1 - w / al
    ne = al / math.sqrt(2 * x) if x > 0 else float("inf")
    d = ne - n
    # leading fine-structure shift of n_eff (Baumann et al. 2019): delta(w/alpha) = -alpha^4/(8 n^4) + alpha^4 (2l - 3n + 1)/(n^4 (l + 1/2))
    fs = -(n ** 3 / al ** 2) * (-al ** 4 / (8 * n ** 4) + al ** 4 * (2 * l - 3 * n + 1) / (n ** 4 * (l + 0.5)))
    dev.append((abs(d - fs), n, l, m, a, al))
    if abs(d) > thr or ok != 1: bad.append((n, l, m, a, al, ne, w, ok))
dev.sort()
print("%d level energies; |n_eff - n - fine structure|: median %.4f 99%% %.4f; entries with |n_eff - n| > %g: %d" % (len(E), dev[len(dev) // 2][0], dev[int(.99 * len(dev))][0], thr, len(bad)))
print("by spin:", dict(sorted(collections.Counter(b[3] for b in bad).items())))
print("by alpha:", dict(sorted(collections.Counter(b[4] for b in bad).items())))
lv = collections.Counter((b[0], b[1], b[2]) for b in bad)
print("levels affected: %d;" % len(lv), ", ".join("%d-%d-%d: %d" % (k + (v,)) for k, v in sorted(lv.items())))
print("examples (n l m spin alpha n_eff):")
for b in sorted(bad)[:12]: print("  %d %d %d  a=%g alpha=%g  n_eff=%.2f" % b[:6])
gray = [d for d in dev if 0.05 < d[0] and not (abs(d[0]) > thr)]
print("entries with 0.05 < |n_eff - n - fine structure| (not flagged): %d" % sum(1 for d in dev if d[0] > 0.05) , "; largest unflagged:", [("%d-%d-%d a=%g al=%g: %.3f" % (d[1], d[2], d[3], d[4], d[5], d[0])) for d in dev[-len(bad) - 4:-len(bad)]] if bad else "")
if len(sys.argv) > 2:
    with open(sys.argv[2], "w") as o:
        o.write("# level energies of levels_out that sit on a wrong root (level_check.py, |n_eff - n| > %g): n l m spin alpha\n" % thr)
        for b in sorted(bad): o.write("%d %d %d %g %g\n" % b[:5])
