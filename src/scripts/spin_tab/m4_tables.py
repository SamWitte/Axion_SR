#!/usr/bin/env python3
"""Corrected rate tables for the 57 horizon channels with m4 != 0 (a = 0.95), in the normalisation of the stored tables
(spin_fast.jl ... kg=0 hmax=0.25 tag=kg0, spin 0.95).   python3 m4_tables.py <dir with kg0_*.txt>
Writes (nothing in the production tables is touched):
  m4fix/tables/<ch>_LvrHc_SG.dat   the rate_sve table with its 14 grid rows replaced (rows above alpha 0.28 kept)
  m4fix/pts/<ch>.pts               the same points as rows for rates_out (Compute_sg_points format; h_fac 1000 so that
                                   merge_rates.py prefers them) - to be appended only when the tables are to be replaced
  overlay_n88_m4fix/               NL_SG_OVERLAY directory: links to every file of ../n88_overlay plus the 57 tables"""
import glob, os, sys, math
X = "/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext"; T = X + "/spin_tab"; RS = "/data/dust/user/samwitte/Axion_SR_sg/src/rate_sve/"
new = {}
for f in glob.glob(sys.argv[1] + "/kg0_*.txt"):
    spins = None
    for ln in open(f):
        p = ln.split()
        if ln.startswith("# spins"):
            spins = [float(x) for x in p[2:p.index("hmax")]]; assert "kg false" in ln; continue
        j = spins.index(0.95); q = p[6 + 4 * j: 10 + 4 * j]
        if q[0] == "1": new.setdefault(tuple(p[:4]), {})[float(p[4])] = tuple(float(x) for x in q[1:])
os.makedirs(T + "/m4fix/tables", exist_ok=True); os.makedirs(T + "/m4fix/pts", exist_ok=True)
nrow = 0
for c, d in sorted(new.items()):
    name = "_".join(c) + "_LvrHc_SG.dat"
    old = [l.split() for l in open(RS + name) if l.strip()]
    out = []; 
    for r in old:
        al = float(r[0]); k = next((a for a in d if abs(a - al) < 1e-9), None)
        if k is None: out.append("\t".join(r)); continue
        si, sg, x = d[k]; kap = x / (2 * math.sqrt(si * sg)) if si > 0 and sg > 0 else 0.0
        out.append("\t".join(repr(v) for v in (k, si, sg, x, kap))); nrow += 1
    open(T + "/m4fix/tables/" + name, "w").write("\n".join(out) + "\n")
    st = [tuple(int(v) for v in s.split("-")) if "-" in s else tuple(int(ch) for ch in s) for s in c[:3]]
    mask = sum((1 << i) for i, s in enumerate(st) if s[0] > 8)
    with open(T + "/m4fix/pts/" + "_".join(c) + ".pts", "w") as o:
        for k in sorted(d):
            si, sg, x = d[k]; kap = x / (2 * math.sqrt(si * sg)) if si > 0 and sg > 0 else 0.0
            o.write("\t".join(str(v) for v in (k, si, sg, x, kap, 0.0, mask, 1000.0, 0)) + "\n")
ov = T + "/overlay_n88_m4fix"; os.makedirs(ov, exist_ok=True)
n = 0
for f in os.listdir(X + "/n88_overlay"):
    dst = ov + "/" + f
    if not os.path.lexists(dst): os.symlink(X + "/n88_overlay/" + f, dst); n += 1
for c in new:
    name = "_".join(c) + "_LvrHc_SG.dat"; dst = ov + "/" + name
    if os.path.lexists(dst): os.remove(dst)
    os.symlink(T + "/m4fix/tables/" + name, dst)
print("%d channels, %d grid rows replaced; overlay %s: %d links to n88_overlay + %d corrected tables" % (len(new), nrow, ov, n, len(new)))
