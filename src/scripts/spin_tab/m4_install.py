#!/usr/bin/env python3
"""Install the corrected tables of the 57 horizon channels with m4 != 0 (23 Nov 2026).
    python3 m4_install.py            report only
    python3 m4_install.py --apply    write rate_sve (and allW/tables), old files copied to <dir>/backup_pre_m4fix_20261123
Source: m4fix/tables/<ch>_LvrHc_SG.dat (m4_tables.py: the 14 grid rows from spin_fast.jl with the radial step that resolves
the near-horizon wavenumber |w4 - m4 Omega_H|, kg=0, i.e. the flat normalisation int |R|^2 r^2 dr = 1 of the stored points;
rows above alpha = 0.28 as they were).  Installed tables are in the Klein-Gordon normalisation like every table since
22 Nov 2026: grid rows x 1 / (Q1 Q2 Q3) with the charges of kg_norm/kg_levels.txt, exactly as merge_allW.py does for
horizon channels (which reads these source tables instead of the old production ones, so a new merge keeps them).
The SI-only tables <ch>_LvrHc_.dat of these channels are first-campaign archive tables on another alpha grid and are not
touched (the runs take SI from the SG table, NL_SI_FROM_SG=1)."""
import os, sys, math, shutil
X = "/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext"; SRC = X + "/spin_tab/m4fix/tables"
TARGETS = ["/data/dust/user/samwitte/Axion_SR_sg/src/rate_sve", X + "/allW/tables"]
apply = "--apply" in sys.argv
Q = {}
for l in open(X + "/kg_norm/kg_levels.txt"):
    p = l.split(); Q[(p[0], round(float(p[1]), 10))] = float(p[2])
def table(f):
    c = f.replace("_LvrHc_SG.dat", "").split("_"); out = []
    for l in open(os.path.join(SRC, f)):
        r = [float(x) for x in l.split()[:5]]
        if len(r) < 5: continue
        try: fac = 1.0 / (Q[(c[0], round(r[0], 10))] * Q[(c[1], round(r[0], 10))] * Q[(c[2], round(r[0], 10))])
        except KeyError: out.append(r); continue                      # off the grid: unchanged, as in merge_allW.py
        out.append([r[0], r[1] * fac, r[2] * fac, r[3] * fac, r[4]])
    return out
files = sorted(f for f in os.listdir(SRC) if f.endswith("_LvrHc_SG.dat"))
for T in TARGETS:
    bak = os.path.join(T, "backup_pre_m4fix_20261123"); n = 0; big = 0; rows = 0
    for f in files:
        new = table(f); dst = os.path.join(T, f)
        if not os.path.isfile(dst): print("  no table to replace:", dst); continue
        cur = {round(float(l.split()[0]), 10): float(l.split()[2]) for l in open(dst) if len(l.split()) >= 5}
        for r in new:
            a = round(r[0], 10)
            if a in cur and cur[a] > 0 and r[2] > 0:
                rows += 1; big += abs(math.log(r[2] / cur[a])) > 0.0953
        if apply:
            os.makedirs(bak, exist_ok=True)
            if not os.path.exists(os.path.join(bak, f)): shutil.copy2(dst, os.path.join(bak, f))
            with open(dst + ".tmp", "w") as o:
                for r in new: o.write("\t".join(repr(x) for x in r) + "\n")
            os.replace(dst + ".tmp", dst)
        n += 1
    print("%s: %d tables %s; %d of %d rows differ from the current table by > 10 %% (SG)%s" % (
        T, n, "written" if apply else "would be written", big, rows, ("; backup in " + bak) if apply else ""))
