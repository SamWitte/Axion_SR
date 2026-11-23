#!/usr/bin/env python3
"""Completeness of the level energies (levels_out/lev_*.txt): exit 1 if any (level, alpha, spin) row is missing."""
import glob, sys, collections
D = "/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext/"
SP = "0,0.3,0.5,0.7,0.8,0.9,0.95,0.99,0.998".split(","); AL = open(D + "sgx_alphas.txt").read().strip().split(",")
lv = [l.strip().replace("-", " ") for l in open(D + "spin_tab/tasks/levels.txt") if l.strip()]
have = {}
for f in glob.glob(D + "spin_tab/levels_out/lev_*.txt"):
    for ln in open(f):
        p = ln.split()
        if len(p) == 9: have[" ".join(p[:5])] = int(p[7])
want = ["%s %g %g" % (l, float(a), float(al)) for l in lv for a in SP for al in AL]
miss = [k for k in want if k not in have]; bad = collections.Counter(k.split()[3] for k in want if have.get(k) == 0)
print("level energies: %d of %d present, %d missing; eigenvalue search failed (hydrogenic fallback) per spin: %s" % (len(want) - len(miss), len(want), len(miss), dict(bad)))
if miss: print("first missing:", miss[:5])
sys.exit(1 if miss else 0)
