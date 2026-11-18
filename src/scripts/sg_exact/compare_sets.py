#!/usr/bin/env python3
"""Run-by-run comparison of run sets against a reference set (all with the runsum layout).
    python3 compare_sets.py <reference set> <set> [<set> ...]
Per set: largest change of final spin and final M/M_0, runs above 0.001, runs whose (n,4,4) tower extent
(highest n with peak > 1e-10) differs, and the changes of the peaks of 211, 322, 544."""
import sys
import numpy as np
import runsum
ref = runsum.summary(sys.argv[1])
def tower(r, thr=1e-10):
    ns = [s[0] for s, p in r["pk"].items() if s[1] == 4 and s[2] == 4 and p > thr]
    return max(ns) if ns else 0
for name in sys.argv[2:]:
    S = runsum.summary(name)
    keys = sorted(k for k in ref if k in S and ref[k]["t"] >= 0.99e9 and S[k]["t"] >= 0.99e9)
    da = np.array([S[k]["a"] - ref[k]["a"] for k in keys]); dM = np.array([S[k]["M"] - ref[k]["M"] for k in keys])
    print(f"\n=== {name} against {sys.argv[1]}: {len(keys)} runs")
    print(f"  final spin:  max |change| {abs(da).max():.4f}; runs above 0.001: {(abs(da) > 1e-3).sum()}, above 0.0001: {(abs(da) > 1e-4).sum()}")
    print(f"  final M/M_0: max |change| {abs(dM).max():.4f}; runs above 0.001: {(abs(dM) > 1e-3).sum()}")
    for i in np.argsort(-abs(da))[:8]:
        k = keys[i]
        if abs(da[i]) < 5e-5: break
        print(f"    {k[0]:6s} M={k[1]:6s} {k[2]:15s} alpha {k[3]:5s}  a_f {ref[k]['a']:.4f} -> {S[k]['a']:.4f} ({da[i]:+.4f})   M_f/M_0 {ref[k]['M']:.4f} -> {S[k]['M']:.4f}")
    for lev in ((2, 1, 1), (3, 2, 2), (5, 4, 4)):
        r = np.array([S[k]["pk"].get(lev, 0) / ref[k]["pk"][lev] for k in keys if ref[k]["pk"].get(lev, 0) > 1e-12])
        if len(r): print(f"  peak of {lev[0]}{lev[1]}{lev[2]}: ratio median {np.median(r):.3f}, range {r.min():.3g} - {r.max():.3g}; changed by more than 10 % in {(abs(np.log(np.maximum(r, 1e-300))) > np.log(1.1)).sum()} of {len(r)} runs")
    tw = [(k, tower(ref[k]), tower(S[k])) for k in keys if tower(ref[k]) != tower(S[k])]
    print(f"  (n,4,4) tower extent differs in {len(tw)} runs ({sum(1 for t in tw if str(t[0][1]) in ('10', '10.0'))} at M = 10; {sum(1 for t in tw if abs(t[2] - t[1]) > 1)} by more than one level)")
    for k, a, b in tw: print(f"    {k[0]:6s} M={k[1]:6s} {k[2]:15s} alpha {k[3]:5s}  {a} -> {b}")
    nl = [(k, sum(1 for p in ref[k]['pk'].values() if p > 1e-10), sum(1 for p in S[k]['pk'].values() if p > 1e-10)) for k in keys]
    print(f"  number of levels above 1e-10 differs in {sum(1 for _, a, b in nl if a != b)} runs (largest change {max((abs(b - a) for _, a, b in nl), default=0)})")
