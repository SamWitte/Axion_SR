#!/usr/bin/env python3
"""Completeness of the rate campaign: channels complete on the alpha grid, partial, or without points; which level
families the gaps involve.  python3 completeness.py > completeness.txt ; also writes incomplete_channels.txt"""
import os, math, collections, statistics
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
ALPH = [float(x) for x in open("sgx_alphas.txt").read().strip().split(",")]
def st(tok): return tuple(int(x) for x in tok.split("-")) if "-" in tok else tuple(int(c) for c in tok)
def amax(c): return min(0.95 * min(st(s)[2] for s in c[:3]) / (2 * (1 + math.sqrt(1 - 0.95**2))) * 1.03, 1.7)
ch = [tuple(ln.split()[:4]) for ln in open("sgx_channels.txt") if ln.strip()]
ready = set(tuple(ln.split()[:4]) for ln in open("sgx_channels_ready.txt") if ln.strip())
full = part = none = 0; missing_pts = 0; fam = collections.Counter(); partial_amax = []; lowmiss = []; inc = []
for c in ch:
    need = [a for a in ALPH if a <= amax(c)]
    f = f"rates_out/{'_'.join(c)}.pts"; have = set()
    if os.path.isfile(f):
        for ln in open(f):
            p = ln.split()
            if len(p) >= 9: have.add(round(float(p[0]), 10))
    n = len([a for a in need if round(a, 10) in have])
    if n == len(need): full += 1; continue
    missing_pts += len(need) - n; inc.append((c, n, len(need)))
    if n == 0: none += 1
    else: part += 1; partial_amax.append(max(have))
    lv = [st(s) for s in c[:3]]
    fam[(c[3], "n>20" if max(s[0] for s in lv) > 20 else "n<=20", tuple(sorted(set(s[1:] for s in lv if s[1] >= 10))))] += 1
    if max(s[0] for s in lv) <= 20 and max(s[1] for s in lv) <= 8: lowmiss.append(c)
print(f"channels {len(ch)}: complete on all alphas {full}, partial {part}, no points {none}; missing points {missing_pts}; in run list {len(ready)}")
print("incomplete channels with every level n <= 20 and l <= 8:", len(lowmiss), lowmiss[:10])
print("incomplete channels by (final, n range, l >= 10 families involved):")
for k, v in fam.most_common(15): print("   ", v, k)
if partial_amax: print("partial channels: highest alpha present min/median/max:", min(partial_amax), statistics.median(partial_amax), max(partial_amax))
with open("incomplete_channels.txt", "w") as f:
    for c, n, need in inc: f.write(" ".join(c) + f"   {n}/{need} points\n")
