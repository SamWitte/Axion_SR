#!/usr/bin/env python3
"""Highest populated level of the (n,4,4) tower per run and table set.  usage: tower_extent.py <set> [<set> ...]
Prints, per run, the highest n with peak occupation > 1e-10 and (in brackets) > 1e-20."""
import sys, runsum
names = sys.argv[1:] or ['full', 'x2_gr', 'x2_kerr_nr', 'x2_kerr_gr']
S = {n: runsum.summary(n) for n in names}
def tower(r, thr):
    if r is None or r['t'] < 0.99e9: return None
    ns = [s[0] for s, p in r['pk'].items() if s[1] == 4 and s[2] == 4 and p > thr]
    return max(ns) if ns else 0
print('# highest n of the (n,4,4) tower with peak occupation > 1e-10 [> 1e-20]; sets: ' + ' | '.join(names))
diff = tot = 0
for k in sorted(S[names[0]], key=lambda k: (k[0], float(k[1]), k[2], float(k[3]))):
    if 'noSG' in k[2]: continue
    t = [(tower(S[n].get(k), 1e-10), tower(S[n].get(k), 1e-20)) for n in names]
    if all(x[0] is None for x in t[1:]): continue
    tot += 1
    flag = len({x[0] for x in t if x[0] is not None}) > 1
    diff += flag
    print(f"{k[0]:6s} M={k[1]:<6s} {k[2]:14s} alpha {k[3]:<5s} " + ' | '.join(' -     ' if x[0] is None else f'{x[0]:2d} [{x[1]:2d}]' for x in t) + ('   <-- differs' if flag else ''))
print(f'# {diff} of {tot} runs differ in the highest tower level above 1e-10')
