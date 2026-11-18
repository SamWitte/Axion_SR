"""Job list for the Kerr-state Newtonian (+ recoil) and contact amplitudes (scan_gen.jl ... newton_only): every channel
of sgx_channels.txt whose three levels are populated above <threshold> in some self-gravity run at that alpha, at the
run's alpha node and at the next lower node of the table grid (the runs drift down in alpha as the hole loses mass).
Points that already have an exact result are left out.  Ordered by priority (smallest of the three peaks, largest first).
usage: make_jobs_newton.py <threshold> <out file>"""
import os, sys, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
thr = float(sys.argv[1]); outf = sys.argv[2]; sys.argv = sys.argv[:1]
from select_levels import peaks, st, CH
import build_tables as bt
RS = bt.RS
GRID = [0.03, 0.04, 0.05, 0.062, 0.075, 0.09, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.28]
fmt = lambda s: "%d%d%d" % s if s[0] < 10 else "%d-%d-%d" % s
P = peaks()
chans = [(st(p[0]), st(p[1]), st(p[2]), p[3]) for p in (l.split() for l in open(CH)) if len(p) >= 4]
byal = collections.defaultdict(dict)
for (gw, M, var, al), pk in P.items():
    for s, v in pk.items(): byal[al][s] = max(byal[al].get(s, 0.0), v)
H = os.path.dirname(os.path.abspath(__file__))
have = bt.read_results([os.path.join(H, d) for d in ('prod', 'prod_dip', 'prod2')])
jobs = {}
for al, pk in byal.items():
    L = {s for s, v in pk.items() if v > thr}
    lower = max([g for g in GRID if g < al - 1e-9], default=None)
    for c in chans:
        if c[0] in L and c[1] in L and c[2] in L:
            name = '_'.join(fmt(x) for x in c[:3])
            if not os.path.exists(RS + name + '_' + c[3] + '_LvrHc_SG.dat'): continue
            pr = min(pk[c[0]], pk[c[1]], pk[c[2]])
            for a2 in (al, lower):
                if a2 is None: continue
                key = (fmt(c[0]), fmt(c[1]), fmt(c[2]), a2)
                if a2 in have.get(key[:3], {}): continue
                if key not in jobs or pr > jobs[key]: jobs[key] = pr
rows = sorted(jobs.items(), key=lambda kv: (-kv[1], kv[0][3]))
with open(outf, 'w') as f:
    for (s1, s2, s3, al), pr in rows: f.write(f'{s1} {s2} {s3} {al:g}\n')
print(len(rows), 'jobs;', dict(sorted(collections.Counter(k[3] for k, _ in rows).items())))
