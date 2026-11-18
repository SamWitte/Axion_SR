"""Like make_jobs_newton.py, for the levels populated in a given set of run directories (e.g. the SI-only runs):
every channel of sgx_channels.txt among the levels with peak occupation > threshold in a run, at the run's alpha node
and the next lower node, minus the points that already have a result in prod, prod_dip, prod2 or prodN.
usage: make_jobs_newton_runs.py <threshold> <out file> <run dir glob> [...]     (run dirs end in .../alpha_<a>)"""
import glob, os, sys, collections, warnings
import numpy as np
warnings.simplefilter('ignore', DeprecationWarning)
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
thr = float(sys.argv[1]); outf = sys.argv[2]; globs = sys.argv[3:]; sys.argv = sys.argv[:1]
from select_levels import st, CH
import build_tables as bt
RS = bt.RS
GRID = [0.03, 0.04, 0.05, 0.062, 0.075, 0.09, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.28]
fmt = lambda s: "%d%d%d" % s if s[0] < 10 else "%d-%d-%d" % s
byal = collections.defaultdict(dict)
for gl in globs:
    for d in sorted(glob.glob(gl)):
        al = float(os.path.basename(d).split('_')[1])
        g = lambda k: (glob.glob(d + '/' + k + '_*.dat') or [None])[0]
        if g('States') is None: continue
        modes = [tuple(int(float(x)) for x in l.split()[:3]) for l in open(g('Modes'))]
        with open(g('States')) as f:
            for i, line in enumerate(f):
                p = np.fromstring(line, sep='\t').max()
                if p > thr: byal[al][modes[i]] = max(byal[al].get(modes[i], 0.0), p)
have = set()
for d in ('prod', 'prod_dip', 'prod2', 'prodN'):
    have |= {f[:-4] for f in os.listdir(os.path.join(H, d, 'res')) if f.endswith('.txt')}
chans = [(st(p[0]), st(p[1]), st(p[2]), p[3]) for p in (l.split() for l in open(CH)) if len(p) >= 4]
jobs = {}
for al, pk in byal.items():
    near = min(GRID, key=lambda g: abs(g - al))
    nodes = [near] if abs(near - al) < 1e-9 else [max(g for g in GRID if g < al), min(g for g in GRID if g > al)]
    lower = max([g for g in GRID if g < min(nodes) - 1e-9], default=None)
    for c in chans:
        if c[0] in pk and c[1] in pk and c[2] in pk:
            name = '_'.join(fmt(x) for x in c[:3])
            if not os.path.exists(RS + name + '_' + c[3] + '_LvrHc_SG.dat'): continue
            pr = min(pk[c[0]], pk[c[1]], pk[c[2]])
            for a2 in nodes + [lower]:
                if a2 is None: continue
                if f'{name}_{a2:g}' in have: continue
                key = (fmt(c[0]), fmt(c[1]), fmt(c[2]), a2)
                if key not in jobs or pr > jobs[key]: jobs[key] = pr
    print(f'alpha {al}: {len(pk)} levels above {thr:g}, highest (n,4,4): {max([s[0] for s in pk if s[1:] == (4, 4)] or [0])}')
rows = sorted(jobs.items(), key=lambda kv: (-kv[1], kv[0][3]))
with open(outf, 'w') as f:
    for (s1, s2, s3, al), pr in rows: f.write(f'{s1} {s2} {s3} {al:g}\n')
print(len(rows), 'new jobs;', dict(sorted(collections.Counter(k[3] for k, _ in rows).items())))
