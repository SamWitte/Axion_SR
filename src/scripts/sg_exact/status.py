"""Which target points of the exact-exchange campaign have a valid result (as build_tables.py judges), and where the
missing ones stand.  usage: status.py [--list-missing file]"""
import os, sys, collections, time
import build_tables as bt
H = os.path.dirname(os.path.abspath(__file__))
def jobs(fn):
    return [tuple(l.split()[:4]) for l in open(os.path.join(H, fn)) if l.strip() and not l.startswith('#')]
target = {}
for fn in ('jobs_0base.txt', 'jobs_A.txt', 'jobs_B.txt', 'jobs_C.txt', 'jobs_dipole_all.txt'):
    for j in jobs(fn): target.setdefault(j, fn)
res = bt.read_results([os.path.join(H, d) for d in ('prod', 'prod_dip', 'prod2')])
have = {(s1, s2, s3, a) for (s1, s2, s3), pts in res.items() for a in pts}
def norm(j): return (j[0], j[1], j[2], float(j[3]))
todoA = set(map(norm, jobs('jobs_todoA.txt'))); todoB = set(map(norm, jobs('jobs_todoB.txt'))); dip = set(map(norm, jobs('jobs_dipole_all.txt')))
by = collections.Counter(); miss = []
now = time.time()
for j in target:
    k = norm(j); ok = k in have
    by[(k[3], ok)] += 1
    if not ok:
        idn = '_'.join(j)
        lock = [d for d in ('prod', 'prod_dip', 'prod2') if os.path.isdir(os.path.join(H, d, 'locks', idn))]
        age = [int((now - os.path.getmtime(os.path.join(H, d, 'locks', idn))) / 60) for d in lock]
        own = [open(os.path.join(H, d, 'locks', idn, 'owner')).read().split()[0] if os.path.exists(os.path.join(H, d, 'locks', idn, 'owner')) else '?' for d in lock]
        fail = [d for d in ('prod', 'prod_dip', 'prod2') if os.path.exists(os.path.join(H, d, 'fail', idn + '.txt'))]
        miss.append((j, 'dip' if k in dip else '', 'A' if k in todoA else ('B' if k in todoB else '-'), list(zip(lock, own, age)), fail))
for a in sorted({k[0] for k in by}):
    print(f'alpha {a}: valid {by[(a, True)]}, missing {by[(a, False)]}')
print('total valid', sum(v for (a, ok), v in by.items() if ok), 'missing', len(miss))
c = collections.Counter()
for j, d, t, lock, fail in miss:
    st = 'failed' if fail else ('running' if any(o != 'held_for_Lextra' and ag < 50 for _, o, ag in lock) else 'free')
    c[(d or 'plain', 'todo' + t, st)] += 1
for k, v in sorted(c.items()): print(k, v)
if '--list-missing' in sys.argv:
    with open(sys.argv[sys.argv.index('--list-missing') + 1], 'w') as f:
        for j, d, t, lock, fail in miss: f.write(' '.join(j) + f'\t{d}\t{t}\t{lock}\t{fail}\n')
