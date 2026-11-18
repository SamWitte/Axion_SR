"""Kerr level energies G M omega of the exact-exchange jobs, parsed from the worker logs ('=== S1 S2 S3 alpha' followed
by 'modes: w = w1 w2 w3').  energies() -> {(state, alpha): w}"""
import glob, os, re, pickle
H = os.path.dirname(os.path.abspath(__file__))
def energies(refresh=False):
    cache = os.path.join(H, 'energies.pkl')
    if os.path.exists(cache) and not refresh:
        return pickle.load(open(cache, 'rb'))
    E = {}
    for fn in glob.glob(os.path.join(H, 'prod*', 'logs', '*.log')) + glob.glob(os.path.join(H, 'test', '**', '*.log'), recursive=True):
        cur = None
        for line in open(fn, errors='replace'):
            if line.startswith('=== '):
                p = line.split(); cur = p[1:5] if len(p) >= 5 else None
            elif line.startswith('modes: w =') and cur:
                w = line.split()[3:6]
                try:
                    al = float(cur[3]); ws = [float(x) for x in w]
                except ValueError:
                    cur = None; continue
                for s, x in zip(cur[:3], ws):
                    if 0.5 * al < x < al: E[(s, al)] = x
                cur = None
    pickle.dump(E, open(cache, 'wb'))
    return E
if __name__ == '__main__':
    E = energies(refresh=True)
    import collections
    print(len(E), 'level energies;', dict(collections.Counter(a for _, a in E)))
