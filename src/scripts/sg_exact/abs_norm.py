"""Absolute normalisation of the exact-exchange amplitudes against the production tables.
Production (gf_radial): Gamma = c(alpha) * Kp * sum_l4 |amp|^2 with
   Inf: Kp = 2 alpha k r_end^2,  r_end = 30 / k,  k = sqrt(w4^2 - alpha^2)
   BH : Kp = 4 alpha (1 + sqrt(1 - a^2)) * (w4 if m4 == 0 else 1)
c(alpha) is the same for every channel.  Here c = Gamma_SI(table) / (Kp sum |a_SI|^2) per computed node."""
import glob, os, sys, math, collections
import numpy as np
import build_tables as bt
from energies import energies
H = os.path.dirname(os.path.abspath(__file__)); RS = bt.RS; SPIN = 0.95
n_of = lambda s: int(s.split('-')[0]) if '-' in s else int(s[0])
def nodes(dirs=('prod', 'prod_dip', 'prod2'), need_gr=True):
    """Computed channel points.  need_gr=False also returns the points of the Newtonian-only runs (scan_gen.jl ...
    newton_only, no graviton solve: a['GR'] is NaN, has_gr False); the first directory that has a point wins."""
    E = energies(refresh=not need_gr); out = []
    seen = set()
    for d in dirs:
        for fn in glob.glob(os.path.join(H, d, 'res', '*.txt')):
            rows = [l.split() for l in open(fn) if l.strip()]
            if not rows: continue
            r0 = rows[0]
            if (float(r0[5]), int(r0[6]), float(r0[7])) != (0.02, 6, 1.0): continue
            m = [int(x.split('-')[2]) if '-' in x else int(x[2]) for x in r0[:3]]
            has_gr = all(math.isfinite(float(r[17])) and math.isfinite(float(r[18])) for r in rows)
            if has_gr:      # validity of the graviton solve (see build_tables.py)
                if tuple(r0[:4]) in bt.AFFECTED and 'rc2' not in r0: continue
                if (abs(m[0] - m[2]) == 1 or abs(m[1] - m[2]) == 1) and not (len(r0) >= 25 and float(r0[24]) >= 10): continue
            elif need_gr: continue
            al = float(r0[3]); key = (r0[0], r0[1], r0[2], al)
            if key in seen: continue
            try: w4 = E[(r0[0], al)] + E[(r0[1], al)] - E[(r0[2], al)]
            except KeyError: continue
            z = lambda r, i: complex(float(r[i]), float(r[i + 1]))
            a = {k: np.array([z(r, i) for r in rows]) for k, i in (('SI', 11), ('N', 13), ('NR', 15), ('GR', 17))}
            if not np.all(np.isfinite(a['NR'])) or not np.all(np.isfinite(a['SI'])): continue
            dest = r0[8]; m4 = int(r0[9])
            if dest == 'Inf':
                if w4 <= al: continue
                k = math.sqrt(w4 * w4 - al * al); Kp = 2 * al * k * (30 / k) ** 2
            else:
                Kp = 4 * al * (1 + math.sqrt(1 - SPIN ** 2)) * (w4 if m4 == 0 else 1.0)
            f = RS + f'{r0[0]}_{r0[1]}_{r0[2]}_{dest}_LvrHc_SG.dat'
            if not os.path.exists(f): continue
            t = np.loadtxt(f, ndmin=2); i = np.argmin(abs(t[:, 0] - al))
            if abs(t[i, 0] - al) > 1e-6 * al: continue
            seen.add(key)
            out.append(dict(key=key, dest=dest, m4=m4, w4=w4, Kp=Kp, a=a, tab=t[i], nmax=max(n_of(s) for s in r0[:3]), has_gr=has_gr,
                            dip=abs(m[0] - m[2]) == 1 or abs(m[1] - m[2]) == 1))
    return out
if __name__ == '__main__':
    N = nodes()
    print(len(N), 'nodes')
    for al in sorted({x['key'][3] for x in N}):
        for dest in ('BH', 'Inf'):
            s = [x for x in N if x['key'][3] == al and x['dest'] == dest and x['nmax'] <= 8 and x['tab'][1] > 0]
            if not s: continue
            c = np.array([x['tab'][1] / (x['Kp'] * np.sum(abs(x['a']['SI']) ** 2)) for x in s])
            print(f'alpha {al:<6g} {dest:3s} n<=8: {len(s):3d} nodes  c median {np.median(c):.6e}  min/med {c.min() / np.median(c):.4f}  max/med {c.max() / np.median(c):.4f}')
