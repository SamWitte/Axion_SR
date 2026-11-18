"""Channel points whose graviton source was cut too early by the first version of setup_gen (exchange_gen.jl).

The source T_b3 was kept to rcut = max(34.5/(k_b + k_3), 50/k_a), which assumes pure exponentials.  For two high-n
levels b, 3 and a low-n level a this ends inside the region where psi_b psi_3* is still large (the levels oscillate
out to ~2 n^2 a0).  Here the needed radius is estimated with hydrogenic functions: the last radius where
r^2 |R_b R_3| (source) or r^2 |R_a| (acted-on level) exceeds 1e-10 of its maximum.
usage: affected_rcut.py <job files...>   -> prints the affected lines (S1 S2 S3 alpha)
"""
import sys
import numpy as np
st = lambda s: tuple(int(x) for x in s.split('-')) if '-' in s else (int(s[0]), int(s[1]), int(s[2]))
R = np.logspace(-1, np.log10(6000.0), 6000)          # radii in a0 = GM / alpha^2
_cache = {}
def radial(n, l):
    if (n, l) not in _cache:
        x = 2 * R / n; k = n - l - 1; a = 2 * l + 1
        L0 = np.ones_like(x); L1 = 1 + a - x
        if k == 0: L = L0
        elif k == 1: L = L1
        else:
            for j in range(1, k):
                L0, L1 = L1, ((2 * j + 1 + a - x) * L1 - (j + a) * L0) / (j + 1)
            L = L1
        _cache[(n, l)] = np.abs(R ** l * np.exp(-R / n) * L)
    return _cache[(n, l)]
def last_above(v, tol=1e-10):
    w = v * R ** 2
    return R[min(np.where(w > tol * w.max())[0][-1] + 1, len(R) - 1)]
def pieces(s1, s2, s3):
    return [(s2, s3, s1)] + ([] if s1 == s2 else [(s1, s3, s2)])
def affected(s1, s2, s3):
    worst = 1.0
    for b, t, a in pieces(s1, s2, s3):
        old = max(34.5 * b[0] * t[0] / (b[0] + t[0]), 50.0 * a[0])
        need = max(last_above(radial(b[0], b[1]) * radial(t[0], t[1])), last_above(radial(a[0], a[1])))
        worst = max(worst, need / old)
    return worst
if __name__ == '__main__':
    seen = set()
    for f in sys.argv[1:]:
        for l in open(f):
            p = l.split()
            if not p or p[0].startswith('#') or tuple(p[:4]) in seen: continue
            seen.add(tuple(p[:4]))
            w = affected(st(p[0]), st(p[1]), st(p[2]))
            if w > 1.0: print(' '.join(p[:4]), f'{w:.2f}')
