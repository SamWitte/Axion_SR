"""Levels populated in the SG runs and the channels among them.  usage: select_levels.py [threshold ...]"""
import glob, os, sys, re, collections, pickle
import numpy as np
BASE = '/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext/'
CH = '/data/dust/user/samwitte/Axion_SR_sg/src/rate_sve/sgx_channels.txt'
def st(s): return tuple(int(x) for x in s.split('-')) if '-' in s else (int(s[0]), int(s[1]), int(s[2]))
def peaks():
    cache = BASE + 'sg_exact/peaks.pkl'
    if os.path.exists(cache): return pickle.load(open(cache, 'rb'))
    out = {}
    for d in sorted(glob.glob(BASE + 'runs/full/gw_*/M_*/SG_*/alpha_*')):
        gw, M, var, al = d.split('/')[-4:]
        g = lambda k: glob.glob(d + '/' + k + '_*')[0]
        modes = np.loadtxt(g('Modes')).astype(int); U = np.loadtxt(g('States'))
        if U.shape[0] != len(modes): U = U.T
        pk = U.max(axis=1)
        out[(gw, M, var, float(al.split('_')[1]))] = {tuple(modes[i, :3]): pk[i] for i in range(len(modes)) if pk[i] > 1e-30}
    pickle.dump(out, open(cache, 'wb'))
    return out
if __name__ == '__main__':
    P = peaks()
    chans = [(st(p[0]), st(p[1]), st(p[2]), p[3]) for p in (l.split() for l in open(CH)) if len(p) >= 4]
    print('runs', len(P))
    for thr in [float(x) for x in sys.argv[1:]] or [1e-6, 1e-8, 1e-10]:
        byal = collections.defaultdict(set)
        for (gw, M, var, al), pk in P.items():
            byal[al] |= {s for s, v in pk.items() if v > thr}
        allv = set().union(*byal.values())
        tot = 0; pts = 0; uni = set()
        print(f'threshold {thr:g}: levels in union {len(allv)}, max n {max(s[0] for s in allv)}')
        for al in sorted(byal):
            L = byal[al]
            sel = [c for c in chans if c[0] in L and c[1] in L and c[2] in L]
            uni |= set(sel); pts += len(sel)
            print(f'   alpha {al}: levels {len(L)}  channels {len(sel)}  identical-in {sum(c[0]==c[1] for c in sel)}  max n {max(s[0] for s in L)}  max l {max(s[1] for s in L)}')
        selU = [c for c in chans if c[0] in allv and c[1] in allv and c[2] in allv]
        print(f'   per-alpha channel points {pts}; distinct channels {len(uni)}; channels among the union of levels {len(selU)}')
