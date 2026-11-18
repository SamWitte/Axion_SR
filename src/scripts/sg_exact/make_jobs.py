"""Job list for the exact self-gravity exchange: every channel of sgx_channels.txt whose three levels are populated
above a threshold in some self-gravity run at that alpha, ordered by priority (the smallest of the three peak
occupations, largest first).  Lines: S1 S2 S3 alpha.   usage: make_jobs.py <threshold> <out prefix>"""
import os, sys, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from select_levels import peaks, st, CH
RS = '/data/dust/user/samwitte/Axion_SR_sg/src/rate_sve/'
fmt = lambda s: "%d%d%d" % s if s[0] < 10 else "%d-%d-%d" % s
ERG = lambda n: 1 - 0.1 ** 2 / (2 * n * n)          # hydrogenic, alpha = 0.1 (the rule of the channel lists)

def nc_piece(b, s3, a, l4max, Lextra=6):
    Lmax = min(b[1] + s3[1], a[1] + l4max) + Lextra
    mg = abs(b[2] - s3[2])
    return sum(max(0, Lmax - max(s, mg) + 1) for s in (0, 0, 0, 0, 1, 1, 1, 1, 2, 2)), Lmax

def cost(c):
    s1, s2, s3, dest = c
    l4max = 2 if dest == 'BH' else max(s1[1] + s2[1] - s3[1], abs(s1[2] + s2[2] - s3[2]))
    pcs = [nc_piece(s2, s3, s1, l4max)] + ([] if s1 == s2 else [nc_piece(s1, s3, s2, l4max)])
    ncmax = max(p[0] for p in pcs)
    return ncmax, sum(p[0] ** 2 for p in pcs)

if __name__ == '__main__':
    thr = float(sys.argv[1]); pre = sys.argv[2]
    P = peaks()
    chans = [(st(p[0]), st(p[1]), st(p[2]), p[3]) for p in (l.split() for l in open(CH)) if len(p) >= 4]
    byal = collections.defaultdict(dict)
    for (gw, M, var, al), pk in P.items():
        for s, v in pk.items():
            byal[al][s] = max(byal[al].get(s, 0.0), v)
    jobs = {}
    for al, pk in byal.items():
        L = {s for s, v in pk.items() if v > thr}
        for c in chans:
            if c[0] in L and c[1] in L and c[2] in L:
                name = '_'.join(fmt(x) for x in c[:3])
                if not os.path.exists(RS + name + '_' + c[3] + '_LvrHc_SG.dat'):
                    continue                                    # no Newtonian table (zero or incomplete channel)
                key = (c[0], c[1], c[2], al)
                pr = min(pk[c[0]], pk[c[1]], pk[c[2]])
                if key not in jobs or pr > jobs[key][0]:
                    jobs[key] = (pr, c)
    rows = sorted(jobs.items(), key=lambda kv: -kv[1][0])
    tiers = collections.defaultdict(list)
    for (s1, s2, s3, al), (pr, c) in rows:
        ncmax, nc2 = cost(c)
        t = 'A' if ncmax <= 130 else 'B' if ncmax <= 190 else 'C'
        tiers[t].append((fmt(s1), fmt(s2), fmt(s3), al, pr, ncmax, nc2))
    for t, v in sorted(tiers.items()):
        with open(f'{pre}_{t}.txt', 'w') as f:
            f.write('# S1 S2 S3 alpha   (priority, NC_max of the two graviton solves as trailing comment)\n')
            for r in v:
                f.write(f'{r[0]} {r[1]} {r[2]} {r[3]:g}\n')
        print(f'tier {t}: {len(v)} jobs, sum NC^2 / 110^2 = {sum(r[6] for r in v) / 110 ** 2:.0f} (baseline-solve units); NC range {min(r[5] for r in v)}-{max(r[5] for r in v)}')
    with open(f'{pre}_all.tsv', 'w') as f:
        for t, v in sorted(tiers.items()):
            for r in v:
                f.write('\t'.join(str(x) for x in (t,) + r) + '\n')
    print('total', len(rows), 'by alpha', dict(collections.Counter(k[3] for k, _ in rows)))
