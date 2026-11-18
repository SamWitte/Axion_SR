"""Overlay self-gravity tables from the exact-exchange results (scan_gen.jl), in three variants that keep apart the two
things in which the exact code differs from the production tables of rate_sve:

  (1) the exchange: Newtonian potential + BH recoil (production)  ->  linearised Einstein equations on Kerr;
  (2) the bound states: production uses hydrogenic radial functions for n > 8, the exact code Kerr states for all n.

  --mode gr       only (1).  Gamma_SI as tabulated;  Gamma_SG = Gamma_SG(table) * sum|a_GR|^2 / sum|a_NR|^2;
                  kappa = kappa(table) * kappa_GR / kappa_NR  (a_NR, a_GR from the same code with the same states).
  --mode kerr_nr  only (2).  All three columns from the code with the Newtonian + recoil amplitudes, in absolute units.
  --mode kerr_gr  (1) and (2).  All three columns from the code with the exact amplitudes.

Absolute units (abs_norm.py): Gamma = c * sym * Kp * sum_l4 |a|^2, Kp from gf_radial's rate formula, sym = 4 for two
different incoming states, c = median of Gamma_SI(table) / (sym Kp sum|a_SI|^2) over the computed channels with all
n <= 8 at the same alpha and destination (those agree to 0.5 %).  This holds for Inf channels and BH channels with
m4 = 0.  For BH channels with m4 != 0 the absolute normalisation is not reproduced (step-size dependent horizon
amplitude), so there the kerr modes keep Gamma_SI(table) and set Gamma_SG = Gamma_SI(table) * 64 alpha^4 sum|a|^2/sum|a_SI|^2.

Nodes that have the Newtonian amplitudes only (scan_gen.jl ... newton_only; directory prodN) enter kerr_nr as they are;
in kerr_gr their Gamma_SG is multiplied by the exact/Newtonian factor of the same channel's nodes with a graviton solve,
carried over as 1 + c alpha^2 (limited to the range spanned by 1 and the computed factors), if the channel has any.

Rows at alphas that were not computed: every column is multiplied by f(alpha) = 1 + c(alpha) alpha^2 (c interpolated
linearly between the computed nodes, constant outside; f limited to the range spanned by 1 and the computed factors) if all computed factors of that column lie in (0.5, 2); otherwise
the row is left as tabulated (large factors come from zeros of an amplitude, which are local in alpha).  kappa of such
rows is left as tabulated.  Only the computed nodes are exact.

usage: build_tables2.py --mode gr|kerr_nr|kerr_gr --out <overlay dir> [res dirs ...]      (default: prod prod_dip prod2)
"""
import os, sys, math, collections
import numpy as np
import abs_norm as an

RS = an.RS


def node_values(x, mode, cal, grmodel=None):
    """(Gamma_SI, Gamma_SG, kappa, absolute units?) of one computed node, or None.
    grmodel = (f_SG(alpha), kappa ratio): exact / Newtonian of the same channel at its nodes with a graviton solve; used
    in mode kerr_gr for nodes that have the Newtonian amplitudes only."""
    al = x['key'][3]; a = x['a']; tab = x['tab']
    gsi_t, gsg_t, kap_t = tab[1], tab[2], tab[4]
    S_SI = np.sum(abs(a['SI']) ** 2); S_NR = np.sum(abs(a['NR']) ** 2)
    S_GR = np.sum(abs(a['GR']) ** 2) if x['has_gr'] else float('nan')
    if not (S_SI > 0 and S_NR > 0 and gsi_t > 0 and gsg_t > 0) or (x['has_gr'] and not S_GR > 0):
        return None
    kap = lambda S, k: float(np.sum((a['SI'] * np.conj(a[k])).real) / math.sqrt(S_SI * S))
    if mode == 'gr':
        if not x['has_gr']:
            return None
        kn, kg = kap(S_NR, 'NR'), kap(S_GR, 'GR')
        return gsi_t, gsg_t * S_GR / S_NR, (kap_t * kg / kn if abs(kn) > 1e-6 else kg), False
    c = cal.get((al, x['dest']), cal.get('any', 0.0) / al)
    absolute = c > 0 and (x['dest'] == 'Inf' or x['m4'] == 0)
    gsi = c * x['sym'] * x['Kp'] * S_SI if absolute else gsi_t
    if mode == 'kerr_nr' or (not x['has_gr'] and grmodel is None):
        return gsi, gsi * 64 * al ** 4 * S_NR / S_SI, kap(S_NR, 'NR'), absolute
    if x['has_gr']:
        return gsi, gsi * 64 * al ** 4 * S_GR / S_SI, kap(S_GR, 'GR'), absolute
    fsg, kr = grmodel
    return gsi, gsi * 64 * al ** 4 * S_NR / S_SI * fsg(al), kap(S_NR, 'NR') * kr, absolute


def factor_fn(alphas, vals, always=False):
    alphas = np.array(alphas); vals = np.array(vals)
    if not always and not np.all((vals > 0.5) & (vals < 2.0)):
        return lambda al: 1.0
    o = np.argsort(alphas); alphas = alphas[o]; cc = (vals[o] - 1) / alphas ** 2
    lo, hi = min(1.0, vals.min()), max(1.0, vals.max())          # never extrapolate beyond the computed factors
    return lambda al: min(hi, max(lo, 1 + float(np.interp(al, alphas, cc)) * al ** 2))


def main():
    args = sys.argv[1:]
    mode = args[args.index('--mode') + 1]; out = args[args.index('--out') + 1]
    assert mode in ('gr', 'kerr_nr', 'kerr_gr')
    dirs = [a for i, a in enumerate(args) if not a.startswith('--') and (i == 0 or args[i - 1] not in ('--mode', '--out'))] or ['prod', 'prod_dip', 'prod2']
    N = an.nodes(dirs, need_gr=(mode == 'gr'))
    for x in N:
        x['sym'] = 1.0 if x['key'][0] == x['key'][1] else 4.0
    calv = collections.defaultdict(list)
    for x in N:
        if x['nmax'] <= 8 and x['tab'][1] > 0 and (x['dest'] == 'Inf' or x['m4'] == 0):
            calv[(x['key'][3], x['dest'])].append(x['tab'][1] / (x['sym'] * x['Kp'] * np.sum(abs(x['a']['SI']) ** 2)))
    cal = {k: float(np.median(v)) for k, v in calv.items()}
    cal['any'] = float(np.median([v * k[0] for k, vv in calv.items() if k[0] >= 0.09 for v in vv]))      # c * alpha, universal
    chan = collections.defaultdict(list)
    for x in N:
        chan[x['key'][:3]].append(x)
    by = collections.defaultdict(dict); n_model = 0
    for ch, xs in chan.items():
        grmodel = None
        ex = [x for x in xs if x['has_gr']]
        if mode == 'kerr_gr' and ex and len(ex) < len(xs):
            S = lambda x, k: np.sum(abs(x['a'][k]) ** 2)
            ex = [x for x in ex if S(x, 'NR') > 0 and S(x, 'GR') > 0]
            if ex:
                fsg = factor_fn([x['key'][3] for x in ex], [S(x, 'GR') / S(x, 'NR') for x in ex], always=True)
                kp = lambda x, k: float(np.sum((x['a']['SI'] * np.conj(x['a'][k])).real) / math.sqrt(S(x, 'SI') * S(x, k)))
                x0 = max(ex, key=lambda x: x['key'][3])
                kn = kp(x0, 'NR'); kr = kp(x0, 'GR') / kn if abs(kn) > 1e-6 else 1.0
                grmodel = (fsg, kr if 0.5 < kr < 2 else 1.0); n_model += 1
        for x in xs:
            v = node_values(x, mode, cal, grmodel)
            if v is not None:
                by[ch][x['key'][3]] = (v, x)
    os.makedirs(out, exist_ok=True)
    lines = []; nfiles = 0
    for (s1, s2, s3), pts in sorted(by.items()):
        name = f'{s1}_{s2}_{s3}'
        files = [RS + f'{name}_{d}_LvrHc_SG.dat' for d in ('BH', 'Inf') if os.path.exists(RS + f'{name}_{d}_LvrHc_SG.dat')]
        als = sorted(pts)
        fSI = factor_fn(als, [pts[a][0][0] / pts[a][1]['tab'][1] for a in als])
        fSG = factor_fn(als, [pts[a][0][1] / pts[a][1]['tab'][2] for a in als])
        for fn in files:
            t = np.loadtxt(fn, ndmin=2).copy()
            for i in range(len(t)):
                al, gsi, gsg, gx, kap = t[i]
                if gsg <= 0 or gsi <= 0:
                    continue
                k = [a for a in als if abs(a - al) <= 1e-6 * al]
                if k:
                    nsi, nsg, nk, _ = pts[k[0]][0]
                else:
                    nsi, nsg, nk = gsi * fSI(al), gsg * fSG(al), kap
                nk = max(-1.0, min(1.0, nk))
                t[i] = [al, nsi, nsg, 2 * nk * math.sqrt(nsi * nsg), nk]
            tmp = os.path.join(out, os.path.basename(fn) + '.tmp')
            np.savetxt(tmp, t, delimiter='\t')
            os.replace(tmp, os.path.join(out, os.path.basename(fn)))
            nfiles += 1
        for al in als:
            (nsi, nsg, nk, absolute), x = pts[al]
            tab = x['tab']
            lines.append('\t'.join(str(v) for v in (s1, s2, s3, x['dest'], al, f'{nsi / tab[1]:.6g}', f'{nsg / tab[2]:.6g}', f'{tab[4]:.4f}',
                                                    f'{nk:.4f}', int(absolute), x['nmax'], int(x['dip']), len(files))))
    with open(os.path.join(out, 'SUMMARY.tsv'), 'w') as f:
        f.write(f'# mode {mode}\n# S1 S2 S3 dest alpha Gamma_SI(new/table) Gamma_SG(new/table) kappa(table) kappa(new) absolute_units n_max dipole n_files\n')
        f.write('\n'.join(lines) + '\n')
    print(f'mode {mode}: {len(by)} channels, {len(lines)} computed nodes ({sum(1 for x in N if x["has_gr"])} with a graviton solve, '
          f'{n_model} channels with the exact/Newtonian factor carried to Newtonian-only nodes), {nfiles} tables written to {out}; '
          f'calibration c*alpha = {cal["any"]:.4e}')


if __name__ == '__main__':
    main()
