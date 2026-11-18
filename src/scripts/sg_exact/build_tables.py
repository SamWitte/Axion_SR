"""Fully relativistic self-gravity tables from the exact-exchange results (scan_gen.jl) and the Newtonian + recoil
tables of rate_sve.

For a channel computed at alpha0 (amplitudes a_SI, a_NR (Newton + recoil), a_GR (exact) per multipole l4 of mode 4):
    C_GR      = 64 alpha^4 sum|a_GR|^2 / sum|a_SI|^2            (Gamma_SG / Gamma_SI at f_a = M_pl)
    Gamma_SG  = C_GR * Gamma_SI(table)
    Gamma_x   = 16 alpha^2 sum Re(a_SI conj a_GR) / sum|a_SI|^2 * Gamma_SI(table)
    kappa     = Gamma_x / (2 sqrt(Gamma_SI Gamma_SG))
Gamma_SI is the table's (unchanged).  The other alpha nodes of the table are rescaled with the factors
F = Gamma_SG(exact) / Gamma_SG(table) and G = Gamma_x(exact) / Gamma_x(table) of the computed nodes:
    F(alpha) = 1 + c alpha^2, c interpolated linearly in alpha between computed nodes and constant outside
(relative corrections are O(alpha^2)); if a computed F or G is outside (0.5, 2) the factor of the nearest computed
node is used unchanged.  Only the nodes at the computed alphas are exact.

usage: build_tables.py <res dir> [<res dir> ...] --out <overlay dir> [--summary file]
"""
import glob, os, sys, math, collections
import numpy as np

RS = '/data/dust/user/samwitte/Axion_SR_sg/src/rate_sve/'
_aff = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'affected.txt')
AFFECTED = {tuple(l.split()[:4]) for l in open(_aff)} if os.path.exists(_aff) else set()


def read_results(dirs):
    """{(S1, S2, S3): {alpha: dict}} for the default numerical settings (drho 0.02, Lextra 6, rmax_fac 1)."""
    out = collections.defaultdict(dict)
    for d in dirs:
        for fn in glob.glob(os.path.join(d, 'res', '*.txt')):
            rows = [l.split() for l in open(fn) if l.strip()]
            if not rows:
                continue
            r0 = rows[0]
            if (float(r0[5]), int(r0[6]), float(r0[7])) != (0.02, 6, 1.0):
                continue
            al = float(r0[3])
            # points whose graviton source was cut too early by the first version of setup_gen (affected.txt, written by
            # affected_rcut.py) are only valid when recomputed with the corrected cut (version token rc2)
            if tuple(r0[:4]) in AFFECTED and 'rc2' not in r0:
                continue
            # channels with a dipole exchange (|m_b - m_3| = 1) are only converged with the larger angular basis of the
            # dipole solver (trailing field Lextra_dip >= 10, written by scan_gen.jl since 2026-11-16 16:30)
            m = [int(x.split('-')[2]) if '-' in x else int(x[2]) for x in r0[:3]]
            if (abs(m[0] - m[2]) == 1 or abs(m[1] - m[2]) == 1) and not (len(r0) >= 25 and float(r0[24]) >= 10):
                continue
            z = lambda r, i: complex(float(r[i]), float(r[i + 1]))
            aSI = np.array([z(r, 11) for r in rows]); aN = np.array([z(r, 13) for r in rows])
            aNR = np.array([z(r, 15) for r in rows]); aGR = np.array([z(r, 17) for r in rows])
            if not np.all(np.isfinite(aGR)):
                continue
            S_SI = np.sum(abs(aSI) ** 2)
            e = dict(alpha=al, dest=r0[8],
                     C_NR=64 * al ** 4 * np.sum(abs(aNR) ** 2) / S_SI, C_GR=64 * al ** 4 * np.sum(abs(aGR) ** 2) / S_SI,
                     X_NR=16 * al ** 2 * np.sum((aSI * np.conj(aNR)).real) / S_SI,
                     X_GR=16 * al ** 2 * np.sum((aSI * np.conj(aGR)).real) / S_SI,
                     nr_over_n=math.sqrt(np.sum(abs(aNR) ** 2) / np.sum(abs(aN) ** 2)),
                     rho0=aGR[0] / aNR[0], l4share=abs(aGR[0]) ** 2 / np.sum(abs(aGR) ** 2),
                     zmax=float(r0[19]), mode_res=float(r0[20]), mismatch=float(r0[21]), shift=float(r0[22]), sec=float(r0[23]))
            out[(r0[0], r0[1], r0[2])][al] = e
    return out


def factor_model(alphas, vals):
    """Factor at any alpha from the computed (alpha, factor) points, see the module docstring."""
    alphas = np.array(alphas); vals = np.array(vals)
    o = np.argsort(alphas); alphas = alphas[o]; vals = vals[o]
    pert = np.all((vals > 0.5) & (vals < 2.0))
    def f(al):
        if pert:
            c = np.interp(al, alphas, (vals - 1) / alphas ** 2)
            return 1 + c * al ** 2
        return vals[np.argmin(abs(alphas - al))]
    return f


def main():
    args = sys.argv[1:]
    out = args[args.index('--out') + 1]
    summ = args[args.index('--summary') + 1] if '--summary' in args else os.path.join(out, 'SUMMARY.tsv')
    dirs = [a for i, a in enumerate(args) if not a.startswith('--') and (i == 0 or args[i - 1] not in ('--out', '--summary'))]
    res = read_results(dirs)
    os.makedirs(out, exist_ok=True)
    lines = []; nfiles = 0
    for (s1, s2, s3), pts in sorted(res.items()):
        name = f'{s1}_{s2}_{s3}'
        files = [RS + f'{name}_{d}_LvrHc_SG.dat' for d in ('BH', 'Inf') if os.path.exists(RS + f'{name}_{d}_LvrHc_SG.dat')]
        if not files:
            continue
        tab = np.loadtxt(files[0], ndmin=2)
        comp = {}
        for al, e in pts.items():
            k = np.argmin(abs(tab[:, 0] - al))
            if abs(tab[k, 0] - al) > 1e-6 * al or tab[k, 2] <= 0:
                continue
            gsi, gsg, gx = tab[k, 1], tab[k, 2], tab[k, 3]
            F = e['C_GR'] * gsi / gsg
            G = e['X_GR'] * gsi / gx if gx != 0 else float('nan')
            comp[al] = (F, G, e, gsg / gsi)
        if not comp:
            continue
        als = sorted(comp)
        fF = factor_model(als, [comp[a][0] for a in als])
        Gs = [comp[a][1] for a in als]
        fG = factor_model(als, Gs) if np.all(np.isfinite(Gs)) else None
        for fn in files:
            t = np.loadtxt(fn, ndmin=2).copy()
            for i in range(len(t)):
                al, gsi, gsg, gx, kap = t[i]
                if gsg <= 0:
                    continue
                k = [a for a in als if abs(a - al) <= 1e-6 * al]
                if k:                                   # exact node
                    e = comp[k[0]][2]
                    nsg = e['C_GR'] * gsi; nx = e['X_GR'] * gsi
                else:
                    nsg = gsg * fF(al); nx = gx * (fG(al) if fG else 1.0)
                nk = nx / (2 * math.sqrt(gsi * nsg)) if gsi > 0 and nsg > 0 else kap
                t[i] = [al, gsi, nsg, nx, max(-1.0, min(1.0, nk))]
            tmp = os.path.join(out, os.path.basename(fn) + '.tmp')
            np.savetxt(tmp, t, delimiter='\t')
            os.replace(tmp, os.path.join(out, os.path.basename(fn)))
            nfiles += 1
        for al in als:
            F, G, e, ctab = comp[al]
            lines.append('\t'.join(str(x) for x in (
                s1, s2, s3, e['dest'], al, f'{F:.6g}', f'{G:.6g}', f"{e['C_GR'] / e['C_NR']:.6g}", f"{e['C_NR'] / ctab:.6g}",
                f"{e['nr_over_n']:.4g}", f"{e['rho0'].real:.6g}", f"{e['rho0'].imag:.3g}", f"{e['l4share']:.4g}", f"{e['zmax']:.2g}",
                f"{e['mismatch']:.2g}", f"{e['sec']:.0f}", len(files))))
    with open(summ, 'w') as f:
        f.write('# S1 S2 S3 dest alpha F_SG(exact/table) F_x(exact/table) exact/NR(same code) NR(code)/table |A_NR/A_N| '
                'Re,Im A_GR/A_NR(leading l4) share_of_leading_l4 lorenz_resid mode_mismatch seconds n_files\n')
        f.write('\n'.join(lines) + '\n')
    print(f'{len(res)} channels with results, {len(lines)} exact nodes, {nfiles} tables written to {out}')


if __name__ == '__main__':
    main()
