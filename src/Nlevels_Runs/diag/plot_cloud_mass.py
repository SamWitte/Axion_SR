#!/usr/bin/env python3
"""Cloud-mass figure for a set of finished runs (one column per run):
spin a(t) and M_BH(t)/M_0 in two small panels on top, and the fractional cloud mass
M_cloud/M_0 = alpha_0 * u_i per level (thin, coloured by level) and in total (black) on a
log-log panel below (the total is drawn underneath, so it shows as a black edge around the dominant level). u_i = N_i/(G M_0^2) are the saved occupations, so alpha_0 * sum(u)
equals the BH mass loss while the cloud only grows (checked on acc0_gwoff_fullbasis_v2).

  python3 plot_cloud_mass.py OUT.pdf NMAX LEAF[:LABEL] [LEAF[:LABEL] ...]
LEAF = output/<tag>/BH_x/fa_y/alpha_z holding NL_output_Nmax_<NMAX>.tar.gz; alpha_0 is read from the path.
All columns share the time axis and the spin and BH-mass ranges.
Also writes OUT with .png. No legend (user request); the 3 dominant levels of each run are labelled
at their peaks.
"""
import os
import re
import sys
import tarfile
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["cmr10", "Computer Modern Roman", "DejaVu Serif"],
    "mathtext.fontset": "cm", "axes.formatter.use_mathtext": True, "axes.unicode_minus": False,
    "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
    "xtick.major.size": 8, "ytick.major.size": 8, "xtick.minor.size": 4, "ytick.minor.size": 4,
    "xtick.major.width": 1.1, "ytick.major.width": 1.1, "xtick.minor.width": 0.8, "ytick.minor.width": 0.8,
    "xtick.labelsize": 13, "ytick.labelsize": 13, "axes.linewidth": 1.1,
})
LEVEL_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
OTHER = "#c9c8c2"
LABEL_FS = 17


def load(leaf, nmax):
    with tarfile.open(os.path.join(leaf, f"NL_output_Nmax_{nmax}.tar.gz")) as t:
        g = lambda k: [m for m in t.getmembers() if os.path.basename(m.name).startswith(k + "_")][0]
        T = np.array(t.extractfile(g("Time")).read().split(), float)
        S = np.array(t.extractfile(g("Spin")).read().split(), float)
        M = np.array(t.extractfile(g("MassBH")).read().split(), float)
        modes = [tuple(int(float(x)) for x in ln.split()[:3]) for ln in t.extractfile(g("Modes")).read().decode().splitlines() if ln.strip()]
        U = np.array([np.array(ln.split(), float) for ln in t.extractfile(g("States")) if ln.strip()])
    return T, S, M, modes, U


def main():
    out, nmax, specs = sys.argv[1], int(sys.argv[2]), sys.argv[3:]
    runs = []
    for sp in specs:
        leaf, _, lab = sp.partition(":")
        a0 = float(re.search(r"alpha_([0-9.]+)", leaf).group(1))
        runs.append((leaf, lab or rf"$\alpha = {a0:g}$", a0) + load(leaf, nmax))
    # colour follows the level: the levels that dominate any run, in order of first appearance by peak
    order = []
    for r in runs:
        modes, U = r[6], r[7]
        for i in np.argsort(-U.max(axis=1)):
            if modes[i] not in order and U[i].max() > 1e-6 * U.max():
                order.append(modes[i])
    color = {md: LEVEL_COLORS[k] for k, md in enumerate(order[:len(LEVEL_COLORS)])}
    ncol = len(runs)
    fig, axes = plt.subplots(3, ncol, figsize=(4.6 * ncol + 0.6, 8.2), sharex="col", squeeze=False,
                             gridspec_kw=dict(height_ratios=[1, 1, 3], hspace=0.06, wspace=0.28))
    # common axes across the columns: time from a decade before the earliest growth (the cloud first reaches
    # 1e-9 of its peak) to the end; spin and BH-mass ranges spanning all runs
    t_lo, t_hi, s_rng, m_rng = np.inf, 0.0, [np.inf, -np.inf], [np.inf, -np.inf]
    for (leaf, lab, a0, T, S, M, modes, U) in runs:
        m = T > 0
        tot = a0 * U[:, m].sum(axis=0)
        k = np.argmax(tot > 1e-9 * tot.max())
        t_lo = min(t_lo, T[m][k] / 10); t_hi = max(t_hi, T[m][-1])
        s_rng = [min(s_rng[0], S[m].min()), max(s_rng[1], S[m].max())]
        mm = M[m] / M[m][0]
        m_rng = [min(m_rng[0], mm.min()), max(m_rng[1], mm.max())]
    pad = lambda r: (r[0] - 0.05 * (r[1] - r[0]), r[1] + 0.05 * (r[1] - r[0]))
    for c, (leaf, lab, a0, T, S, M, modes, U) in enumerate(runs):
        m = T > 0
        T, S, M, U = T[m], S[m], M[m], U[:, m]
        axS, axM, axC = axes[0, c], axes[1, c], axes[2, c]
        axS.plot(T, S, color="k", lw=1.6)
        axM.plot(T, M / M[0], color="k", lw=1.6)
        frac = a0 * U
        tot = frac.sum(axis=0)
        peak = tot.max()
        floor = max(peak * 1e-9, 1e-14)
        for i in np.argsort(U.max(axis=1)):          # draw small levels first
            col = color.get(modes[i], OTHER)
            if frac[i].max() < floor:
                continue
            axC.plot(T, np.where(frac[i] > floor * 1e-2, frac[i], np.nan), color=col,
                     lw=1.6 if col != OTHER else 0.8, zorder=3 if col != OTHER else 2)
        axC.plot(T, tot, color="k", lw=3.4, zorder=1)          # total underneath: shows as a black edge where one level dominates
        placed = []                                     # (log10 t, log10 y) of labels already drawn
        for i in np.argsort(-U.max(axis=1))[:3]:      # direct labels on the 3 dominant levels, nudged apart
            k = np.argmax(frac[i])
            lx, ly = np.log10(T[k]), np.log10(frac[i][k])
            while any(abs(lx - px) < 1.0 and abs(ly - py) < 0.35 for px, py in placed):
                ly -= 0.4
            placed.append((lx, ly))
            n_, l_, m_ = modes[i]
            right = lx > np.log10(t_hi) - 0.12 * (np.log10(t_hi) - np.log10(t_lo))     # near the right edge: label to the left
            axC.annotate(rf"${n_}{l_}{m_}$", (10**lx, 10**ly), xytext=(-5 if right else 5, 3), textcoords="offset points",
                         ha="right" if right else "left", fontsize=12, color=color.get(modes[i], "#52514e"))
        for ax in (axS, axM, axC):
            ax.set_xscale("log")
            ax.set_xlim(t_lo, t_hi)
        axS.set_ylim(*pad(s_rng))
        axM.set_ylim(*pad(m_rng))
        axC.set_yscale("log")
        axC.set_ylim(floor, peak * 3)
        axS.set_title(lab, fontsize=LABEL_FS)
        axC.set_xlabel(r"$t$ [yrs]", fontsize=LABEL_FS)
        if c == 0:
            axS.set_ylabel(r"$a$", fontsize=LABEL_FS)
            axM.set_ylabel(r"$M_{\rm BH}/M_0$", fontsize=LABEL_FS)
            axC.set_ylabel(r"$M_{\rm cloud}/M_0$", fontsize=LABEL_FS)
        axM.ticklabel_format(axis="y", useOffset=False)
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(os.path.splitext(out)[0] + ".png", dpi=160, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
