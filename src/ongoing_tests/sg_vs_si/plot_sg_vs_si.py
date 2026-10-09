"""Plot Gamma_SG / Gamma_SI for the channels computed by sg_vs_si.jl.

Gamma_SG / Gamma_SI = C(alpha) (f_a / M_pl)^4, with C read from out/<channel>/REL_*.dat (Kerr
eigenfunctions) and out/<channel>/NR_*.dat (hydrogenic eigenfunctions, same Green's function).
Writes one two-panel figure + summary table per channel, and an overview of the crossover f_a.
"""
import glob
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

HERE = os.path.dirname(os.path.abspath(__file__))
M_PL = 1.22e19  # GeV, as in Core/constants.jl
CHANNELS = {
    "211_211_322_BH": r"$211\times211\to322\times{\rm BH}$",
    "322_322_544_BH": r"$322\times322\to544\times{\rm BH}$",
    "322_322_211_Inf": r"$322\times322\to211\times\infty$",
    "544_544_322_Inf": r"$544\times544\to322\times\infty$",
}
# columns: alpha a ReE1 ImE1 ReE3 ImE3 ReE rate_SI rate_SI_table C C_low phase xfrac_SI xfrac_SG
COL_C, COL_RATE, COL_TAB = 9, 7, 8


def load(ch, pattern):
    files = sorted(glob.glob(os.path.join(HERE, "out", ch, pattern)))
    if not files:
        return None
    d = np.array([np.loadtxt(f) for f in files])
    d = d[np.isfinite(d[:, COL_C]) & (d[:, COL_C] > 0)]
    return d[np.argsort(d[:, 0])]


INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
BLUES = ["#5598e7", "#2a78d6", "#1c5cab", "#0d366b"]          # sequential ramp, ordered by f_a
CATEG = ["#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"]          # categorical slots (channel identity)
CMAP = LinearSegmentedColormap.from_list(
    "bluered", ["#104281", "#3987e5", "#b7d3f6", "#f0efec", "#f6c3c0", "#e34948", "#8f1f1f"])
plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "cm", "font.size": 10.5,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.linewidth": 0.8, "xtick.direction": "in", "ytick.direction": "in",
})


def channel_figure(ch, label, rel, nr, spin):
    alpha, C = rel[:, 0], rel[:, COL_C]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)

    for fa, col in zip([1e15, 1e16, 1e17, 1e18], BLUES):
        ratio = C * (fa / M_PL) ** 4
        ax1.plot(alpha, ratio, color=col, lw=2, solid_capstyle="round")
        ax1.text(alpha[-1] + 0.015 * (alpha[-1] - alpha[0]), ratio[-1], rf"$f_a=10^{{{int(np.log10(fa))}}}$ GeV",
                 color=INK, va="center", ha="left", fontsize=9.5)
    ax1.axhline(1.0, color=INK2, lw=1, ls=(0, (4, 3)))
    ax1.text(0.02, 1.4, r"$\Gamma_{\rm SG}=\Gamma_{\rm SI}$", color=INK2, fontsize=9.5, va="bottom",
             transform=ax1.get_yaxis_transform())
    ax1.set_yscale("log")
    ax1.set_xlim(alpha[0], alpha[-1] + 0.34 * (alpha[-1] - alpha[0]))
    ax1.set_xlabel(r"$\alpha = G M \mu$")
    ax1.set_ylabel(r"$\Gamma_{\rm SG}\,/\,\Gamma_{\rm SI}$")
    ax1.set_title(rf"(a)  {label},  $\tilde a=0.95$", color=INK, fontsize=10.5, loc="left")
    ax1.grid(True, which="major", color=GRID, lw=0.6)
    ax1.set_axisbelow(True)
    ax1.set_xticks([t for t in ax1.get_xticks() if alpha[0] <= t <= alpha[-1]])

    ag = np.linspace(alpha[0], alpha[-1], 300)
    lfg = np.linspace(14, 19, 300)
    Z = np.interp(ag, alpha, np.log10(C))[None, :] + 4 * (lfg[:, None] - np.log10(M_PL))
    lim = 12
    cf = ax2.contourf(ag, lfg, np.clip(Z, -lim, lim), levels=np.arange(-lim, lim + 0.01, 2), cmap=CMAP,
                      norm=TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim), extend="both")
    ax2.contour(ag, lfg, Z, levels=[0], colors=INK, linewidths=2)
    if nr is not None:
        Znr = np.interp(ag, nr[:, 0], np.log10(nr[:, COL_C]))[None, :] + 4 * (lfg[:, None] - np.log10(M_PL))
        ax2.contour(ag, lfg, Znr, levels=[0], colors=INK, linewidths=1.2, linestyles=[(0, (4, 3))])
    if spin is not None:
        ax2.plot(spin[:, 0], np.log10(M_PL * spin[:, COL_C] ** -0.25), "o", ms=8, mfc="white", mec=INK, mew=1.4,
                 label=r"$\tilde a=0.5$ (check)")
    ax2.plot([], [], color=INK, lw=2, label=r"$\Gamma_{\rm SG}=\Gamma_{\rm SI}$ (Kerr eigenfunctions)")
    ax2.plot([], [], color=INK, lw=1.2, ls=(0, (4, 3)), label="same, hydrogenic eigenfunctions")
    ax2.legend(loc="lower right", fontsize=8.5, frameon=True, facecolor="white", edgecolor=GRID)
    ax2.text(0.06, 0.88, "self-gravity dominates", color=INK, fontsize=10, transform=ax2.transAxes)
    ax2.text(0.02, 0.05, "self-interactions dominate", color="white", fontsize=9.5, transform=ax2.transAxes)
    ax2.set_xlabel(r"$\alpha = G M \mu$")
    ax2.set_ylabel(r"$\log_{10}\,(f_a\,/\,{\rm GeV})$")
    ax2.set_title(r"(b)  which effect dominates", color=INK, fontsize=10.5, loc="left")
    cb = fig.colorbar(cf, ax=ax2, pad=0.02)
    cb.set_label(r"$\log_{10}(\Gamma_{\rm SG}/\Gamma_{\rm SI})$", color=INK)
    cb.outline.set_edgecolor(INK2)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(HERE, f"sg_vs_si_{ch}.{ext}"), dpi=200)
    plt.close(fig)


def write_summary(ch, rel, nr):
    rows = []
    for r in rel:
        Cnr = np.interp(r[0], nr[:, 0], nr[:, COL_C]) if nr is not None else np.nan
        rows.append([r[0], r[COL_C], Cnr, M_PL * r[COL_C] ** -0.25, r[COL_RATE], r[COL_TAB], r[COL_C] * r[COL_RATE]])
    hdr = ("alpha  C_kerr  C_hydrogenic  f_a_equal[GeV]  GammaSI(f_a=M_pl, this calc)  "
           "GammaSI(table LvrHc)  GammaSG (f_a-independent)   [rates in the units of rate_sve tables]")
    np.savetxt(os.path.join(HERE, f"sg_vs_si_{ch}_summary.dat"), np.array(rows), header=hdr, fmt="%.6e")


fig, ax = plt.subplots(figsize=(6.4, 4.4), constrained_layout=True)
for (ch, label), col in zip(CHANNELS.items(), CATEG):
    rel, nr, spin = load(ch, "REL_[0-9][0-9].dat"), load(ch, "NR_[0-9][0-9].dat"), load(ch, "REL_a0.5_*.dat")
    if rel is None:
        continue
    channel_figure(ch, label, rel, nr, spin)
    write_summary(ch, rel, nr)
    ax.plot(rel[:, 0], M_PL * rel[:, COL_C] ** -0.25, color=col, lw=2, label=label)
    if nr is not None:
        ax.plot(nr[:, 0], M_PL * nr[:, COL_C] ** -0.25, color=col, lw=1, ls=(0, (4, 3)))
ax.plot([], [], color=INK2, lw=1, ls=(0, (4, 3)), label="hydrogenic eigenfunctions")
ax.set_yscale("log")
ax.set_xlabel(r"$\alpha = G M \mu$")
ax.set_ylabel(r"$f_a$ [GeV] at which $\Gamma_{\rm SG}=\Gamma_{\rm SI}$")
ax.set_title(r"Self-gravity dominates above each curve  ($\tilde a=0.95$)", color=INK, fontsize=10.5, loc="left")
ax.grid(True, which="major", color=GRID, lw=0.6)
ax.set_axisbelow(True)
ax.legend(fontsize=9, frameon=True, facecolor="white", edgecolor=GRID, loc="lower right")
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(HERE, f"sg_vs_si_overview.{ext}"), dpi=200)
print("saved")
