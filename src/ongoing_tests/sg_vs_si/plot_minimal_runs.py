"""Compare two run_Nlevels.jl outputs (SI only vs SI + self-gravity) for the minimal 5-state channel list.

usage: python3 plot_minimal_runs.py <dir_SI> <dir_SG> <out.png>
"""
import glob
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

STATES = [(2, 1, 1), (3, 2, 2), (4, 3, 3), (5, 4, 4), (7, 6, 6)]
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 10,
                     "axes.edgecolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "xtick.direction": "in", "ytick.direction": "in"})


def load(d):
    f = lambda p: glob.glob(os.path.join(d, p + "_*.dat"))[0]
    t = np.loadtxt(f("Time")); sp = np.loadtxt(f("Spin")); md = np.loadtxt(f("Modes")); st = np.loadtxt(f("States"))
    st = st if st.shape[0] == md.shape[0] else st.T
    idx = {(int(m[0]), int(m[1]), int(m[2])): i for i, m in enumerate(md)}
    return t, sp, st, idx


runs = [("SI only", load(sys.argv[1]), "-"), ("SI + self-gravity", load(sys.argv[2]), (0, (4, 2)))]
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.5, 7), sharex=True, constrained_layout=True,
                               gridspec_kw={"height_ratios": [3, 1.2]})
for label, (t, sp, st, idx), ls in runs:
    for s, c in zip(STATES, COLORS):
        y = np.maximum(st[idx[s]], 1e-30)
        ax1.plot(t, y, color=c, ls=ls, lw=1.8)
    ax2.plot(t, sp, color=INK, ls=ls, lw=1.6, label=label)
for s, c in zip(STATES, COLORS):
    ax1.plot([], [], color=c, lw=2, label=r"$|%d%d%d\rangle$" % s)
ax1.plot([], [], color=INK2, lw=1.6, ls="-", label="SI only")
ax1.plot([], [], color=INK2, lw=1.6, ls=(0, (4, 2)), label="SI + self-gravity")
ax1.set_xscale("log"); ax1.set_yscale("log"); ax1.set_ylim(1e-12, 2)
ax1.set_ylabel(r"occupation  $N/(GM_0^2)$")
ax1.legend(ncol=2, fontsize=8.5, frameon=True, facecolor="white", edgecolor=GRID, loc="lower left")
ax1.grid(True, color=GRID, lw=0.6); ax1.set_axisbelow(True)
ax2.set_xscale("log"); ax2.set_xlabel("t  [yr]"); ax2.set_ylabel(r"$\tilde a$")
ax2.grid(True, color=GRID, lw=0.6); ax2.legend(fontsize=8.5, frameon=False)
ax1.set_title(sys.argv[4] if len(sys.argv) > 4 else "minimal 5-state network", color=INK, fontsize=10.5, loc="left")
fig.savefig(sys.argv[3], dpi=170)
for label, (t, sp, st, idx), _ in runs:
    print(label, " final spin %.4f" % sp[-1], "  ".join("%d%d%d: max %.2e final %.2e" % (s + (st[idx[s]].max(), st[idx[s]][-1])) for s in STATES))
