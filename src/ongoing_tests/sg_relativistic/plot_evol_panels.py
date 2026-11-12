"""Level occupations (top) and BH spin (bottom) of run_Nlevels.jl outputs, one panel per (alpha, f_a).

Generalisation of sg_vs_si/plot_grid_runs.py: several line sets (e.g. SI only / SI + exact gravity / SI + Newtonian
fixed-BH gravity), and either an alpha x f_a grid or a single f_a with the alphas wrapped into columns.
Run directories are <root>/a<alpha>_f<fa>_<mode> as written by run_one.sh.

usage:
  python3 plot_evol_panels.py out.pdf --alphas 0.05,0.1,0.2 --fas 1e16,1e17,1e18 \
      --set "SI only|runs/grid|si|-" --set "SI + gravity|runs/grid_rel|sg|--"
  python3 plot_evol_panels.py out.pdf --alphas 0.05,0.1,...,0.8 --fas 1e30 --ncols 3 --title '$f_a\\to\\infty$' ...
"""
import argparse
import glob
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

STATES = [(2, 1, 1), (3, 2, 2), (4, 3, 3), (5, 4, 4), (7, 6, 6)]
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
YMIN = 1e-10
TMIN, TMAX = 1e-4, 1e10
STYLES = {"-": "-", "--": (0, (3.5, 1.8)), ":": (0, (1, 1.5)), "-.": (0, (5, 1.5, 1, 1.5))}
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 8,
                     "axes.edgecolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.labelcolor": INK, "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.labelsize": 7, "ytick.labelsize": 7})

ap = argparse.ArgumentParser()
ap.add_argument("out")
ap.add_argument("--alphas", required=True)
ap.add_argument("--fas", required=True)
ap.add_argument("--set", action="append", required=True, help="label|root|mode|style")
ap.add_argument("--ncols", type=int, default=0, help="wrap the alphas into this many columns (single f_a)")
ap.add_argument("--title", default="")
ap.add_argument("--tmin", type=float, default=TMIN)
ap.add_argument("--width", type=float, default=7.0)
ap.add_argument("--legend_ncol", type=int, default=0)
ap.add_argument("--panel_label", action="store_true", help="write alpha and f_a inside each panel")
args = ap.parse_args()
alphas, fas = args.alphas.split(","), args.fas.split(",")
sets = [dict(zip(("label", "root", "mode", "style"), s.split("|"))) for s in args.set]


def load(d):
    f = lambda p: glob.glob(os.path.join(d, p + "_*.dat"))[0]
    t = np.loadtxt(f("Time")); sp = np.loadtxt(f("Spin")); md = np.loadtxt(f("Modes")); st = np.loadtxt(f("States"))
    st = st if st.shape[0] == md.shape[0] else st.T
    idx = {(int(m[0]), int(m[1]), int(m[2])): i for i, m in enumerate(md)}
    return t, sp, st, idx


def fa_label(fa):
    if fa >= 1e25:
        return r"$f_a\to\infty$"
    e = int(np.floor(np.log10(fa) + 1e-9)); m = fa / 10 ** e
    return r"$f_a=10^{%d}$ GeV" % e if abs(m - 1) < 1e-6 else r"$f_a=%.3g\times10^{%d}$ GeV" % (m, e)


if args.ncols:
    assert len(fas) == 1
    panels = [(a, fas[0]) for a in alphas]
    ncols = args.ncols
else:
    panels = [(a, fa) for a in alphas for fa in fas]
    ncols = len(fas)
nrows = (len(panels) + ncols - 1) // ncols
W = args.width
nleg_rows = 1 if not args.legend_ncol else 2
col_titles = not (args.ncols or args.panel_label)
top_in = 0.22 + 0.14 * nleg_rows + (0.14 if col_titles else 0.0) + (0.2 if args.title else 0.0)
fig = plt.figure(figsize=(W, 2.55 * nrows + top_in))
gs = GridSpec(2 * nrows, ncols, figure=fig, height_ratios=[2.3, 1.0] * nrows, hspace=0.0, wspace=0.07,
              left=0.6 / W, right=1 - 0.1 / W, top=1 - top_in / fig.get_figheight(), bottom=0.42 / fig.get_figheight())
seen = set()
summary = []
for k, (a, fa) in enumerate(panels):
    i, j = divmod(k, ncols)
    ax = fig.add_subplot(gs[2 * i, j]); axs = fig.add_subplot(gs[2 * i + 1, j], sharex=ax)
    failed = []
    for S in sets:
        d = os.path.join(S["root"], f"a{a}_f{fa}_{S['mode']}")
        if not glob.glob(os.path.join(d, "Time_*.dat")):
            continue
        t, sp, st, idx = load(d)
        if t[-1] < 0.99 * TMAX:
            failed.append((S["label"], t[-1]))
        ls = STYLES.get(S["style"], S["style"])
        for s, c in zip(STATES, COLORS):
            y = st[idx[s]]
            if y.max() > YMIN:
                seen.add(s)
            ax.plot(t, np.maximum(y, 1e-30), color=c, ls=ls, lw=1.1)
        axs.plot(t, sp, color=INK, ls=ls, lw=1.0)
        summary.append((a, fa, S["label"], {s: st[idx[s]].max() for s in STATES},
                        {T: np.interp(T, t, sp) for T in (1e6, 1e8, 1e9, 1e10)}))
    if failed:
        ax.text(0.5, 0.5, "solver stopped at\n" + ", ".join(r"$t=%.2g$ yr (%s)" % (te, m) for m, te in failed),
                transform=ax.transAxes, ha="center", va="center", fontsize=7, color=INK2)
    for x in (ax, axs):
        x.set_xscale("log"); x.set_xticks([1e-3, 1e0, 1e3, 1e6, 1e9])
        x.set_xlim(args.tmin, TMAX)
        x.grid(True, color=GRID, lw=0.5); x.set_axisbelow(True)
    ax.set_yscale("log"); ax.set_ylim(YMIN, 2); ax.set_yticks([1e-9, 1e-6, 1e-3, 1])
    axs.set_ylim(0, 1); axs.set_yticks([0.2, 0.5, 0.8])
    plt.setp(ax.get_xticklabels(), visible=False)
    last_row = (k + ncols >= len(panels))
    if not last_row:
        plt.setp(axs.get_xticklabels(), visible=False)
    else:
        axs.set_xlabel(r"$t$ [yr]")
    if j == 0:
        ax.set_ylabel(r"$N/GM_0^2$"); axs.set_ylabel(r"$\tilde a$")
    else:
        plt.setp(ax.get_yticklabels(), visible=False); plt.setp(axs.get_yticklabels(), visible=False)
    if args.panel_label:
        ax.text(0.03, 0.05, r"$\alpha=%s$" % a + "\n" + fa_label(float(fa)), transform=ax.transAxes, ha="left", va="bottom",
                fontsize=7.5, color=INK, bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=GRID, lw=0.6))
    elif args.ncols:
        ax.text(0.97, 0.93, r"$\alpha=%s$" % a, transform=ax.transAxes, ha="right", va="top", fontsize=8.5,
                color=INK, bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=GRID, lw=0.6))
    else:
        if i == 0:
            ax.set_title(fa_label(float(fa)), color=INK, fontsize=8.5)
        if j == ncols - 1:
            ax.text(0.97, 0.93, r"$\alpha=%s$" % a, transform=ax.transAxes, ha="right", va="top", fontsize=8.5,
                    color=INK, bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=GRID, lw=0.6))

handles = [Line2D([], [], color=c, lw=1.6, label=r"$|%d%d%d\rangle$" % s) for s, c in zip(STATES, COLORS) if s in seen]
set_handles = [Line2D([], [], color=INK2, lw=1.2, ls=STYLES.get(S["style"], S["style"]), label=S["label"]) for S in sets]
if args.legend_ncol and max(len(handles), len(set_handles)) <= args.legend_ncol:
    # two rows (states / line styles); the legend fills column-major, so interleave
    blank = lambda: Line2D([], [], alpha=0, label=" ")
    row1 = handles + [blank() for _ in range(args.legend_ncol - len(handles))]
    row2 = set_handles + [blank() for _ in range(args.legend_ncol - len(set_handles))]
    handles = [h for pair in zip(row1, row2) for h in pair]
else:
    handles += set_handles
narrow = W < 4
fig.legend(handles=handles, loc="upper center", ncol=args.legend_ncol or len(handles), frameon=False,
           fontsize=7 if narrow else 8, bbox_to_anchor=(0.5 if narrow else 0.55, 1.0), handlelength=1.8 if narrow else 2.2,
           columnspacing=0.8 if narrow else 1.0, labelspacing=0.3, handletextpad=0.5)
if args.title:
    fig.text(0.55, 1 - (0.12 + 0.14 * nleg_rows) / fig.get_figheight(), args.title, ha="center", va="top", fontsize=9, color=INK)
fig.savefig(args.out)
fig.savefig(os.path.splitext(args.out)[0] + ".png", dpi=170)

print("alpha    f_a   set                 a(1e6)  a(1e8)  a(1e9)  a(1e10)  " + "  ".join("max%d%d%d" % s for s in STATES))
for a, fa, lab, mx, aT in summary:
    print(f"{a:5s} {fa:6s} {lab[:18]:18s} " + "  ".join("%.3f " % aT[T] for T in (1e6, 1e8, 1e9, 1e10)) + "  "
          + "  ".join("%.1e" % mx[s] for s in STATES))
