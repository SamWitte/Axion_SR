#!/usr/bin/env python3
"""Relativistic (Teukolsky, rate_sve/gw_kerr/ann_*.dat) vs non-relativistic annihilation rates, log y.
  python3 plot_gw_kerr_annihilation.py OUT.png [stem ...]     stem e.g. ann_2-1-1_2-1-1
Per channel: Γ_∞ (Kerr, flux to infinity, solver units) against α at four spins (one sequential hue,
light = low spin; filled markers where both levels grow, hollow where a level is absorbed), and the
non-relativistic rate C α^p used by gw_model=:nonrel (dashed, spin independent; l = 1 pairs carry the
BH-potential α^14 override). Only rows with status ok are drawn.
"""
import os
import re
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
TAB = os.path.join(HERE, "..", "rate_sve", "gw_kerr")
SPIN_TARGETS = [0.01, 0.516, 0.82, 0.998]
SPIN_COLORS = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]   # sequential blue 250 / 400 / 550 / 700
NR_COLOR = "#52514e"
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
DEFAULT = ["ann_2-1-1_2-1-1", "ann_2-1-1_3-2-2", "ann_3-2-2_3-2-2", "ann_4-3-3_4-3-3", "ann_5-4-4_5-4-4", "ann_6-5-5_6-5-5"]


def read(stem):
    path = os.path.join(TAB, stem + ".dat")
    C = p = None
    rows = []
    for ln in open(path):
        if ln.startswith("#"):
            m = re.search(r"\(C, p\) = \(([^,]+), (\d+)\)", ln)
            if m:
                C, p = float(m.group(1)), int(m.group(2))
            continue
        f = ln.split()
        if len(f) < 14 or f[13] != "ok":
            continue
        rows.append(dict(alpha=float(f[0]), a=float(f[1]), gk=float(f[2]), gh=float(f[10]),
                         ima=float(f[8]), imb=float(f[9])))
    return C, p, rows


def label(stem):
    _, a, b = stem.split("_")
    return f"{a.replace('-', '')} × {b.replace('-', '')}"


def main():
    out = sys.argv[1]
    stems = sys.argv[2:] or [s for s in DEFAULT if os.path.isfile(os.path.join(TAB, s + ".dat"))]
    ncol = 3
    nrow = (len(stems) + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.5 * nrow), squeeze=False, facecolor=SURF)
    for k, s in enumerate(stems):
        ax = axes[k // ncol, k % ncol]
        C, p, rows = read(s)
        spins = sorted({r["a"] for r in rows})
        amin, amax = np.inf, 0.0
        for target, col in zip(SPIN_TARGETS, SPIN_COLORS):
            if not spins:
                break
            a = min(spins, key=lambda x: abs(x - target))
            if abs(a - target) > 0.06:
                continue
            sel = sorted((r for r in rows if r["a"] == a and r["gk"] > 0), key=lambda r: r["alpha"])
            if not sel:
                continue
            al = np.array([r["alpha"] for r in sel]); g = np.array([r["gk"] for r in sel])
            grow = np.array([r["ima"] > 0 and r["imb"] > 0 for r in sel])
            ax.plot(al, g, "-", color=col, lw=2)
            ax.plot(al[grow], g[grow], "o", color=col, ms=4)
            ax.plot(al[~grow], g[~grow], "o", mfc=SURF, mec=col, mew=1.2, ms=4)
            amin, amax = min(amin, al.min()), max(amax, al.max())
        if C is not None and amax > 0:
            x = np.linspace(amin, amax, 300)
            ax.plot(x, C * x ** p, "--", color=NR_COLOR, lw=1.6)
            ax.text(0.03, 0.95, f"NR: {C:.3g} α^{p}", transform=ax.transAxes, va="top", fontsize=8, color=INK2)
        ax.set_yscale("log")
        ax.set_facecolor(SURF)
        ax.grid(True, color=GRID, lw=0.6)
        ax.tick_params(colors=INK2, labelsize=8)
        for sp in ax.spines.values():
            sp.set_color(GRID)
        ax.set_title(label(s), color=INK, fontsize=10)
        ax.set_xlabel("α", color=INK2, fontsize=9)
        if k % ncol == 0:
            ax.set_ylabel("Γ  (solver units, per u_a u_b)", color=INK2, fontsize=9)
    for k in range(len(stems), nrow * ncol):
        axes[k // ncol, k % ncol].axis("off")
    handles = [Line2D([], [], color=c, lw=2, marker="o", ms=4, label=f"Kerr, a ≈ {t:g}") for t, c in zip(SPIN_TARGETS, SPIN_COLORS)]
    handles += [Line2D([], [], color=NR_COLOR, lw=1.6, ls="--", label="non-relativistic"),
                Line2D([], [], color=INK2, lw=0, marker="o", mfc=SURF, mec=INK2, ms=4, label="a level is absorbed")]
    fig.legend(handles=handles, loc="lower center", ncol=len(handles), frameon=False, fontsize=8.5, labelcolor=INK2)
    fig.suptitle("GW annihilation rates: Kerr (flux to infinity) vs non-relativistic", color=INK, fontsize=11)
    fig.tight_layout(rect=(0, 0.05, 1, 0.97))
    fig.savefig(out, dpi=150, facecolor=SURF)
    print("wrote", out)


if __name__ == "__main__":
    main()
