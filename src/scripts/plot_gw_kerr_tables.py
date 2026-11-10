#!/usr/bin/env python3
"""Inspect the relativistic GW correction tables (rate_sve/gw_kerr/*.dat, gw_kerr_build.jl).
  python3 plot_gw_kerr_tables.py OUT.png [channel ...]      channel = file stem, e.g. tr_3-2-2_2-1-1
Prints per-table row status and the ratio / horizon fraction ranges, and plots, per channel,
Γ_Kerr/Γ_NR (top) and Γ_H/Γ_∞ (bottom) against α for four spins (one sequential hue, light = low spin).
"""
import glob
import os
import sys
import collections
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker

HERE = os.path.dirname(os.path.abspath(__file__))
TAB = os.path.join(HERE, "..", "rate_sve", "gw_kerr")
SPIN_TARGETS = [0.01, 0.516, 0.82, 0.998]
SPIN_COLORS = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]   # sequential blue 250 / 400 / 550 / 700
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"


def read(path):
    rows = []
    for ln in open(path):
        if ln.startswith("#") or not ln.strip():
            continue
        p = ln.split()
        if len(p) < 14:
            continue
        rows.append(dict(alpha=float(p[0]), a=float(p[1]), gk=float(p[2]), ratio=float(p[4]),
                         gh=float(p[10]), status=p[13], chk=float(p[14]) if len(p) > 14 else np.nan))
    return rows


def label(stem):
    kind, a, b = stem.split("_")
    a, b = a.replace("-", ""), b.replace("-", "")
    return f"{a} → {b}" if kind == "tr" else f"{a} × {b}"


def main():
    out = sys.argv[1]
    files = sorted(glob.glob(os.path.join(TAB, "*.dat")))
    stems = [os.path.basename(f)[:-4] for f in files]
    print(f"{len(files)} tables in {os.path.abspath(TAB)}")
    tot = collections.Counter()
    for f, s in zip(files, stems):
        rows = read(f)
        c = collections.Counter(r["status"].split(":")[0] if not r["status"].startswith("fail") else "fail" for r in rows)
        tot.update(c)
        ok = [r for r in rows if r["status"] == "ok"]
        if ok:
            rat = np.array([r["ratio"] for r in ok]); hf = np.array([r["gh"] / r["gk"] for r in ok if np.isfinite(r["gh"])])
            print(f"  {s:28s} rows {len(rows):4d}  ok {c['ok']:4d} noisy {c['noisy']:3d} fail {c['fail']:3d}   "
                  f"ratio [{rat.min():.3g}, {rat.max():.3g}]   Γ_H/Γ_∞ [{hf.min() if hf.size else np.nan:.2g}, {hf.max() if hf.size else np.nan:.2g}]")
        else:
            print(f"  {s:28s} rows {len(rows):4d}  no ok row  {dict(c)}")
    print("total:", dict(tot))
    pick = sys.argv[2:] or [s for s in ("tr_3-2-2_2-1-1", "tr_4-3-3_3-2-2", "tr_4-3-3_2-1-1", "tr_5-4-4_4-3-3",
                                         "ann_2-1-1_2-1-1", "ann_3-2-2_3-2-2") if s in stems]
    if not pick:
        return
    fig, axes = plt.subplots(2, len(pick), figsize=(3.1 * len(pick), 5.6), squeeze=False, facecolor=SURF)
    for j, s in enumerate(pick):
        rows = [r for r in read(os.path.join(TAB, s + ".dat")) if r["status"] == "ok"]
        spins = sorted({r["a"] for r in rows})
        for target, col in zip(SPIN_TARGETS, SPIN_COLORS):
            if not spins:
                break
            a = min(spins, key=lambda x: abs(x - target))
            if abs(a - target) > 0.06:
                continue
            sel = sorted((r for r in rows if r["a"] == a), key=lambda r: r["alpha"])
            al = np.array([r["alpha"] for r in sel])
            axes[0, j].plot(al, [r["ratio"] for r in sel], "-o", color=col, lw=2, ms=3.5, label=f"a = {a:.2f}")
            hf = np.array([r["gh"] / r["gk"] for r in sel])
            axes[1, j].plot(al, hf, "-o", color=col, lw=2, ms=3.5)
        for i in (0, 1):
            ax = axes[i, j]
            ax.set_facecolor(SURF); ax.set_xscale("log")
            ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
            ax.grid(True, color=GRID, lw=0.6); ax.tick_params(colors=INK2, labelsize=8)
            for sp in ax.spines.values():
                sp.set_color(GRID)
        axes[0, j].set_yscale("log")
        axes[1, j].set_yscale("symlog", linthresh=1e-4)
        axes[0, j].set_title(label(s), color=INK, fontsize=10)
        axes[1, j].set_xlabel("α", color=INK2, fontsize=9)
        axes[0, j].axhline(1, color=INK2, lw=0.8, ls=":")
        axes[1, j].axhline(0, color=INK2, lw=0.8, ls=":")
    axes[0, 0].set_ylabel("Γ_Kerr / Γ_NR", color=INK2, fontsize=9)
    axes[1, 0].set_ylabel("Γ_H / Γ_∞  (symlog)", color=INK2, fontsize=9)
    axes[0, 0].legend(fontsize=7.5, frameon=False, labelcolor=INK2)
    fig.suptitle("Relativistic GW corrections (Teukolsky tables)", color=INK, fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=150, facecolor=SURF)
    print("wrote", out)


if __name__ == "__main__":
    main()
