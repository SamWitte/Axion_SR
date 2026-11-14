#!/usr/bin/env python
"""plot_sg_states.py -- spin (top, small) and mode occupations epsilon(t) (bottom, large) for two variants of the same
(BH mass, alpha) run, in the style of arXiv:2412.03655 Fig. 5: variant A solid, variant B dashed, one colour per |nlm>.

  python plot_sg_states.py --alpha 0.2 --solid SG_only --dashed noSG_noSI [--eps-min 1e-20] [--tmin T] [--out f.png]
"""
import argparse, glob, pathlib, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import AutoMinorLocator, LogLocator

ap = argparse.ArgumentParser()
ap.add_argument("--root", default="output"); ap.add_argument("--bh", default="10"); ap.add_argument("--alpha", required=True)
ap.add_argument("--solid", required=True); ap.add_argument("--dashed", required=True)
ap.add_argument("--eps-min", type=float, default=1e-10); ap.add_argument("--tmin", type=float, default=None); ap.add_argument("--tmax", type=float, default=5e9)
ap.add_argument("--solid-dir", default=None); ap.add_argument("--dashed-dir", default=None)   # explicit run dirs (override root/BH_/variant)
ap.add_argument("--quantity", default="eps", choices=["eps", "mc"], help="eps: occupation eps_x = N/(G M_0^2); mc: cloud mass M_c/M_0 = alpha_0 eps_x")
ap.add_argument("--out", default=None); ap.add_argument("--max-modes", type=int, default=20); ap.add_argument("--label-pos", default="auto", choices=["auto", "lower left", "upper left"])
a = ap.parse_args()
# eps_x = N_x/(G M_0^2) with M_0 the initial BH mass (solver normalisation), so M_c = N m_a = alpha_0 eps_x M_0
sc = float(a.alpha) if a.quantity == "mc" else 1.0
cw = 0.115 if a.quantity == "eps" else 0.10          # label column width (axes fraction)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "TeX Gyre Heros", "Nimbus Sans", "Arial", "Liberation Sans", "DejaVu Sans"],
    "mathtext.fontset": "dejavusans",
    "font.size": 15, "axes.labelsize": 18, "xtick.labelsize": 14, "ytick.labelsize": 14,
    "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
    "xtick.major.size": 6, "ytick.major.size": 6, "xtick.minor.size": 3.5, "ytick.minor.size": 3.5,
    "xtick.major.width": 1.0, "ytick.major.width": 1.0, "axes.linewidth": 1.2, "axes.grid": False,
})
PRETTY = {"noSG_noSI": "no SG, no SI", "SG_only": "self-gravity", "SG_SI_fa_1e16": r"SG + SI, $f_a=10^{16}\,$GeV",
          "SG_SI_fa_5e16": r"SG + SI, $f_a=5\times10^{16}\,$GeV", "SG_SI_fa_1e17": r"SG + SI, $f_a=10^{17}\,$GeV"}
COLORS = {"211": "#8b1a1a", "322": "#1f3a68", "411": "#2a9d8f", "422": "#c9a227", "433": "#7fb3d5", "533": "#6a3d9a",
          "544": "#e07b00", "633": "#d6336c", "644": "#4b8b3b", "733": "#00a2c7", "744": "#7a4a1a"}
FALLBACK = ["#5c5c5c", "#808000", "#000000", "#7fc41c", "#fa8072", "#2f4f4f", "#ff69b4", "#bdb76b", "#8fbc8f", "#4b0082"]

def one(path, tag): return glob.glob(str(path / f"{tag}_*.dat"))[0]

def load(variant, d=None):
    d = pathlib.Path(d) if d else pathlib.Path(a.root) / f"BH_{a.bh}" / variant / f"alpha_{a.alpha}"
    t = np.loadtxt(one(d, "Time")); spin = np.loadtxt(one(d, "Spin")); mass = np.loadtxt(one(d, "MassBH"))
    modes = np.loadtxt(one(d, "Modes"))[:, :3].astype(int); keep = {}
    with open(one(d, "States")) as f:                       # stream: one tab-separated row per mode
        for i, line in enumerate(f):
            row = np.fromstring(line, sep="\t")
            if row.max() > a.eps_min: keep[tuple(int(x) for x in modes[i])] = row
    n = min([len(t), len(spin), len(mass)] + [len(v) for v in keep.values()])
    return dict(t=t[:n], spin=spin[:n], mass=mass[:n], eps={k: v[:n] for k, v in keep.items()}, name=variant)

A = load(a.solid, a.solid_dir); B = None if a.dashed == "none" else load(a.dashed, a.dashed_dir)   # --dashed none: one run
RUNS = [A] + ([B] if B else [])
if a.tmin is None: a.tmin = max(r["t"][0] for r in RUNS)
peak = lambda k: max(r["eps"].get(k, np.zeros(1)).max() for r in RUNS)
modes = sorted(set().union(*[set(r["eps"]) for r in RUNS]), key=lambda k: -peak(k))[:a.max_modes]; modes.sort()
lab = lambda m: "".join(map(str, m)) if m[0] < 10 else ",".join(map(str, m))
import itertools; fb = itertools.cycle(FALLBACK); col = {m: (COLORS[lab(m)] if lab(m) in COLORS else next(fb)) for m in modes}

fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(8.0, 7.6), sharex=True, gridspec_kw=dict(height_ratios=[1, 3.0], hspace=0.05))
for run, ls, lw in zip(RUNS, ("-", "--"), (1.8, 1.8)):
    ax0.plot(run["t"], run["spin"], color="k", ls=ls, lw=lw)
    for m in modes:
        if m in run["eps"]:
            sel = run["t"] >= a.tmin
            ax1.plot(run["t"][sel], np.maximum(run["eps"][m][sel] * sc, 1e-300), color=col[m], ls=ls, lw=lw)
# top panel: spin
ax0.set_ylim(0, 1.08); ax0.set_yticks([0, 0.5, 1]); ax0.yaxis.set_minor_locator(AutoMinorLocator(5))
ax0.set_ylabel(r"$\tilde a$")
# bottom panel: occupations
ax1.set_xscale("log"); ax1.set_yscale("log"); ax1.set_ylim(a.eps_min * sc, 3.0 * sc)
ax1.set_xlim(a.tmin, a.tmax if a.tmax else max(r["t"][-1] for r in RUNS))
ax1.yaxis.set_major_locator(LogLocator(base=10, numticks=12)); ax1.yaxis.set_minor_locator(LogLocator(base=10, subs="auto", numticks=100))
ax1.xaxis.set_minor_locator(LogLocator(base=10, subs="auto", numticks=100))
ax1.set_xlabel("t  [yrs]"); ax1.set_ylabel(r"$\epsilon_x$" if a.quantity == "eps" else r"$M_{c,\,nlm}\,/\,M_0$")
ncol = 1 if len(modes) <= 8 else 2; nrow = -(-len(modes) // ncol)
if a.label_pos == "auto":                  # pick the emptier corner: count curve points inside each candidate label block
    x0, x1 = np.log10(a.tmin), np.log10(a.tmax if a.tmax else max(r["t"][-1] for r in RUNS)); y0, y1 = np.log10(a.eps_min), np.log10(3.0)
    bw, bh = 0.03 + cw * ncol + 0.02, 0.04 + 0.055 * nrow + 0.02
    cnt = {"lower left": 0, "upper left": 0}
    for run in RUNS:
        for m in modes:
            if m in run["eps"]:
                e = run["eps"][m]; sel = (run["t"] >= a.tmin) & (e > a.eps_min)
                fx = (np.log10(run["t"][sel]) - x0) / (x1 - x0); fy = (np.log10(e[sel]) - y0) / (y1 - y0)
                inx = fx < bw; cnt["lower left"] += int(np.sum(inx & (fy < bh))); cnt["upper left"] += int(np.sum(inx & (fy > 1 - bh)))
    a.label_pos = "lower left" if cnt["lower left"] <= cnt["upper left"] else "upper left"
    print("label corner:", a.label_pos, cnt)
for i, m in enumerate(modes):            # coloured text labels, paper style
    c, r = divmod(i, nrow)
    y = 0.04 + 0.055 * (nrow - 1 - r) if a.label_pos == "lower left" else 0.96 - 0.055 * r
    ax1.text(0.03 + cw * c, y, rf"$\epsilon_{{{lab(m)}}}$" if a.quantity == "eps" else lab(m), transform=ax1.transAxes, color=col[m], fontsize=15,
             va="bottom" if a.label_pos == "lower left" else "top")
for ax in (ax0, ax1):
    ax.tick_params(which="both", top=True, right=True, direction="in")
out = a.out or f"figs/{'states' if a.quantity == 'eps' else 'mc'}_BH{a.bh}_alpha{a.alpha}_{a.solid}_vs_{a.dashed}.png"
pathlib.Path(out).parent.mkdir(exist_ok=True); fig.savefig(out, dpi=150, bbox_inches="tight")
print("font:", pathlib.Path(font_manager.findfont(font_manager.FontProperties(family=["sans-serif"]))).name)
print("modes shown:", [lab(m) for m in modes])
for run in RUNS:
    print(f"{run['name']}: {len(run['t'])} points, final spin {run['spin'][-1]:.4f}, final mass {run['mass'][-1]:.3f}, modes>eps_min: {[lab(k) for k in sorted(run['eps'])]}")
print("wrote", out)
