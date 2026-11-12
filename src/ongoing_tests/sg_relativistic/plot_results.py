"""Relativistic one-graviton exchange vs Newtonian self-gravity, per channel (results_prod.dat from scan.jl).

usage: python3 plot_results.py results_prod.dat out.pdf
"""
import sys
import numpy as np
import matplotlib.pyplot as plt

INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
C_GR, C_NR = "#2a78d6", "#eb6834"
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 8,
                     "axes.edgecolor": INK2, "xtick.color": INK2, "ytick.color": INK2, "axes.labelcolor": INK,
                     "xtick.direction": "in", "ytick.direction": "in"})

rows = [l.split() for l in open(sys.argv[1]) if l.strip() and not l.startswith("#")]
chans = [("211", "322", "BH"), ("322", "544", "BH"), ("322", "211", "Inf"), ("544", "322", "Inf")]
label = {"BH": r"\mathrm{BH}", "Inf": r"\infty"}
fig, axs = plt.subplots(2, 4, figsize=(7.0, 3.6), sharex="col", gridspec_kw=dict(height_ratios=[1.4, 1], hspace=0.08, wspace=0.32))
for j, (s1, s3, d) in enumerate(chans):
    sel = [r for r in rows if (r[0], r[1], r[2]) == (s1, s3, d)]
    sel.sort(key=lambda r: float(r[3]))
    if not sel:
        continue
    a = np.array([float(r[3]) for r in sel])
    gr = np.array([complex(float(r[8]), float(r[9])) for r in sel])
    nr = np.array([complex(float(r[10]), float(r[11])) for r in sel])
    ax, axr = axs[0, j], axs[1, j]
    ax.plot(a, np.abs(gr) ** 2, "-o", ms=2.5, lw=1.2, color=C_GR, label="exact (Kerr)")
    ax.plot(a, np.abs(nr) ** 2, "--", lw=1.1, color=C_NR, label="Newton + recoil")
    ax.axhline(1, color=INK2, lw=0.6, ls=":")
    ax.set_title(r"$%s\times%s\to%s\times%s$" % (s1, s1, s3, label[d]), fontsize=8)
    axr.plot(a, np.abs(gr / nr) ** 2 - 1, "-o", ms=2.5, lw=1.2, color=C_GR)
    axr.axhline(0, color=INK2, lw=0.6, ls=":")
    for x in (ax, axr):
        x.grid(True, color=GRID, lw=0.5); x.set_axisbelow(True)
    axr.set_xlabel(r"$\alpha$")
    if j == 0:
        ax.set_ylabel(r"$\Gamma_{\rm SG}/\Gamma_{\rm SG}^{\rm Newton}$")
        axr.set_ylabel(r"$\Gamma/\Gamma_{\rm N+recoil}-1$")
axs[0, 0].legend(frameon=False, fontsize=7, loc="center left")
fig.savefig(sys.argv[2], bbox_inches="tight")
fig.savefig(sys.argv[2].rsplit(".", 1)[0] + ".png", dpi=170, bbox_inches="tight")

print(f"{'channel':22s} {'alpha':>7s} {'A_GR/A_N':>18s} {'A_NR/A_N':>10s} {'|GR/N|^2':>9s} {'|GR/NR|^2':>9s}")
for s1, s3, d in chans:
    for r in sorted([r for r in rows if (r[0], r[1], r[2]) == (s1, s3, d)], key=lambda r: float(r[3])):
        g = complex(float(r[8]), float(r[9])); n = complex(float(r[10]), float(r[11]))
        print(f"{s1}x{s1}->{s3}x{d:4s}       {float(r[3]):7.4f} {g.real:+9.5f}{g.imag:+8.5f}i {n.real:+10.5f} {abs(g)**2:9.4f} {abs(g/n)**2:9.4f}")
