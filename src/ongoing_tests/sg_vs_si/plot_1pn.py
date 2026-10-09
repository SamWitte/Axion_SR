"""Plot the O(alpha^2) corrections delta_X = A_X / A_N to the self-gravity 2->2 amplitude (output of sg_1pn.jl).

columns of out_1pn/<channel>_<i>.dat:  alpha, Re/Im delta for GM (h0i), hij, T00, ret, BHc
Also writes note_1pn/coeffs.tex (table of delta / alpha^2) and note_1pn/fig_1pn.pdf.
"""
import glob
import os

import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
NOTE = os.path.join(HERE, "note_1pn")
os.makedirs(NOTE, exist_ok=True)
CHANNELS = {
    "211_211_322_BH": r"$211\times211\to322\times{\rm BH}$",
    "322_322_544_BH": r"$322\times322\to544\times{\rm BH}$",
    "322_322_211_Inf": r"$322\times322\to211\times\infty$",
    "544_544_322_Inf": r"$544\times544\to322\times\infty$",
}
PIECES = [("GM", r"$h_{0i}$ (gravitomagnetic)"), ("hij", r"$h_{ij}$ (stress)"), ("T00", r"$T_{00}$ kinetic"),
          ("ret", "retardation"), ("BHc", r"BH$\times$cloud (EIH $G^2$)")]
COLORS = {"GM": "#2a78d6", "hij": "#eb6834", "T00": "#1baf7a", "ret": "#eda100", "BHc": "#e87ba4"}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "cm", "font.size": 10,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.linewidth": 0.8, "xtick.direction": "in", "ytick.direction": "in",
})


def load(ch):
    files = sorted(glob.glob(os.path.join(HERE, "out_1pn", f"{ch}_[0-9][0-9].dat")))
    d = np.array([np.loadtxt(f) for f in files])
    d = d[np.argsort(d[:, 0])]
    out = {"alpha": d[:, 0]}
    for j, (k, _) in enumerate(PIECES):
        out[k] = d[:, 1 + 2 * j] + 1j * d[:, 2 + 2 * j]
    out["hij_full"] = 2 * out["hij"]                          # h_ij + equal stress-sourced h_00 piece
    out["sum"] = out["GM"] + out["hij_full"] + out["T00"] + out["ret"] + out["BHc"]
    return out


def kerr_vs_hyd(ch):
    s = np.loadtxt(os.path.join(HERE, f"sg_vs_si_{ch}_summary.dat"))
    return s[:, 0], np.sqrt(s[:, 1] / s[:, 2]) - 1.0              # amplitude ratio Kerr/hydrogenic - 1


fig, axs = plt.subplots(2, 2, figsize=(10.5, 7.2), constrained_layout=True)
rows = []
for ax, (ch, label) in zip(axs.flat, CHANNELS.items()):
    d = load(ch)
    a = d["alpha"]
    for k, lab in PIECES:
        ax.plot(a, 100 * d[k].real, color=COLORS[k], lw=2.2 if k in ("GM", "hij") else 1.3, label=lab)
    xa, kh = kerr_vs_hyd(ch)
    ax.plot(xa, 100 * kh, color=INK2, lw=1.2, ls=(0, (4, 3)), label="Newtonian: Kerr vs hydrogenic")
    ax.axhline(0, color=INK2, lw=0.6)
    ax.set_title(label, color=INK, fontsize=10.5, loc="left")
    ax.set_xlabel(r"$\alpha = GM\mu$")
    ax.set_ylabel(r"$\delta_X = A_X/A_{\rm N}$  [%]")
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_xlim(0, a[-1] * 1.02)
    # coefficients delta/alpha^2 at the smallest alpha and at alpha ~ 0.3
    i3 = np.argmin(abs(a - 0.3))
    rows.append((ch, label, a[0], {k: d[k].real[0] / a[0] ** 2 for k in ["GM", "hij", "T00", "ret", "BHc", "sum"]},
                 a[i3], {k: d[k].real[i3] for k in ["GM", "hij", "T00", "ret", "BHc", "sum"]}))
axs.flat[0].legend(fontsize=8, frameon=True, facecolor="white", edgecolor=GRID, loc="lower left")
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(NOTE, f"fig_1pn.{ext}"), dpi=200)

# LaTeX table of coefficients
with open(os.path.join(NOTE, "coeffs.tex"), "w") as f:
    for ch, label, a0, c0, a3, c3 in rows:
        f.write(f"{label} & " + " & ".join(f"${c0[k]:+.2f}$" for k in ["GM", "hij", "T00", "ret", "BHc", "sum"]) + r" \\" + "\n")
with open(os.path.join(NOTE, "values_at_03.tex"), "w") as f:
    for ch, label, a0, c0, a3, c3 in rows:
        rate = abs(1 + c3["GM"]) ** 2 - 1
        rate_h = abs(1 + c3["hij"]) ** 2 - 1
        f.write(f"{label} & ${a3:.2f}$ & " + " & ".join(f"${100*c3[k]:+.1f}$" for k in ["GM", "hij"])
                + f" & ${100*rate:+.0f}$ & ${100*rate_h:+.0f}$ & ${100*c3['sum']:+.1f}$" + r" \\" + "\n")
for ch, label, a0, c0, a3, c3 in rows:
    print(ch, "coeff(alpha^2) at alpha0:", {k: round(v, 3) for k, v in c0.items()}, " values at alpha=", round(a3, 3), {k: round(100 * v, 2) for k, v in c3.items()})
print("saved")
