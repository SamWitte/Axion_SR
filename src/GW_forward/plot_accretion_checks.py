"""Plot the output of accretion_checks.jl -> test_plots/accretion_checks<sfx>.png

    python plot_accretion_checks.py [Nmax=3]
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_plots")
Nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 3
sfx = "" if Nmax == 3 else f"_Nmax{Nmax}"
fig, ax = plt.subplots(2, 3, figsize=(16, 9))

# (a) Bardeen spin-up, no SR
b = np.loadtxt(os.path.join(d, f"bardeen{sfx}.dat"))
x = b[:, 2] / b[0, 2]
xs = np.linspace(1, np.sqrt(6), 300)
ax[0, 0].plot(xs, np.sqrt(2 / 3) / xs * (4 - np.sqrt(18 / xs**2 - 2)), "k-", lw=3, alpha=0.3, label="Bardeen (1970), $a_0=0$")
ax[0, 0].plot(x, b[:, 1], "C3--", label=r"solver, $a_0=0.01$, $f_{\rm Edd}=1$, $\alpha=0.005$")
ax[0, 0].set(xlabel=r"$M/M_0$", ylabel=r"$a_*$", title="(a) accretion only (SR negligible)")
ax[0, 0].legend(fontsize=8)
ax[1, 0].plot(b[:, 0], b[:, 2], "C3-", label="solver")
ax[1, 0].plot(b[:, 0], b[0, 2] * (1 + b[:, 0] / (0.1 * 4.5133e8)), "k:", label=r"$M_0 + \dot M t$")
ax[1, 0].set(xlabel="t [yr]", ylabel=r"$M\ [M_\odot]$", title="(a) mass growth")
ax[1, 0].legend(fontsize=8)

# (b) SR + accretion
runs = [("fedd0", r"$f_{\rm Edd}=0$ (old default)", "k", "-"),
        ("fedd0_track", r"$f_{\rm Edd}=0$, track_alpha", "C0", "--"),
        ("fedd1e-6", r"$f_{\rm Edd}=10^{-6}$", "C1", ":"),
        ("fedd0.05", r"$f_{\rm Edd}=0.05$", "C2", "-"),
        ("fedd0.5", r"$f_{\rm Edd}=0.5$", "C3", "-")]
for name, lab, c, ls in runs:
    r = np.loadtxt(os.path.join(d, f"sr_{name}{sfx}.dat"))
    t = r[:, 0]
    ax[0, 1].plot(t, r[:, 1], color=c, ls=ls, label=lab)
    ax[0, 2].plot(t, r[:, 2], color=c, ls=ls, label=lab)
    ax[1, 1].plot(t, r[:, 3], color=c, ls=ls, label=lab)   # 211
    ax[1, 2].plot(t, r[:, 6], color=c, ls=ls, label=lab)   # 322
for a, yl, ttl in ((ax[0, 1], r"$a_*$", "(b) spin"), (ax[0, 2], r"$M\ [M_\odot]$", "(b) BH mass"),
                   (ax[1, 1], r"$u_{211}$", "(b) 211 occupation"), (ax[1, 2], r"$u_{322}$", "(b) 322 occupation")):
    a.set(xscale="log", xlabel="t [yr]", ylabel=yl, title=ttl, xlim=(1e1, 1e8))
for a in (ax[1, 1], ax[1, 2]):
    a.set(yscale="log", ylim=(1e-6, 1))
ax[0, 1].legend(fontsize=8)
fig.suptitle(r"Accretion checks.  (b): $M_0=10\,M_\odot$, $a_0=0.9$, $\alpha_0=0.1$, $N_{\max}=$" + str(Nmax))
fig.tight_layout()
fig.savefig(os.path.join(d, f"accretion_checks{sfx}.png"), dpi=130)
print("saved", os.path.join(d, f"accretion_checks{sfx}.png"))
