"""
Plot the GW lines of an ensemble written by run_forward.jl, weighting each
realization by how well its BH mass and spin *today* match a measurement.

    python plot_realizations.py <outdir> <tag> [--M_obs 10 --M_obs_sigma 1]
                                [--a_obs 0.35 --a_obs_sigma 0.05]
                                [--n_labels 8] [--h_min 1e-35]

Weights: if <outdir>/weights_<tag>.dat exists (written by run_zoom.jl, which
combines its zoom rounds by importance sampling) and no --M_obs/--a_obs is
given, those weights and that measurement are used. Otherwise the realizations
are taken to be draws from the birth prior, so the posterior weight is the
likelihood of today's values,
    w = N(M_today; M_obs, M_obs_sigma) N(a_today; a_obs, a_obs_sigma)
        × [a_obs_min < a_today < a_obs_max],
for whichever factors are given (w = 1 if none). Realizations that
did not reach t = age get w = 0. Marker opacity is proportional to w / max w;
histograms are weighted by w. The effective sample size
(sum w)^2 / sum w^2 is printed and shown in the title; with tight measurements
only a few prior draws survive, so run more realizations.

Panels:
  (a) birth mass and spin;   (b) mass and spin today, with the 1 and 2 sigma box;
  (c) h(iota) vs f for every line, coloured by line (the n_labels lines that
      occur in most realizations; colours do not depend on the weights);
  (d) weighted distribution of log10 h(iota) for each of those lines;
  (e) each line's frequency about its weighted median (points; bar = weighted
      16-84% range);
  (f) h(iota) vs spin today;
  (g), (h) weighted distribution of log10 age and log10 f_Edd (posterior),
      with the unweighted round-1 (prior) draws for comparison;
  (i) h(iota) vs age.
With a zoom run, (a) colours the birth parameters by round and draws each
round's box.
Lines with h(iota) < h_min (levels left at their seed occupation) are not
shown. Saves <outdir>/realizations_<tag>[_<suffix>].png.
"""
import argparse
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def table(fn):
    # headers start with "# "; read the column names explicitly
    with open(fn) as f:
        names = f.readline().lstrip("#").split()
    return np.atleast_1d(np.genfromtxt(fn, names=names, dtype=None, encoding=None, comments="#"))


def label(r):
    a = f"{r['na']}{r['la']}{r['ma']}"
    b = f"{r['nb']}{r['lb']}{r['mb']}"
    return f"{a}x{b}" if r["kind"] == "annihilation" else f"{a}->{b}"


def weighted_median(x, w):
    i = np.argsort(x)
    c = np.cumsum(w[i])
    return x[i][np.searchsorted(c, 0.5 * c[-1])]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("tag")
    ap.add_argument("--M_obs", type=float, default=None)
    ap.add_argument("--M_obs_sigma", type=float, default=1.0)
    ap.add_argument("--a_obs", type=float, default=None)
    ap.add_argument("--a_obs_sigma", type=float, default=0.05)
    ap.add_argument("--a_obs_min", type=float, default=None)
    ap.add_argument("--a_obs_max", type=float, default=None)
    ap.add_argument("--n_labels", type=int, default=8)
    ap.add_argument("--h_min", type=float, default=1e-35)
    ap.add_argument("--suffix", default="", help="appended to the output file name")
    args = ap.parse_args()

    samp = table(os.path.join(args.outdir, f"samples_{args.tag}.dat"))
    lines = table(os.path.join(args.outdir, f"lines_{args.tag}.dat"))

    ok = samp["status"] == "ok"
    ids = samp["id"].astype(int)
    rnd = np.ones(len(samp), dtype=int)
    boxes = None
    fw = os.path.join(args.outdir, f"weights_{args.tag}.dat")
    given = [args.M_obs, args.a_obs, args.a_obs_min, args.a_obs_max]
    if os.path.isfile(fw) and all(g is None for g in given):
        with open(fw) as f:
            hdr = f.readline()
        m = re.search(r"M_obs = (\S+) \+- (\S+), a_obs = (\S+) \+- (\S+)[,;]", hdr)
        val = lambda x: None if x.lower() == "nan" else float(x)
        args.M_obs, args.M_obs_sigma, args.a_obs, args.a_obs_sigma = (val(m[i]) for i in range(1, 5))
        m = re.search(r"a_range = (\S+) (\S+);", hdr)
        if m:
            args.a_obs_min, args.a_obs_max = val(m[1]), val(m[2])
        wt = np.loadtxt(fw, comments="#", ndmin=2)
        wmap = {int(r[0]): (int(r[1]), r[2]) for r in wt}
        rnd = np.array([wmap[i][0] if i in wmap else 0 for i in ids])
        w = np.array([wmap[i][1] if i in wmap else 0.0 for i in ids])
        w = np.where(ok, w, 0.0)
        fr = os.path.join(args.outdir, f"rounds_{args.tag}.dat")
        if os.path.isfile(fr):
            boxes = np.loadtxt(fr, comments="#", ndmin=2)
    else:
        logw = np.zeros(len(samp))
        if args.M_obs is not None:
            logw += -0.5 * ((samp["M_today_Msun"] - args.M_obs) / args.M_obs_sigma) ** 2
        if args.a_obs is not None:
            logw += -0.5 * ((samp["a_today"] - args.a_obs) / args.a_obs_sigma) ** 2
        if args.a_obs_min is not None:
            logw[samp["a_today"] <= args.a_obs_min] = -np.inf
        if args.a_obs_max is not None:
            logw[samp["a_today"] >= args.a_obs_max] = -np.inf
        good = ok & np.isfinite(logw)
        w = np.where(good, np.exp(logw - (np.max(logw[good]) if good.any() else 0.0)), 0.0)
    ess = w.sum() ** 2 / np.sum(w ** 2) if w.sum() > 0 else 0.0
    print(f"{len(samp)} realizations, {ok.sum()} ok, effective sample size {ess:.1f}")
    wid = dict(zip(samp["id"].astype(int), w))
    aid = dict(zip(samp["id"].astype(int), samp["a_today"]))
    ageid = dict(zip(samp["id"].astype(int), samp["age_yr"]))

    lines = lines[lines["h_iota"] >= args.h_min]
    lw = np.array([wid[int(i)] for i in lines["id"]])
    keep = lw > 0
    lines, lw = lines[keep], lw[keep]
    labs = np.array([label(r) for r in lines])
    # fixed colours: ordered by number of realizations with that line (unweighted, all draws)
    alllabs = np.array([label(r) for r in table(os.path.join(args.outdir, f"lines_{args.tag}.dat"))
                        if r["h_iota"] >= args.h_min])
    u, cnt = np.unique(alllabs, return_counts=True)
    top = [l for l in u[np.argsort(-cnt, kind="stable")][:args.n_labels] if l in set(labs)]
    cols = plt.cm.tab10(np.arange(len(top)) % 10)

    def rgba(c, ws):
        out = np.tile(c, (len(ws), 1))
        out[:, 3] = 0.05 + 0.95 * ws / max(w.max(), 1e-300)
        return out

    fig, ax = plt.subplots(3, 3, figsize=(17, 14))
    ax = ax.ravel()

    # (a), (b): birth and today
    sel = ok & (w > 0)
    ax[0].scatter(samp["M0_Msun"][~ok], samp["a0"][~ok], marker="x", c="r", s=15, label="incomplete")
    for r in np.unique(rnd[ok & (rnd > 0)]):     # rnd = 0: run still in progress, no weight yet
        s_ = ok & (rnd == r)
        c = np.array([0, 0, 0, 1.0]) if boxes is None else np.array(plt.cm.viridis((r - 1) / max(rnd.max() - 1, 1)))
        ax[0].scatter(samp["M0_Msun"][s_], samp["a0"][s_], s=14, color=rgba(c, w[s_]))
        if boxes is not None and (boxes[:, 0] == r).any():
            b = boxes[boxes[:, 0] == r][0]
            ax[0].add_patch(plt.Rectangle((b[2], b[4]), b[3] - b[2], b[5] - b[4], fill=False, ec=c, lw=1.5,
                                          label=f"round {r} ({int(b[1])})"))
    if boxes is not None:
        ax[0].legend(fontsize=7, loc="lower right")
    ax[0].set(xlabel=r"$M_0\ [M_\odot]$", ylabel=r"$a_{*,0}$", title="(a) birth")
    ax[1].scatter(samp["M_today_Msun"][ok], samp["a_today"][ok], s=14, color=rgba(np.array([0, 0, 0, 1.0]), w[ok]))
    for k, ls in ((1, "-"), (2, "--")):
        if args.M_obs is not None and args.a_obs is not None:
            ax[1].add_patch(plt.Rectangle((args.M_obs - k * args.M_obs_sigma, args.a_obs - k * args.a_obs_sigma),
                                          2 * k * args.M_obs_sigma, 2 * k * args.a_obs_sigma,
                                          fill=False, ec="C3", ls=ls))
        elif args.M_obs is not None:
            ax[1].axvspan(args.M_obs - k * args.M_obs_sigma, args.M_obs + k * args.M_obs_sigma, color="C3", alpha=0.1)
        elif args.a_obs is not None:
            ax[1].axhspan(args.a_obs - k * args.a_obs_sigma, args.a_obs + k * args.a_obs_sigma, color="C3", alpha=0.1)
    for lim in (args.a_obs_min, args.a_obs_max):
        if lim is not None:
            ax[1].axhline(lim, color="C3", lw=1.5)
    ax[1].set(xlabel=r"$M\ [M_\odot]$ today", ylabel=r"$a_*$ today", title="(b) today (red: measurement)")

    # (c)-(f): lines
    other = ~np.isin(labs, top)
    ax[2].scatter(lines["f_Hz"][other], lines["h_iota"][other], s=6, color=rgba(np.array([0.6, 0.6, 0.6, 1.0]), lw[other]))
    for c, l in zip(cols, top):
        s = labs == l
        f, h, ws = lines["f_Hz"][s], lines["h_iota"][s], lw[s]
        ax[2].scatter(f, h, s=12, color=rgba(c, ws))
        ax[2].scatter([], [], s=12, color=c, label=l)
        ax[3].hist(np.log10(h), bins=30, weights=ws, histtype="step", color=c, label=l)
        fm = weighted_median(f, ws)
        k = top.index(l)
        df = (f - fm) / fm
        jit = k + 0.15 * np.random.default_rng(k).uniform(-1, 1, len(df))
        ax[4].scatter(jit, df, s=12, color=rgba(c, ws))
        if ws.sum() > 0:
            i = np.argsort(df)
            cw = np.cumsum(ws[i]) / ws.sum()
            lo_, hi_ = df[i][np.searchsorted(cw, 0.16)], df[i][min(np.searchsorted(cw, 0.84), len(df) - 1)]
            ax[4].plot([k + 0.3, k + 0.3], [lo_, hi_], color=c, lw=3)
        ax[4].text(k, 1.0, f"{fm:.4g} Hz", rotation=90, ha="center", va="bottom", fontsize=7,
                   transform=ax[4].get_xaxis_transform())
        ax[5].scatter([aid[int(i)] for i in lines["id"][s]], h, s=12, color=rgba(c, ws))
        ax[8].scatter([ageid[int(i)] for i in lines["id"][s]], h, s=12, color=rgba(c, ws))
    ax[2].set(xscale="log", yscale="log", xlabel="f [Hz]", ylabel=r"$h(\iota)$", title="(c) lines today")
    ax[2].legend(fontsize=7, ncol=2)
    ax[3].set(xlabel=r"$\log_{10} h(\iota)$", ylabel="weight", title="(d) strain")
    ax[4].set_yscale("symlog", linthresh=1e-6)
    ax[4].set_xticks(range(len(top)))
    ax[4].set_xticklabels(top, rotation=45, fontsize=8)
    ax[4].set(ylabel=r"$(f - \tilde f)/\tilde f$")
    ax[4].set_title(r"(e) frequency spread ($\tilde f$ = weighted median)", pad=45)
    ax[5].set(yscale="log", xlabel=r"$a_*$ today", ylabel=r"$h(\iota)$", title="(f) strain vs spin today")
    ax[8].set(xscale="log", yscale="log", xlabel="age [yr]", ylabel=r"$h(\iota)$", title="(i) strain vs age")

    # (g), (h): age and accretion rate, posterior vs prior (round-1 draws)
    pri = rnd == rnd[rnd > 0].min() if (rnd > 0).any() else np.ones(len(samp), bool)
    for a_, col, xl, ttl in ((ax[6], "age_yr", r"$\log_{10}$ age [yr]", "(g) age"),
                             (ax[7], "f_edd", r"$\log_{10} f_{\rm Edd}$", "(h) accretion rate")):
        v = np.log10(np.maximum(samp[col], 1e-300))
        if np.ptp(v[ok]) == 0:
            a_.text(0.5, 0.5, f"fixed: {samp[col][0]:.3g}", transform=a_.transAxes, ha="center")
        else:
            bins = np.linspace(v[ok].min(), v[ok].max(), 21)
            a_.hist(v[pri & ok], bins=bins, density=True, histtype="stepfilled", color="0.85", label="prior (round 1)")
            if w.sum() > 0:
                a_.hist(v[ok], bins=bins, weights=w[ok], density=True, histtype="step", color="k", lw=2, label="posterior")
            a_.legend(fontsize=8)
        a_.set(xlabel=xl, ylabel="density", title=ttl)

    obs = []
    if args.M_obs is not None:
        obs.append(f"M = {args.M_obs} ± {args.M_obs_sigma}")
    if args.a_obs is not None:
        obs.append(f"a = {args.a_obs} ± {args.a_obs_sigma}")
    if args.a_obs_min is not None:
        obs.append(f"a > {args.a_obs_min}")
    if args.a_obs_max is not None:
        obs.append(f"a < {args.a_obs_max}")
    fig.suptitle(f"{args.tag}: {len(samp)} realizations ({ok.sum()} ok), "
                 f"weighted by {', '.join(obs) if obs else 'nothing (prior)'}; ESS = {ess:.1f}")
    fig.tight_layout()
    out = os.path.join(args.outdir, f"realizations_{args.tag}{'_' + args.suffix if args.suffix else ''}.png")
    fig.savefig(out, dpi=130)
    print("saved", out)


if __name__ == "__main__":
    main()
