#!/usr/bin/env python3
"""
Combine AGN companion-resonance runs at different boson masses (run_agn_scan.sh).

    python agn_scan_summary.py --out scan_runs scan_runs/alpha* [--M_ref 1e6] [--tag quad]

Reads <rundir>/params.json and <rundir>/first_resonances[_tag].csv of every run
directory (alpha_ref = G M_ref mu) and writes
  <out>/scan_kinds[_tag].csv        alpha_ref, mu_eV, channel, N, fraction of every outcome
  <out>/scan_transitions[_tag].csv  alpha_ref, mu_eV, channel, first transition, count, fraction
  <out>/scan_summary[_tag].png      left: outcome fractions vs alpha_ref; right: fraction of
                                    companions whose first line is each transition (heatmap)
"""
import argparse
import collections
import csv
import json
import os

import numpy as np

import agn_cloud_resonances as A

GROUPS = [("hyperfine", ("hyperfine",)), ("fine", ("fine",)), ("Bohr", ("bohr",)),
          ("no line reached", ("stalled", "none_ahead")), ("no cloud", ("no_SR_state",))]
GROUP_COL = {"hyperfine": A.KIND_COL["hyperfine"], "fine": A.KIND_COL["fine"],
             "Bohr": A.KIND_COL["bohr"], "no line reached": "#c3c2b7", "no cloud": A.MUTED}
BLUES = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
         "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]


def load_runs(dirs, M_ref, tag):
    runs = []
    for d in dirs:
        fn = os.path.join(d, f"first_resonances{'_' + tag if tag else ''}.csv")
        if not (os.path.exists(fn) and os.path.exists(os.path.join(d, "params.json"))):
            print(f"[skip] {d}: no {os.path.basename(fn)}")
            continue
        mu = json.load(open(os.path.join(d, "params.json")))["mu_eV"]
        rows = list(csv.DictReader(open(fn)))
        if not rows:
            print(f"[skip] {d}: no companions")
            continue
        runs.append(dict(dir=d, mu=mu, alpha=A.alpha_of(M_ref * A.Msun, mu), rows=rows))
    return sorted(runs, key=lambda r: r["alpha"])


def kind_of(row):
    return row["outcome"] if row["outcome"] in ("hyperfine", "fine", "bohr") else None


def write_tables(runs, out, tag):
    sfx = f"_{tag}" if tag else ""
    with open(os.path.join(out, f"scan_kinds{sfx}.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["alpha_ref", "mu_eV", "channel", "N"] + A.OUTCOMES)
        for r in runs:
            for ch in ["all"] + sorted({x["channel"] for x in r["rows"]}):
                sub = [x for x in r["rows"] if ch == "all" or x["channel"] == ch]
                w.writerow([r["alpha"], r["mu"], ch, len(sub)]
                           + [np.mean([x["outcome"] == k for x in sub]) for k in A.OUTCOMES])
    with open(os.path.join(out, f"scan_transitions{sfx}.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["alpha_ref", "mu_eV", "channel", "first", "kind", "count", "fraction"])
        for r in runs:
            for ch in ["all"] + sorted({x["channel"] for x in r["rows"]}):
                sub = [x for x in r["rows"] if ch == "all" or x["channel"] == ch]
                for lab, n in collections.Counter(x["first"] for x in sub).most_common():
                    kind = next(x["outcome"] for x in sub if x["first"] == lab)
                    w.writerow([r["alpha"], r["mu"], ch, lab, kind, n, n / len(sub)])


def tex_label(lab):
    """'322->320/32m2' -> mathtext kets."""
    def ket(s):
        out, i = "", 0
        while i < len(s):
            if s[i] == "m":
                out += r"\bar{%s}" % s[i + 1]
                i += 2
            else:
                out += s[i]
                i += 1
        return r"$|" + out + r"\rangle$"
    init, finals = lab.split("->")
    return ket(init) + r" $\rightarrow$ " + " / ".join(ket(f) for f in finals.split("/"))


def plot(runs, out, tag, M_ref, top=15):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.patches import Patch
    al = np.array([r["alpha"] for r in runs])
    # transitions ranked by their largest fraction over the scan
    frac = {}
    for j, r in enumerate(runs):
        c = collections.Counter(x["first"] for x in r["rows"] if kind_of(x))
        for lab, n in c.items():
            frac.setdefault(lab, np.zeros(len(runs)))[j] = n / len(r["rows"])
    kind = {x["first"]: x["outcome"] for r in runs for x in r["rows"] if kind_of(x)}
    labs = sorted(frac, key=lambda l: -frac[l].max())[:top]
    n_rows = max(len(labs), 6)
    fig = plt.figure(figsize=(12.5, 1.6 + 0.32 * n_rows))
    fig.patch.set_facecolor("#fcfcfb")
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.25], wspace=0.55)
    ax, hx = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    for a_ in (ax, hx):
        a_.set_facecolor("#fcfcfb")

    # left: outcome fractions vs alpha_ref
    for name, outs in GROUPS:
        f = np.array([np.mean([x["outcome"] in outs for x in r["rows"]]) for r in runs])
        ax.plot(al, f, "-o", color=GROUP_COL[name], lw=2, ms=6.5, mec="#fcfcfb", mew=2,
                label=name, zorder=3, solid_capstyle="round")
    ax.set_xscale("log")
    ax.set_xticks(al)
    ax.set_xticklabels([f"{a:.3g}" for a in al])
    ax.minorticks_off()
    ax.set_ylim(0, 1.02)
    ax.set_xlabel(rf"$\alpha$ at $M = 10^{{{np.log10(M_ref):.0f}}}\,M_\odot$", fontsize=9, color=A.INK2)
    ax.set_ylabel("fraction of companions", fontsize=9, color=A.INK2)
    ax.set_title("first outcome", fontsize=9.5, color=A.INK, loc="left")
    ax.grid(axis="y", color=A.GRID, lw=0.8, zorder=0)
    ax.tick_params(colors=A.MUTED, labelsize=8, length=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(A.AXIS)
    ax.legend(fontsize=8, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0), labelcolor=A.INK2)

    # right: transition x alpha heatmap (one hue; empty cells neutral)
    if labs:
        Z = np.array([frac[l] for l in labs])
        cmap = LinearSegmentedColormap.from_list("blues", BLUES)
        cmap.set_bad("#f0efec")
        vmax = Z.max()
        hx.imshow(np.ma.masked_equal(Z, 0), cmap=cmap, vmin=0, vmax=vmax, aspect="auto",
                  interpolation="nearest")
        for i in range(Z.shape[0]):
            for j in range(Z.shape[1]):
                if Z[i, j] >= 0.01:
                    hx.text(j, i, f"{100 * Z[i, j]:.0f}", ha="center", va="center", fontsize=7.5,
                            color="white" if Z[i, j] > 0.55 * vmax else A.INK)
        hx.set_yticks(range(len(labs)))
        hx.set_yticklabels([tex_label(l) for l in labs], fontsize=8.5, color=A.INK)
        # kind as a fixed-size swatch between label and map (x in axes units, y in rows)
        from matplotlib.transforms import blended_transform_factory
        hx.scatter(np.full(len(labs), -0.025), np.arange(len(labs)), s=55, marker="s",
                   c=[A.KIND_COL[kind[l]] for l in labs], clip_on=False, zorder=3,
                   transform=blended_transform_factory(hx.transAxes, hx.transData))
        hx.set_xticks(range(len(runs)))
        hx.set_xticklabels([f"{a:.3g}" for a in al])
        hx.set_xlabel(r"$\alpha$ at $M_{\rm ref}$", fontsize=9, color=A.INK2)
        hx.tick_params(colors=A.MUTED, labelsize=8, length=0, pad=6)
        hx.tick_params(axis="y", colors=A.INK, pad=22)
        for s in hx.spines.values():
            s.set_visible(False)
        hx.set_title("first line crossed  (% of companions)", fontsize=9.5, color=A.INK, loc="left")
        hx.legend(handles=[Patch(fc=A.KIND_COL[k], label="Bohr" if k == "bohr" else k)
                           for k in ("hyperfine", "fine", "bohr")],
                  fontsize=8, frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0),
                  labelcolor=A.INK2, handlelength=1)
    else:
        hx.axis("off")
    N = [len(r["rows"]) for r in runs]
    fig.text(0.01, 0.985, "Companion resonances across boson masses", fontsize=10.5, weight="bold",
             color=A.INK, va="top")
    fig.text(0.01, 0.985 - 0.24 / fig.get_figheight(),
             f"{len(runs)} runs, {min(N)}-{max(N)} companions each{'  (' + tag + ')' if tag else ''}",
             fontsize=8.5, color=A.INK2, va="top")
    fig.subplots_adjust(left=0.07, right=0.9, top=1 - 0.75 / fig.get_figheight(), bottom=0.5 / fig.get_figheight() + 0.06)
    fig.savefig(os.path.join(out, f"scan_summary{'_' + tag if tag else ''}.png"), dpi=150,
                facecolor=fig.get_facecolor())
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dirs", nargs="+")
    ap.add_argument("--out", default=".")
    ap.add_argument("--M_ref", type=float, default=1e6, help="reference mass for alpha [M_sun]")
    ap.add_argument("--tag", default="", help="read first_resonances_<tag>.csv (e.g. quad)")
    args = ap.parse_args()
    runs = load_runs(args.dirs, args.M_ref, args.tag)
    if not runs:
        print("no analyzed runs found")
        return
    os.makedirs(args.out, exist_ok=True)
    write_tables(runs, args.out, args.tag)
    plot(runs, args.out, args.tag, args.M_ref)
    for r in runs:
        f = {k: np.mean([x["outcome"] == k for x in r["rows"]]) for k in A.OUTCOMES}
        print(f"  alpha_ref = {r['alpha']:.3g} (mu = {r['mu']:.3g} eV), {len(r['rows'])} companions: "
              + ", ".join(f"{k} {v:.2f}" for k, v in f.items()))
    sfx = f"_{args.tag}" if args.tag else ""
    print(f"wrote scan_kinds{sfx}.csv, scan_transitions{sfx}.csv, scan_summary{sfx}.png to {args.out}")


if __name__ == "__main__":
    main()
