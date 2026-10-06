#!/usr/bin/env python3
"""
AGN central BH + axion cloud (Julia solver) -> companions -> first resonance.

    1. python agn_julia_pipeline.py export  --rundir runs --N 20 [--seed 7] [--mu 3e-17] [--logM_today 5.5 6.5]
         Params (PARAMS below, with any overrides) -> runs/params.json, runs/hist_XXXX.{json,txt}
    2. julia evolve_agn_history.jl --rundir runs --mu <mu_eV> [--f_a 1e18] [--Nmax 4] [--ids ...]
         (prints the exact command after step 1; run histories in parallel on a
         server by splitting --ids)
    3. python agn_julia_pipeline.py analyze --rundir runs [--n_inj 10]
         companions dropped into each evolved history (run_external), first
         resonance per companion -> runs/first_resonances.csv, crossings.csv,
         runs/first_resonance_summary.png, runs/gallery.png, runs/cloud_XXXX.png

Step 3 always reloads runs/params.json, so the boson mass, occupation threshold,
injection and disk settings are those the histories were exported with.
"""
import argparse
import csv
import dataclasses
import json
import os

import numpy as np

import agn_cloud_resonances as A

# Edit for a new run (only read by "export"; "analyze" uses runs/params.json).
PARAMS = A.Params(mu_eV=3e-17, l_stars=(2, 3), logM_today=(5.5, 6.5), z_obs_max=1.0,
                  occ_threshold=1e-6)


def save_params(p, path):
    with open(path, "w") as fh:
        json.dump(dataclasses.asdict(p), fh, indent=1)


def load_params(path):
    d = json.load(open(path))
    as_tuple = lambda v: tuple(as_tuple(x) for x in v) if isinstance(v, list) else v
    return A.Params(**{k: as_tuple(v) for k, v in d.items()})


def evolution_track(rundir, idx, z_seed):
    """evol_XXXX.csv -> track dict for A.plot_cloud_track."""
    data = np.genfromtxt(os.path.join(rundir, f"evol_{idx:04d}.csv"), delimiter=",", names=True)
    cols = [c for c in data.dtype.names if c.startswith("Mc_")]
    M = data["M_Msun"] * A.Msun
    return dict(t=data["age_yr"] * A.yr + A.t_of_z(z_seed), M=M, a=data["a"],
                Mc=np.array([data[c] * M for c in cols]),
                levels=[A._parse_state(c) for c in cols])


def plot_evolution(track, fname, show_floor=1e-8):
    """M, a and Mc/M vs time; only levels whose Mc/M ever exceeds show_floor."""
    import matplotlib.pyplot as plt
    t = track["t"] / A.Gyr
    frac = track["Mc"] / track["M"]
    keep = [j for j in range(len(track["levels"])) if frac[j].max() > show_floor]
    cols = plt.cm.tab10(np.arange(len(keep)) % 10)
    fig, ax = plt.subplots(3, 1, figsize=(9, 7.5), sharex=True, gridspec_kw=dict(height_ratios=[1, 1, 1.4]))
    ax[0].semilogy(t, track["M"] / A.Msun, color="k")
    ax[0].set_ylabel(r"$M\ [M_\odot]$")
    ax[1].plot(t, track["a"], color="k")
    ax[1].set(ylabel="spin a", ylim=(0, 1))
    for c, j in zip(cols, keep):
        ax[2].semilogy(t, np.maximum(frac[j], 1e-30), color=c, label=A.ket(track["levels"][j]))
    ax[2].set(ylabel=r"$M_c/M$", ylim=(show_floor, 1), xlabel="cosmic time [Gyr]")
    if keep:
        ax[2].legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1))
    fig.tight_layout()
    fig.savefig(fname, dpi=130)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["export", "analyze"])
    ap.add_argument("--rundir", default="runs")
    ap.add_argument("--N", type=int, default=20, help="export: number of histories")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--mu", type=float, default=None, help="export: override PARAMS.mu_eV")
    ap.add_argument("--logM_today", type=float, nargs=2, default=None, help="export: override PARAMS.logM_today")
    ap.add_argument("--n_inj", type=int, default=10, help="analyze: companions per history")
    ap.add_argument("--n_panels", type=int, default=6)
    args = ap.parse_args()
    pfile = os.path.join(args.rundir, "params.json")

    if args.step == "export":
        p = PARAMS
        if args.mu is not None:
            p = dataclasses.replace(p, mu_eV=args.mu)
        if args.logM_today is not None:
            p = dataclasses.replace(p, logM_today=tuple(args.logM_today))
        os.makedirs(args.rundir, exist_ok=True)
        save_params(p, pfile)
        A.export_histories(args.N, p, args.rundir, seed=args.seed)
        print(f"wrote {args.N} histories and params.json to {args.rundir}. Next:\n"
              f"  julia evolve_agn_history.jl --rundir {args.rundir} --mu {p.mu_eV!r} --f_a 1e18 --Nmax 4")
        return

    p = load_params(pfile)
    res = A.run_external(args.rundir, p, n_inj=args.n_inj)
    if not res:
        print("no evolved histories found (evol_XXXX.csv)")
        return
    kinds = ["hyperfine", "fine", "stalled", "none_ahead", "no_SR_state"]
    print(f"{len(res)} companions in {len(set(R['id'] for R in res))} histories")
    for k in kinds:
        print(f"  {k:12s} {np.mean([R['outcome'] == k for R in res]):.3f}")
    with open(os.path.join(args.rundir, "first_resonances.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["history", "injection", "age_inj_Gyr", "M_Msun", "a", "alpha", "sigma", "r0_rg",
                    "populated", "outcome", "first_init", "first_final", "first_f_Hz", "t_cross_yr"])
        for R in res:
            F = R["first"]
            w.writerow([R["id"], R["inj"], R["age_inj"] / A.Gyr, R["M"] / A.Msun, R["a"], R["alpha"],
                        R["sigma"], R["r0"], " ".join(A.ket_ascii(s) for s in R["live_states"]), R["outcome"],
                        A.ket_ascii(F["init"]) if F else "", A.ket_ascii(F["final"]) if F else "",
                        F["f_res"] if F else "", F["t_cross"] / A.yr if F else ""])
    rows = A.crossing_table(res)
    if rows:
        with open(os.path.join(args.rundir, "crossings.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
    A.plot_baseline(res, p, os.path.join(args.rundir, "first_resonance_summary.png"))
    A.plot_gallery(res, p, os.path.join(args.rundir, "gallery.png"), n_panels=args.n_panels)
    for idx in sorted(set(R["id"] for R in res)):
        d = json.load(open(os.path.join(args.rundir, f"hist_{idx:04d}.json")))
        plot_evolution(evolution_track(args.rundir, idx, d["z_seed"]),
                       os.path.join(args.rundir, f"cloud_{idx:04d}.png"))
    print(f"wrote first_resonances.csv, crossings.csv and plots to {args.rundir}")


if __name__ == "__main__":
    main()
