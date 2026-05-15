"""
Three diagnostic visualisations for the Nlevels convergence study.

Usage
-----
  python plot_analysis.py                 # all leaf directories
  python plot_analysis.py --dir <path>    # one directory
  python plot_analysis.py --test          # first directory only (quick check)

Output
------
For each leaf directory three figures are saved:
  <figures_root>/<bh>/<fa>/<alpha>_heatmap.jpg      -- Plot 1
  <figures_root>/<bh>/<fa>/<alpha>_convergence.jpg  -- Plot 2
  <figures_root>/<bh>/<fa>/<alpha>_relevance.jpg    -- Plot 3
"""

import os
import glob
import argparse
import time as _time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ── configuration ──────────────────────────────────────────────────────────────
OUTPUT_ROOT  = "/mnt/users/switte/Axion_SR/src/Nlevels_Runs/output"
FIGURES_ROOT = "/mnt/users/switte/Axion_SR/src/Nlevels_Runs/figures"
NMAX_LIST    = [3, 4, 5, 6, 7, 8, 15]
MIN_MAX_VAL  = 1e-30   # ignore states whose peak occupation never exceeds this


# ── I/O helpers ────────────────────────────────────────────────────────────────
def _log(msg):
    print(msg, flush=True)

def _glob1(pattern):
    m = glob.glob(pattern)
    return m[0] if m else None

def _load(path):
    if path is None:
        return None
    return pd.read_csv(path, sep=r"\s+", header=None, engine="c").values

def load_nmax(directory, nmax):
    """Return (time, spin, modes, states) or None if any file is missing."""
    t      = _load(_glob1(os.path.join(directory, f"Time_*Nmax_{nmax}.dat")))
    sp     = _load(_glob1(os.path.join(directory, f"Spin_*Nmax_{nmax}.dat")))
    modes  = _load(_glob1(os.path.join(directory, f"Modes_*Nmax_{nmax}.dat")))
    states = _load(_glob1(os.path.join(directory, f"States_*Nmax_{nmax}.dat")))
    if any(x is None for x in [t, sp, modes, states]):
        return None
    t  = t.ravel();  sp = sp.ravel()
    if modes.ndim  == 1: modes  = modes[np.newaxis, :]
    if states.ndim == 1: states = states[np.newaxis, :]
    return t, sp, modes, states

def mode_key(row):
    return (int(row[0]), int(row[1]), int(row[2]))

def mode_label(row):
    n, l, m = int(row[0]), int(row[1]), int(row[2])
    return rf"$({n},{l},{m})$"

def figure_path(directory, suffix):
    rel   = os.path.relpath(directory, OUTPUT_ROOT)
    parts = rel.split(os.sep)
    bh    = parts[0] if len(parts) > 0 else "unknown"
    fa    = parts[1] if len(parts) > 1 else "unknown"
    alpha = parts[2] if len(parts) > 2 else "plot"
    out_dir = os.path.join(FIGURES_ROOT, bh, fa)
    os.makedirs(out_dir, exist_ok=True)
    return os.path.join(out_dir, f"{alpha}_{suffix}.jpg")

def suptitle(directory):
    return os.path.relpath(directory, OUTPUT_ROOT)


# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 1 — Heatmap: peak occupation per state vs Nmax (active states only)
#
# Rows   = states (n,l,m) that are active (peak > MIN_MAX_VAL) at ≥1 Nmax,
#          sorted by (n,l,m) and grouped by principal quantum number n.
# Cols   = Nmax values present in this directory.
# Colour = log10(peak occupation); grey = state not yet in mode list at that Nmax.
# "—"    = state is in mode list but occupation below threshold.
#
# Reading guide
# -------------
# • A state turns coloured when Nmax first includes it with significant occupation.
# • Identical colour across all Nmax columns → that state is converged.
# • Horizontal group separators mark transitions between different n values.
# ═══════════════════════════════════════════════════════════════════════════════
def plot_heatmap(directory):
    data_by_nmax = {}
    for nmax in NMAX_LIST:
        res = load_nmax(directory, nmax)
        if res is not None:
            data_by_nmax[nmax] = res
    if not data_by_nmax:
        _log("  [heatmap] no data"); return

    nmaxes = sorted(data_by_nmax.keys())

    # ── find all states that are active (peak > threshold) at any Nmax ──
    active_keys = set()
    for nmax, (t, sp, modes, states) in data_by_nmax.items():
        for i, row in enumerate(modes):
            if np.max(states[i]) > MIN_MAX_VAL:
                active_keys.add(mode_key(row))

    if not active_keys:
        _log("  [heatmap] no active states"); return

    sorted_keys = sorted(active_keys)   # sorted by (n,l,m)

    # ── build peak-occupation matrix ──
    # peak[i,j] = log10(peak occ) for state i at Nmax j, or NaN
    n_states = len(sorted_keys)
    n_cols   = len(nmaxes)
    peak    = np.full((n_states, n_cols), np.nan)
    in_list = np.zeros((n_states, n_cols), dtype=bool)

    for j, nmax in enumerate(nmaxes):
        t, sp, modes, states = data_by_nmax[nmax]
        mk = {mode_key(r): idx for idx, r in enumerate(modes)}
        for i, key in enumerate(sorted_keys):
            if key in mk:
                in_list[i, j] = True
                pv = float(np.max(states[mk[key]]))
                if pv > MIN_MAX_VAL:
                    peak[i, j] = np.log10(pv)

    # ── figure (transposed: rows=Nmax, cols=states) ──
    fig_w = max(5, 0.55 * n_states + 1.5)
    fig, ax = plt.subplots(figsize=(fig_w, max(3, 1.2 * n_cols + 1.5)))

    vmin = np.nanmin(peak) if not np.all(np.isnan(peak)) else -30
    vmax = np.nanmax(peak) if not np.all(np.isnan(peak)) else 1
    vmin = min(vmin, -2);  vmax = max(vmax, 1)

    cmap = plt.cm.turbo.copy()
    cmap.set_bad("0.88")   # grey for "not in list"

    # peak has shape (n_states, n_cols); transpose so rows=Nmax, cols=states
    im = ax.imshow(peak.T, aspect="auto", origin="upper",
                   cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest")

    # mark states in list but below threshold (indices are now transposed)
    for i in range(n_states):
        for j in range(n_cols):
            if in_list[i, j] and np.isnan(peak[i, j]):
                ax.text(i, j, "–", ha="center", va="center",
                        fontsize=9, color="0.5")

    # ── x-axis: state labels, grouped by n ──
    ns = [k[0] for k in sorted_keys]
    col_labels = [f"$({k[0]},{k[1]},{k[2]})$" for k in sorted_keys]

    ax.set_xticks(range(n_states))
    ax.set_xticklabels(col_labels, fontsize=7, rotation=90)
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")

    # vertical separators between n-groups
    n_changes = [i for i in range(1, n_states) if ns[i] != ns[i-1]]
    for ic in n_changes:
        ax.axvline(ic - 0.5, color="white", lw=1.5)

    # ── y-axis: Nmax labels ──
    ax.set_yticks(range(n_cols))
    ax.set_yticklabels([f"$N_{{\\rm max}}={nm}$" for nm in nmaxes], fontsize=9)

    # horizontal grid between rows
    for y in np.arange(0.5, n_cols - 0.5, 1):
        ax.axhline(y, color="white", lw=0.6)

    cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.08, shrink=0.9)
    cb.set_label(r"$\log_{10}(\mathrm{peak\ occ.})$", fontsize=8)

    fig.tight_layout()
    path = figure_path(directory, "heatmap")
    fig.savefig(path, dpi=200, bbox_inches="tight", format="jpeg")
    plt.close(fig)
    _log(f"  [heatmap] saved: {path}")


# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 2 — Convergence: fractional spin deviation from the largest Nmax
#
# Single panel, log-log axes.
# Each line: |spin_N(t) - spin_ref(t)| / |spin_ref(t)|  vs time,
# where spin_ref is the largest available Nmax run.
#
# Reading guide
# -------------
# • Lines dropping toward zero with increasing Nmax → system is converging.
# • The Nmax where the line falls below your tolerance threshold is sufficient.
# ═══════════════════════════════════════════════════════════════════════════════
def plot_convergence(directory):
    data_by_nmax = {}
    for nmax in NMAX_LIST:
        res = load_nmax(directory, nmax)
        if res is not None:
            data_by_nmax[nmax] = res
    if len(data_by_nmax) < 2:
        _log("  [convergence] need ≥2 Nmax"); return

    nmaxes   = sorted(data_by_nmax.keys())
    nmax_ref = nmaxes[-1]
    t_ref, sp_ref, _, _ = data_by_nmax[nmax_ref]

    # reference: positive-time only
    mask_ref  = t_ref > 0
    t_pos_ref = t_ref[mask_ref]
    sp_pos_ref = sp_ref[mask_ref]

    colors = plt.cm.viridis(np.linspace(0, 0.85, len(nmaxes) - 1))

    fig, ax = plt.subplots(figsize=(10, 5))

    from matplotlib.lines import Line2D
    legend_handles = []

    for color, nmax in zip(colors, nmaxes[:-1]):
        t, sp, _, _ = data_by_nmax[nmax]
        mask = t > 0
        t_pos = t[mask];  sp_pos = sp[mask]

        # interpolate onto reference time grid (within overlap)
        sp_interp = np.interp(t_pos_ref, t_pos, sp_pos,
                              left=np.nan, right=np.nan)
        frac_dev = np.abs(sp_interp - sp_pos_ref) / (np.abs(sp_pos_ref) + 1e-300)

        valid = np.isfinite(frac_dev) & (frac_dev >= 1e-6)
        if valid.any():
            ax.plot(t_pos_ref[valid], frac_dev[valid],
                    color=color, lw=1.2)
        # always add a legend entry, note if below threshold
        lbl = f"$N_{{\\rm max}}={nmax}$"
        if not valid.any():
            lbl += r"  $(<10^{-6})$"
        legend_handles.append(Line2D([0], [0], color=color, lw=1.2, label=lbl))

    # threshold reference lines
    ax.axhline(1e-2, ls="--", color="gray", lw=0.8, label="1%")
    ax.axhline(1e-3, ls=":",  color="gray", lw=0.8, label="0.1%")
    legend_handles += [
        Line2D([0], [0], ls="--", color="gray", lw=0.8, label="1%"),
        Line2D([0], [0], ls=":",  color="gray", lw=0.8, label="0.1%"),
    ]

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(1e-6, 1)
    ax.set_xlabel("Time", fontsize=10)
    ax.set_ylabel(r"$|\tilde{a}_N - \tilde{a}_{\rm ref}|\,/\,|\tilde{a}_{\rm ref}|$",
                  fontsize=10)
    ax.set_title(
        f"Fractional spin deviation from $N_{{\\rm max}}={nmax_ref}$",
        fontsize=10
    )
    ax.legend(handles=legend_handles, fontsize=9, framealpha=0.7)

    fig.suptitle(suptitle(directory), fontsize=9)
    fig.tight_layout()
    path = figure_path(directory, "convergence")
    fig.savefig(path, dpi=200, bbox_inches="tight", format="jpeg")
    plt.close(fig)
    _log(f"  [convergence] saved: {path}")


# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 3 — State relevance: occupation vs time for all active states
#          (largest available Nmax only)
#
# All active states are shown on a single log-log panel, coloured by n.
# This directly answers "which states matter and when."
#
# Reading guide
# -------------
# • States that rise steeply and reach large values are dominant.
# • States that appear late or stay small are marginal.
# • States with the same n are shown in the same colour family.
# ═══════════════════════════════════════════════════════════════════════════════
def plot_relevance(directory):
    # use the largest available Nmax
    data = None
    for nmax in reversed(NMAX_LIST):
        res = load_nmax(directory, nmax)
        if res is not None:
            data = res;  nmax_used = nmax;  break
    if data is None:
        _log("  [relevance] no data"); return

    t, sp, modes, states = data
    mask = t > 0
    t_pos = t[mask]

    # collect active states
    active = [(i, mode_key(modes[i]), mode_label(modes[i]))
              for i in range(len(modes))
              if np.max(states[i]) > MIN_MAX_VAL]

    if not active:
        _log("  [relevance] no active states"); return

    # colour by n
    all_n   = sorted({k[0] for _, k, _ in active})
    cmap_n  = plt.cm.tab10
    n_color = {n: cmap_n(i % 10) for i, n in enumerate(all_n)}

    # sort active states by peak occupation (descending) for legend ordering
    active.sort(key=lambda x: np.max(states[x[0]]), reverse=True)

    fig, ax = plt.subplots(figsize=(11, 6))

    for idx, (i, key, label) in enumerate(active):
        occ  = states[i][mask]
        n    = key[0]
        lw   = 1.8 if idx == 0 else 1.0
        ax.plot(t_pos, occ, color=n_color[n], lw=lw, label=label, alpha=0.85)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(bottom=MIN_MAX_VAL)
    ax.set_xlabel("Time", fontsize=10)
    ax.set_ylabel("Occupation number", fontsize=10)
    ax.set_title(
        f"Active state occupations  ($N_{{\\rm max}}={nmax_used}$)\n"
        r"Coloured by principal quantum number $n$",
        fontsize=10
    )

    # legend: states ordered by peak occupation
    ax.legend(fontsize=8, ncol=2, framealpha=0.7,
              loc="upper left", handlelength=1.5)

    # secondary legend for n colours
    from matplotlib.lines import Line2D
    n_handles = [Line2D([0], [0], color=n_color[n], lw=2, label=f"$n={n}$")
                 for n in all_n]
    leg2 = ax.legend(handles=n_handles, fontsize=8, loc="lower right",
                     title="Principal $n$", title_fontsize=8, framealpha=0.7)
    ax.add_artist(leg2)
    # re-add state legend (was overwritten by leg2)
    ax.legend(fontsize=8, ncol=2, framealpha=0.7,
              loc="upper left", handlelength=1.5)

    fig.suptitle(suptitle(directory), fontsize=9)
    fig.tight_layout()
    path = figure_path(directory, "relevance")
    fig.savefig(path, dpi=200, bbox_inches="tight", format="jpeg")
    plt.close(fig)
    _log(f"  [relevance] saved: {path}")


# ═══════════════════════════════════════════════════════════════════════════════
# main
# ═══════════════════════════════════════════════════════════════════════════════
def find_leaf_dirs(root):
    leaves = []
    for dirpath, _, filenames in os.walk(root):
        if any(f.endswith(".dat") for f in filenames):
            leaves.append(dirpath)
    return sorted(leaves)


def process_directory(d):
    t0 = _time.time()
    try:
        plot_heatmap(d)
        plot_convergence(d)
        plot_relevance(d)
        _log(f"  Total: {_time.time()-t0:.1f}s")
    except Exception as e:
        import traceback
        _log(f"  ERROR in {d}: {e}")
        traceback.print_exc()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir",  default=None)
    parser.add_argument("--test", action="store_true")
    args = parser.parse_args()

    if args.dir:
        leaves = [args.dir]
    else:
        leaves = find_leaf_dirs(OUTPUT_ROOT)
        if args.test:
            leaves = leaves[:1]

    _log(f"Processing {len(leaves)} director{'y' if len(leaves)==1 else 'ies'}.")
    for i, d in enumerate(leaves, 1):
        _log(f"\n[{i}/{len(leaves)}] {d}")
        process_directory(d)


if __name__ == "__main__":
    main()
