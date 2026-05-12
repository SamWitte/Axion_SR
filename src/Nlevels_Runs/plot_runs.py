"""
Plot spin and occupation numbers for each parameter combination in the output directory.
One figure per leaf directory (BH / fa / alpha), saved as JPEG at dpi=200.

Layout (shared x-axis = time):
  Panel 0 : spin  ã  in [0,1]  (linear)
  Panels 1-7 : States for Nmax = 3, 4, 5, 6, 7, 8, 15  (log, shared y limits)
"""

import os
import sys
import glob
import gc
import time as _time
import zipfile
import multiprocessing
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MAX_PLOT_POINTS = 50000   # downsample to this many points before passing to matplotlib
N_WORKERS = 40


def _log(msg):
    print(msg, flush=True)

_HERE        = os.path.dirname(os.path.abspath(__file__))
OUTPUT_ROOT  = os.path.join(_HERE, "output")
FIGURES_ROOT = os.path.join(_HERE, "figures")
NMAX_LIST    = [4, 5, 6, 7, 8, 15]
MIN_MAX_VAL = 1e-30   # only plot states whose max exceeds this


def find_leaf_dirs(root):
    """Return all directories that contain .dat files (leaf dirs)."""
    leaves = []
    for dirpath, dirnames, filenames in os.walk(root):
        if any(f.endswith(".dat") for f in filenames):
            leaves.append(dirpath)
    return sorted(leaves)


def glob_one(pattern):
    matches = glob.glob(pattern)
    if not matches:
        return None
    return matches[0]


def _fast_load(path):
    """Load a whitespace/tab-delimited .dat file quickly via pandas."""
    return pd.read_csv(path, sep=r"\s+", header=None, engine="c").values


def load_time(directory, nmax):
    path = glob_one(os.path.join(directory, f"Time_*Nmax_{nmax}.dat"))
    if path is None:
        return None
    return _fast_load(path).ravel()


def load_spin(directory, nmax):
    path = glob_one(os.path.join(directory, f"Spin_*Nmax_{nmax}.dat"))
    if path is None:
        return None
    return _fast_load(path).ravel()


def load_modes(directory, nmax):
    """Return array of shape (n_modes, 4): columns n, l, m, freq."""
    path = glob_one(os.path.join(directory, f"Modes_*Nmax_{nmax}.dat"))
    if path is None:
        return None
    return _fast_load(path)


def load_states_active(directory, nmax, mask, stride, min_val):
    """Read States file in chunks via the C parser (much faster than line-by-line).

    Returns (active_row_indices, data) where data has shape
    (n_active, n_plot_points) — already masked to positive time and
    downsampled.
    """
    path = glob_one(os.path.join(directory, f"States_*Nmax_{nmax}.dat"))
    if path is None:
        return np.array([], dtype=int), np.zeros((0, 0))
    size_mb = os.path.getsize(path) / 1e6
    _log(f"    Streaming States Nmax={nmax} ({size_mb:.0f} MB) ...")
    t0 = _time.time()
    active_idx, active_data = [], []
    row_i = 0
    for chunk in pd.read_csv(path, sep=r"\s+", header=None, engine="c", chunksize=10):
        arr = chunk.values
        row_max = arr.max(axis=1)
        for j in range(len(arr)):
            if row_max[j] > min_val:
                active_idx.append(row_i + j)
                active_data.append(arr[j][mask][::stride])
        row_i += len(arr)
    n_active = len(active_idx)
    _log(f"    ... done in {_time.time()-t0:.1f}s  ({n_active} active modes)")
    if n_active == 0:
        return np.array([], dtype=int), np.zeros((0, 0))
    return np.array(active_idx, dtype=int), np.array(active_data)


def make_label(mode_row):
    """Format quantum numbers (n, l, m) as a label string."""
    n, l, m = int(mode_row[0]), int(mode_row[1]), int(mode_row[2])
    return rf"$({n},{l},{m})$"


def _compute_xlim(directory):
    """Scan only Time files to compute shared x-axis limits without loading States."""
    t_min_val, t_max_val = None, None
    for nmax in NMAX_LIST:
        t = load_time(directory, nmax)
        if t is None or len(t) < 2:
            continue
        pos = t[t > 0]
        if len(pos) == 0:
            continue
        t_min_val = pos.min() if t_min_val is None else min(t_min_val, pos.min())
        t_max_val = t[-1]    if t_max_val is None else max(t_max_val, t[-1])
    return t_min_val, t_max_val


def plot_directory(directory):
    _log(f"  Loading Time/Spin ...")
    # Load time from first available Nmax as reference for the spin panel
    time = None
    for nmax in NMAX_LIST:
        time = load_time(directory, nmax)
        if time is not None and len(time) > 1:
            break
    if time is None:
        _log(f"  Skipping {directory}: no Time file found")
        return

    spin = None
    for nmax in NMAX_LIST:
        spin = load_spin(directory, nmax)
        if spin is not None and len(spin) > 1:
            break
    if spin is None:
        _log(f"  Skipping {directory}: no Spin file found")
        return

    # Compute x limits from Time files only (avoids loading all States up front)
    t_min, t_max = _compute_xlim(directory)
    if t_min is None:
        t_min = time[time > 0].min()
        t_max = time[-1]

    n_panels = 1 + len(NMAX_LIST)   # spin + one per Nmax
    height_ratios = [1.5] + [2] * len(NMAX_LIST)

    fig, axes = plt.subplots(
        n_panels, 1,
        figsize=(14, 3 * n_panels),
        sharex=True,
        gridspec_kw={"height_ratios": height_ratios},
    )
    fig.subplots_adjust(hspace=0.05)

    # ---- Panel 0: Spin ----
    ax0 = axes[0]
    ax0.set_xscale("log")
    mask = time > 0
    t_spin = time[mask]
    sp_spin = spin[mask]
    stride0 = max(1, len(t_spin) // MAX_PLOT_POINTS)
    ax0.plot(t_spin[::stride0], sp_spin[::stride0], color="tab:blue", lw=1.2, label=r"$\tilde{a}$")
    del t_spin, sp_spin, time, spin
    ax0.set_ylim(0, 1)
    ax0.set_ylabel(r"$\tilde{a}$")
    ax0.legend(fontsize=8, loc="upper right", framealpha=0.6)
    ax0.tick_params(labelbottom=False)

    # ---- Panels 1-7: States per Nmax ----
    # Load, plot, and immediately free each Nmax to avoid holding all in memory
    occ_ymin, occ_ymax = 1e-16, 10.0
    prop_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    for idx, nmax in enumerate(NMAX_LIST):
        _log(f"  Processing Nmax={nmax} ...")
        ax = axes[1 + idx]
        ax.set_yscale("log")
        ax.set_ylim(occ_ymin, occ_ymax)
        ax.set_ylabel(f"$N_{{\\rm max}}={nmax}$\nOccupation", fontsize=8)

        t = load_time(directory, nmax)
        modes = load_modes(directory, nmax)

        if t is None or modes is None:
            ax.text(0.5, 0.5, "no data", transform=ax.transAxes,
                    ha="center", va="center", color="gray")
        else:
            if modes.ndim == 1:
                modes = modes[np.newaxis, :]
            mask = t > 0
            t_plot = t[mask]
            stride = max(1, len(t_plot) // MAX_PLOT_POINTS)

            active_idx, active_data = load_states_active(
                directory, nmax, mask, stride, MIN_MAX_VAL)

            if len(active_idx) == 0:
                ax.text(0.5, 0.5, "no data", transform=ax.transAxes,
                        ha="center", va="center", color="gray")
            else:
                t_ds = t_plot[::stride]
                for plotted, (row_i, row) in enumerate(zip(active_idx, active_data)):
                    color = prop_cycle[plotted % len(prop_cycle)]
                    label = make_label(modes[row_i]) if row_i < len(modes) else f"mode {row_i}"
                    ax.plot(t_ds, row, lw=0.8, color=color, label=label)
                ax.legend(fontsize=7, loc="upper right", ncol=2,
                          framealpha=0.6, handlelength=1.2)

            del t, t_plot, modes, active_idx, active_data
            gc.collect()

        if idx < len(NMAX_LIST) - 1:
            ax.tick_params(labelbottom=False)

    # Apply shared x limits to all panels
    ax0.set_xlim(t_min, t_max)

    # x-axis label on bottom panel
    axes[-1].set_xlabel("Time")

    # Title from directory path
    rel = os.path.relpath(directory, OUTPUT_ROOT)
    fig.suptitle(rel, fontsize=10, y=1.002)

    # Build output path: figures/<bh_mass>/<fa>/<alpha>.jpg
    rel = os.path.relpath(directory, OUTPUT_ROOT)
    parts = rel.split(os.sep)
    if len(parts) >= 3:
        bh_label, fa_label, alpha_label = parts[0], parts[1], parts[2]
    elif len(parts) == 2:
        bh_label, fa_label, alpha_label = parts[0], parts[1], "plot"
    else:
        bh_label, fa_label, alpha_label = parts[0], "unknown_fa", "plot"

    out_dir = os.path.join(FIGURES_ROOT, bh_label, fa_label)
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, alpha_label + ".jpg")

    _log(f"  Saving figure ...")
    fig.savefig(out_path, dpi=200, bbox_inches="tight", format="jpeg")
    plt.close(fig)
    gc.collect()
    _log(f"  [DONE] Saved: {out_path}")


def zip_figures():
    zip_path = FIGURES_ROOT + ".zip"
    _log(f"\nZipping {FIGURES_ROOT} -> {zip_path} ...")
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for dirpath, _, filenames in os.walk(FIGURES_ROOT):
            for fname in filenames:
                fpath = os.path.join(dirpath, fname)
                arcname = os.path.relpath(fpath, os.path.dirname(FIGURES_ROOT))
                zf.write(fpath, arcname)
    _log(f"[DONE] Archive saved: {zip_path}")


def _process_one(args):
    i, total, d = args
    _log(f"\n[{i}/{total}] Processing: {d}")
    t0 = _time.time()
    try:
        plot_directory(d)
        _log(f"  Total time for this dir: {_time.time()-t0:.1f}s")
    except Exception as e:
        import traceback
        _log(f"  ERROR in {d}: {e}")
        traceback.print_exc()


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--directory", default=None,
                    help="Process a single leaf directory (skips zip)")
    ap.add_argument("--no-zip", action="store_true",
                    help="Skip zipping the figures directory")
    args = ap.parse_args()

    if args.directory:
        d = os.path.abspath(args.directory)
        _log(f"Single-directory mode: {d}")
        _process_one((1, 1, d))
        return

    leaves = find_leaf_dirs(OUTPUT_ROOT)
    total = len(leaves)
    _log(f"Found {total} directories to process.")
    job_args = [(i, total, d) for i, d in enumerate(leaves, 1)]
    workers = min(N_WORKERS, total) if total > 0 else 1
    _log(f"Using {workers} parallel workers.")
    with multiprocessing.Pool(workers) as pool:
        pool.map(_process_one, job_args)

    if not args.no_zip:
        zip_figures()


if __name__ == "__main__":
    main()
