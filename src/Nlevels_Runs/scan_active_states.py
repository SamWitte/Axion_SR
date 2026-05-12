#!/usr/bin/env python3
"""
For each BH mass, iterate over alpha values. At each alpha, use only the
largest-Nmax run for each fa. Report:
  1. All states where occupation number > THRESHOLD at any point / any fa
  2. For those states, which fa values they appear in

Processing is parallelised over (bh, alpha, fa) tasks via multiprocessing.
Results are collected first then printed in a clean, sorted layout.
"""

import os
import re
import json
import glob
import time as _time
import numpy as np
import pandas as pd
import multiprocessing
from collections import defaultdict

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
THRESHOLD  = 1e-30
BH_MASSES  = ["BH_1e4", "BH_1e8"]
N_WORKERS  = 8


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def parse_nmax(filename):
    m = re.search(r"Nmax_(\d+)", filename)
    return int(m.group(1)) if m else -1

def alpha_sort_key(name):
    m = re.search(r"alpha_(\d+\.?\d*)", name)
    return float(m.group(1)) if m else 0.0

def fa_sort_key(name):
    m = re.search(r"fa_1e(\d+)", name)
    return int(m.group(1)) if m else 0

def find_largest_nmax_files(alpha_dir):
    """Return (nmax, states_path, modes_path) for the largest Nmax in alpha_dir."""
    states_files = [
        f for f in os.listdir(alpha_dir)
        if f.startswith("States_") and f.endswith(".dat")
    ]
    if not states_files:
        return None
    states_files.sort(key=parse_nmax)
    best       = states_files[-1]
    modes_file = best.replace("States_", "Modes_")
    states_path = os.path.join(alpha_dir, best)
    modes_path  = os.path.join(alpha_dir, modes_file)
    if not os.path.exists(modes_path):
        return None
    return parse_nmax(best), states_path, modes_path

def read_modes(modes_path):
    modes = []
    with open(modes_path) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 3:
                modes.append((int(parts[0]), int(parts[1]), int(parts[2])))
    return modes

def max_occ_per_state(states_path, n_states):
    """pandas C-parser chunk reader — same method as plot_runs.py."""
    maxes = np.zeros(n_states)
    row_i = 0
    for chunk in pd.read_csv(
        states_path, sep=r"\s+", header=None, engine="c", chunksize=10
    ):
        arr     = chunk.values
        row_max = arr.max(axis=1)
        for j in range(len(arr)):
            if row_i + j < n_states:
                maxes[row_i + j] = row_max[j]
        row_i += len(arr)
        if row_i >= n_states:
            break
    return maxes


# ---------------------------------------------------------------------------
# worker (must be top-level for multiprocessing pickling)
# ---------------------------------------------------------------------------

def process_task(args):
    """
    Process one (bh, alpha, fa) combination.
    Returns (bh, alpha, fa, nmax, active_states) where active_states is a
    list of (n, l, m) tuples with max occupation > THRESHOLD.
    """
    bh, alpha, fa, alpha_dir = args

    result = find_largest_nmax_files(alpha_dir)
    if result is None:
        print(f"  [SKIP] {bh} / {alpha} / {fa} — no valid files", flush=True)
        return (bh, alpha, fa, None, [])

    nmax, states_path, modes_path = result
    modes    = read_modes(modes_path)
    n_states = len(modes)
    size_gb  = os.path.getsize(states_path) / 1e9

    t0 = _time.time()
    print(f"  [START] {bh} / {alpha} / {fa}  Nmax={nmax}  {n_states} states  {size_gb:.1f} GB", flush=True)

    maxes        = max_occ_per_state(states_path, n_states)
    active_states = [modes[i] for i in range(n_states) if maxes[i] > THRESHOLD]

    print(f"  [DONE]  {bh} / {alpha} / {fa}  ({_time.time()-t0:.0f}s)  {len(active_states)} active", flush=True)
    return (bh, alpha, fa, nmax, active_states)


# ---------------------------------------------------------------------------
# output formatting
# ---------------------------------------------------------------------------

def print_results(all_results, bh_order, fa_order, alpha_order):
    """Print collected results in a clean, sorted layout."""

    # nested dict: bh -> alpha -> fa -> list of active (n,l,m)
    data = defaultdict(lambda: defaultdict(dict))
    for bh, alpha, fa, nmax, active in all_results:
        if nmax is not None:
            data[bh][alpha][fa] = active

    print("\n" + "=" * 70)
    print("RESULTS")
    print("=" * 70)

    for bh in bh_order:
        if bh not in data:
            continue
        print(f"\n{'='*70}")
        print(f"  BH mass: {bh}")
        print(f"{'='*70}")

        for alpha in alpha_order[bh]:
            if alpha not in data[bh]:
                continue
            fa_data = data[bh][alpha]

            # aggregate: state -> sorted list of fa values where it is active
            state_to_fas = defaultdict(list)
            for fa in fa_order[bh]:
                for state in fa_data.get(fa, []):
                    state_to_fas[state].append(fa)

            print(f"\n  {alpha}")
            print(f"  {'-'*40}")

            if not state_to_fas:
                print(f"    (no states with occupation > {THRESHOLD:.0e})")
                continue

            # determine which fa values exist for this alpha
            fas_present = [fa for fa in fa_order[bh] if fa in fa_data]

            # header row
            col_w = 16
            header = f"    {'state (n,l,m)':<16}" + "".join(f"{fa:>{col_w}}" for fa in fas_present)
            print(header)
            print("    " + "-" * (16 + col_w * len(fas_present)))

            for state in sorted(state_to_fas):
                n, l, m = state
                label   = f"({n},{l},{m})"
                active_fas = set(state_to_fas[state])
                row = f"    {label:<16}" + "".join(
                    f"{'yes':>{col_w}}" if fa in active_fas else f"{'':>{col_w}}"
                    for fa in fas_present
                )
                print(row)

    print()


# ---------------------------------------------------------------------------
# single-task mode: process one alpha_dir, write JSON result
# ---------------------------------------------------------------------------

def run_single_task(alpha_dir, output_json):
    """Process one (bh, fa, alpha) directory and write JSON output."""
    parts = alpha_dir.rstrip("/").split(os.sep)
    # expect .../output/BH_xxx/fa_xxx/alpha_xxx
    alpha = parts[-1]
    fa    = parts[-2]
    bh    = parts[-3]

    bh_dummy, alpha_dummy, fa_dummy, _dir = bh, alpha, fa, alpha_dir
    _, _, _, nmax, active_states = process_task((bh_dummy, alpha_dummy, fa_dummy, alpha_dir))

    result = {
        "bh":           bh,
        "alpha":        alpha,
        "fa":           fa,
        "nmax":         nmax,
        "active_states": [list(s) for s in active_states],
    }
    os.makedirs(os.path.dirname(os.path.abspath(output_json)), exist_ok=True)
    with open(output_json, "w") as fh:
        json.dump(result, fh)
    print(f"Wrote {output_json}", flush=True)


# ---------------------------------------------------------------------------
# aggregate mode: read all JSON files and print combined table
# ---------------------------------------------------------------------------

def run_aggregate(results_dir):
    """Read all JSON results and print the combined active-states table."""
    json_files = glob.glob(os.path.join(results_dir, "*.json"))
    if not json_files:
        print(f"No JSON files found in {results_dir}")
        return

    all_results = []
    for fp in json_files:
        with open(fp) as fh:
            r = json.load(fh)
        all_results.append((
            r["bh"], r["alpha"], r["fa"], r["nmax"],
            [tuple(s) for s in r["active_states"]],
        ))

    bhs   = sorted({r[0] for r in all_results})
    fas   = {bh: sorted({r[2] for r in all_results if r[0] == bh}, key=fa_sort_key)   for bh in bhs}
    alphs = {bh: sorted({r[1] for r in all_results if r[0] == bh}, key=alpha_sort_key) for bh in bhs}
    print_results(all_results, bhs, fas, alphs)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--alpha-dir",    default=None,
                    help="Process a single alpha directory (requires --output-json)")
    ap.add_argument("--output-json",  default=None,
                    help="Path to write JSON result for --alpha-dir mode")
    ap.add_argument("--aggregate",    default=None,
                    help="Directory of JSON results to aggregate and print")
    args = ap.parse_args()

    if args.alpha_dir:
        if not args.output_json:
            raise SystemExit("--output-json is required with --alpha-dir")
        run_single_task(args.alpha_dir, args.output_json)
        return

    if args.aggregate:
        run_aggregate(args.aggregate)
        return

    # Default: full scan across all BH masses (original behaviour)
    tasks      = []
    bh_order   = []
    fa_order   = {}
    alpha_order = {}

    for bh in BH_MASSES:
        bh_dir = os.path.join(OUTPUT_DIR, bh)
        if not os.path.isdir(bh_dir):
            continue
        bh_order.append(bh)

        fa_dirs = sorted(
            [d for d in os.listdir(bh_dir) if os.path.isdir(os.path.join(bh_dir, d))],
            key=fa_sort_key,
        )
        fa_order[bh] = fa_dirs

        all_alphas = set()
        for fa in fa_dirs:
            fa_path = os.path.join(bh_dir, fa)
            for name in os.listdir(fa_path):
                if name.startswith("alpha_") and os.path.isdir(os.path.join(fa_path, name)):
                    all_alphas.add(name)
        alphas = sorted(all_alphas, key=alpha_sort_key)
        alpha_order[bh] = alphas

        for alpha in alphas:
            for fa in fa_dirs:
                alpha_dir = os.path.join(bh_dir, fa, alpha)
                if os.path.isdir(alpha_dir):
                    tasks.append((bh, alpha, fa, alpha_dir))

    n_tasks  = len(tasks)
    workers  = min(N_WORKERS, n_tasks)
    print(f"Found {n_tasks} tasks — running with {workers} workers\n", flush=True)

    t0 = _time.time()
    with multiprocessing.Pool(workers) as pool:
        all_results = pool.map(process_task, tasks)
    print(f"\nAll tasks done in {_time.time()-t0:.0f}s", flush=True)

    print_results(all_results, bh_order, fa_order, alpha_order)


if __name__ == "__main__":
    main()
