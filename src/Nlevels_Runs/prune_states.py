#!/usr/bin/env python3
"""Drop inactive modes (States rows with max occupation <= threshold)."""

import glob
import os
import time

import numpy as np

MIN_MAX_VAL = 1e-30


def _log(msg):
    print(msg, flush=True)


def _glob_one(directory, prefix, nmax):
    matches = glob.glob(os.path.join(directory, f"{prefix}_*Nmax_{nmax}.dat"))
    return matches[0] if matches else None


def _line_max(line, min_val=MIN_MAX_VAL):
    # numpy parses a ~10-40 MB mode-row far faster than Python split/float.
    arr = np.fromstring(line, dtype=np.float64, sep=" ")
    if arr.size == 0:
        return 0.0
    return float(arr.max())


def prune_states_modes(directory, nmax, min_val=MIN_MAX_VAL):
    """Keep only States rows (modes) whose peak occupation exceeds min_val."""
    states_path = _glob_one(directory, "States", nmax)
    modes_path = _glob_one(directory, "Modes", nmax)
    if states_path is None:
        return 0, 0

    marker = f"{states_path}.pruned"
    if os.path.isfile(marker):
        return 0, 0

    if os.path.getsize(states_path) == 0:
        open(marker, "w").close()
        return 0, 0

    size_before = os.path.getsize(states_path)
    t0 = time.time()
    tmp_path = f"{states_path}.prtmp"
    kept_idx = []
    row_i = 0
    n_in = 0
    n_out = 0

    with open(states_path, "r") as fin, open(tmp_path, "w") as fout:
        for line in fin:
            if not line.strip():
                continue
            n_in += 1
            if _line_max(line, min_val) > min_val:
                kept_idx.append(row_i)
                fout.write(line if line.endswith("\n") else line + "\n")
                n_out += 1
            row_i += 1

    os.replace(tmp_path, states_path)
    open(marker, "w").close()

    if modes_path is not None:
        if kept_idx:
            kept = set(kept_idx)
            mode_lines = []
            with open(modes_path, "r") as fin:
                for i, line in enumerate(fin):
                    if i in kept:
                        mode_lines.append(line if line.endswith("\n") else line + "\n")
            with open(modes_path, "w") as fout:
                fout.writelines(mode_lines)
        else:
            open(modes_path, "w").close()

    size_after = os.path.getsize(states_path)
    _log(
        f"    pruned Nmax={nmax}: {n_in} -> {n_out} modes, "
        f"{size_before / 1e6:.0f} -> {size_after / 1e6:.1f} MB "
        f"({time.time() - t0:.1f}s)"
    )
    return n_in, n_out
