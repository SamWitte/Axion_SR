#!/usr/bin/env python3
"""Consistency of run_Nlevels output directories: Time, Spin, MassBH have the same length, every States line has that
many columns, States and Modes have the same number of lines, no NUL bytes, final time 1e9 yr.
usage: check_run_files.py <run directory> [...]"""
import glob, os, sys
import numpy as np
bad = 0; n = 0
for d in sys.argv[1:]:
    g = lambda k: (glob.glob(os.path.join(d, k + "_*.dat")) or [None])[0]
    if g("Time") is None: print("NO OUTPUT", d); bad += 1; continue
    n += 1; msg = []
    t = np.loadtxt(g("Time")); a = np.loadtxt(g("Spin")); M = np.loadtxt(g("MassBH"))
    if not (len(t) == len(a) == len(M)): msg.append(f"lengths {len(t)} {len(a)} {len(M)}")
    if abs(t[-1] - 1e9) > 1e3: msg.append(f"t_end {t[-1]:.4g}")
    if not (np.all(np.isfinite(a)) and np.all(np.isfinite(M)) and np.all(np.diff(t) >= 0)): msg.append("non-finite or non-monotonic")
    nm = sum(1 for _ in open(g("Modes"))); ns = 0
    with open(g("States"), "rb") as f:
        for line in f:
            ns += 1
            if b"\0" in line: msg.append(f"NUL in States line {ns}")
            if line.count(b"\t") + 1 != len(t): msg.append(f"States line {ns}: {line.count(b'\t') + 1} columns, {len(t)} times")
    if ns != nm: msg.append(f"States {ns} lines, Modes {nm}")
    for k in ("Time", "Spin", "MassBH", "Modes"):
        if b"\0" in open(g(k), "rb").read(): msg.append("NUL in " + k)
    if glob.glob(os.path.join(d, "*.reducing")): msg.append("leftover .reducing file")
    if msg: bad += 1; print("BAD", d, "; ".join(msg[:4]))
print(f"{n} run directories checked, {bad} with a problem")
