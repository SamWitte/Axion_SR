#!/usr/bin/env python3
"""Task chunks of the recomputation of all "to infinity" rate points with the corrected outgoing condition
(gf_radial inf_wkb, see ../allL for the first pass with all multipoles and the old condition).
Same 64,362 tasks as allL/tasks.  The channels whose final state gave hour-long points in allL (near threshold, small
alpha) are split into one task per alpha for alpha <= 0.125 and go first, in small shuffled chunks, so that a slow point
does not hold up a whole chunk."""
import glob, os, random
A = "/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext"; W = A + "/allW"
ALPHAS = open(A + "/sgx_alphas.txt").read().strip().split(",")
SLOW3 = {"11-10-10", "12-10-10", "13-12-8", "13-12-12", "12-8-8", "12-6-6"}
lines = [l.strip() for f in sorted(glob.glob(A + "/allL/tasks/chunk_*.txt")) for l in open(f) if l.strip()]
slow, rest = [], []
for l in lines:
    p = l.split()
    if p[2] in SLOW3 and not any(x.startswith("alphas=") for x in p):
        lo = [a for a in ALPHAS if float(a) <= 0.125]; hi = [a for a in ALPHAS if float(a) > 0.125]
        slow += [f"{l} alphas={a}" for a in lo] + [f"{l} alphas={','.join(hi)}"]
    else: rest.append(l)
random.seed(20261119); random.shuffle(slow)
chunks = [slow[i:i + 8] for i in range(0, len(slow), 8)] + [rest[i:i + 40] for i in range(0, len(rest), 40)]
for i, c in enumerate(chunks): open(f"{W}/tasks/chunk_{i:05d}.txt", "w").write("\n".join(c) + "\n")
print(len(lines), "channels ->", len(slow), "single-alpha tasks of", len(slow) // 9, "slow-family channels +", len(rest), "channel tasks;", len(chunks), "chunks")
