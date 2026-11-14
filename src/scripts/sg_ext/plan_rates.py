#!/usr/bin/env python3
"""Plan the rate computation of the extended self-gravity system (2026-11-12).

States   S = Table I of the draft (M = 10 Msun column; the 1e4 column is a subset) + the 21 levels of the first SG
             campaign + (n,4,4) for n = 5..20                                             -> levels_sgx.txt
Channels A = every channel of the paper's list (rate_sve/load_rate_input_Nmax_45.txt) and of the first SG campaign
             (sg21_channels.txt) whose three levels are in S (an SI table exists for each)
         B = every rule-allowed channel (rate_sve/gen_custom_input.py rules) among the levels of S with n <= 20
             that contains an (n,4,4) level and has no SI table under either ordering/label
  -> sgx_channels.txt (the --rate_input list of the runs), tasks/chunk_p<prio>_<k>.txt (Compute_sg_points.jl lists)

Every channel is computed on the campaign's alpha grid (SI and SG together), including those that already have an SG
table from the first campaign (those are on a coarse grid and with the default fourth-quantum step).
Priority (chunks are taken in file-name order): p0 = all three levels in SG21 + (n,4,4); p1 = all three levels among
p0's and the levels populated at alpha <= 0.2 in Fig. 2 with n <= 20; p2 = the rest.

    python3 plan_rates.py                 # report
    python3 plan_rates.py --write         # write levels, channel list and task chunks
    python3 plan_rates.py --only-missing --write   # chunks for channels whose .pts is incomplete
"""
import os, re, sys, collections, argparse, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from states_lib import *

ap = argparse.ArgumentParser()
ap.add_argument("--write", action="store_true")
ap.add_argument("--only-missing", action="store_true")
ap.add_argument("--chunk-sec", type=float, default=5400.0, help="target cost of a chunk (s)")
ap.add_argument("--alphas", default="0.03,0.04,0.05,0.062,0.075,0.09,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.28")
a = ap.parse_args()

HERE = os.path.dirname(os.path.abspath(__file__))
RS = "/data/dust/user/samwitte/Axion_SR_sg/src/rate_sve/"
OUT = os.path.join(HERE, "rates_out")
ALPHAS = [float(x) for x in a.alphas.split(",")]

def st(tok): return tuple(int(x) for x in tok.split("-")) if "-" in tok else tuple(int(c) for c in tok)
def name(c): return f"{fmt(c[0])}_{fmt(c[1])}_{fmt(c[2])}_{c[3]}"
def key(c): return (frozenset([c[0], c[1]]), c[2])
def read_list(path):
    out = []
    for ln in open(path):
        p = ln.split()
        if len(p) >= 4 and not p[0].startswith("#"): out.append((st(p[0]), st(p[1]), st(p[2]), p[3]))
    return out

si, sg = {}, {}
for f in os.listdir(RS):
    m = re.match(r"(.+)_(.+)_(.+)_(BH|Inf)_LvrHc_(SG)?\.dat$", f)
    if m: (sg if m.group(5) else si)[(st(m.group(1)), st(m.group(2)), st(m.group(3)), m.group(4))] = f
si_keys = collections.defaultdict(list)
for c in si: si_keys[key(c)].append(c)

S = sorted(set(expand(T1_10)) | set(SG21) | set(N44)); Sset = set(S)
inS = lambda c: c[0] in Sset and c[1] in Sset and c[2] in Sset
A, seenA = [], set()
for c in read_list(RS + "load_rate_input_Nmax_45.txt") + read_list(RS + "sg21_channels.txt"):
    if inS(c) and key(c) not in seenA and c in si:
        A.append(c); seenA.add(key(c))
S20 = [s for s in S if s[0] <= 20]
B = [c for c in channels(S20) if any(s[1:] == (4, 4) for s in c[:3]) and key(c) not in si_keys and key(c) not in seenA]
# Once the channel list has been written, it is the reference: the merge adds SI tables for the B channels, so
# re-deriving B from "no SI table" afterwards would lose them (this made the first --only-missing report too short).
LIST = os.path.join(HERE, "sgx_channels.txt")
if os.path.isfile(LIST) and not a.write or (a.only_missing and os.path.isfile(LIST)):
    allc = read_list(LIST); ka = set(key(c) for c in A)
    A = [c for c in allc if key(c) in ka]; B = [c for c in allc if key(c) not in ka]
used = set(s for c in A + B for s in c[:3])
levels = [s for s in S if s in used]
P0 = set(SG21) | set(N44); P1 = P0 | set(s for s in FIG_A02 if s[0] <= 20)
def prio(c):
    if a.only_missing:   # leftovers: p0 = cheap (escaping fourth quantum, or every level n <= 20 and l <= 8), p1 = heavy
        return 0 if (c[3] == "Inf" or (max(s[0] for s in c[:3]) <= 20 and max(s[1] for s in c[:3]) <= 8)) else 1
    return 0 if all(s in P0 for s in c[:3]) else 1 if all(s in P1 for s in c[:3]) else 2
print(f"states {len(S)} (in channels: {len(levels)})   A: {len(A)} channels with SI table ({sum(1 for c in A if c in sg)} have an old SG table)   "
      f"B: {len(B)} new (n,4,4) channels   total {len(A) + len(B)}")

def amax(c): return min(0.95 * min(s[2] for s in c[:3]) / (2 * (1 + math.sqrt(1 - 0.95**2))) * 1.03, 1.7)
def npts_done(c):
    f = os.path.join(OUT, name(c) + ".pts")
    if not os.path.isfile(f): return 0
    return len(set(round(float(ln.split()[0]), 10) for ln in open(f) if len(ln.split()) >= 9))
def cost(c):   # seconds per channel, from the pilots: ~6 s per point, bound fourth quantum ~4/alpha s more (fine step)
    return sum((6.0 + (4.0 / x if c[3] == "BH" else 0.0)) for x in ALPHAS if x <= amax(c)) + 5.0 * sum(1 for s in c[:3] if s[0] <= 8)
tasks = []
for c in A + B:
    need = sum(1 for x in ALPHAS if x <= amax(c))
    if a.only_missing and npts_done(c) >= need: continue
    tasks.append((prio(c), c, need))
cnt = collections.Counter((p, c[3]) for p, c, _ in tasks)
tot = sum(cost(c) for _, c, _ in tasks)
print(f"tasks: {len(tasks)} channels, {sum(t[2] for t in tasks)} points, est. {tot / 3600:.0f} core-h;  by (priority, final): {dict(sorted(cnt.items()))}")
for p in (0, 1, 2):
    print(f"   p{p}: {sum(1 for q, _, _ in tasks if q == p)} channels, est. {sum(cost(c) for q, c, _ in tasks if q == p) / 3600:.0f} core-h")
if a.write:
    if not a.only_missing:
        with open(os.path.join(HERE, "levels_sgx.txt"), "w") as f:
            f.write("".join(f"{n}-{l}-{m}\n" for n, l, m in levels))
        with open(os.path.join(HERE, "sgx_channels.txt"), "w") as f:
            for c in A + B: f.write("    ".join((fmt(c[0]), fmt(c[1]), fmt(c[2]), c[3])) + " \n")
        with open(os.path.join(HERE, "sgx_alphas.txt"), "w") as f: f.write(a.alphas + "\n")
    td = os.path.join(HERE, "tasks"); os.makedirs(td, exist_ok=True)
    for old in os.listdir(td): os.remove(os.path.join(td, old))
    for old in os.listdir(os.path.join(HERE, "done")) if os.path.isdir(os.path.join(HERE, "done")) else []:
        os.remove(os.path.join(HERE, "done", old))
    nch = 0
    for p in (0, 1, 2):
        # group channels that share their two incoming levels (cache hits), fill chunks to the target cost
        ts = sorted((c for q, c, _ in tasks if q == p), key=lambda c: (c[0], c[1], c[2]))
        cur, cc, k = [], 0.0, 0
        def flush():
            global nch
            if not cur: return
            with open(os.path.join(td, f"chunk_p{p}_{k:05d}.txt"), "w") as f:
                for c in cur: f.write(f"{fmt(c[0])} {fmt(c[1])} {fmt(c[2])} {c[3]}\n")
            nch += 1
        for c in ts:
            cur.append(c); cc += cost(c)
            if cc >= a.chunk_sec:
                flush(); cur, cc, k = [], 0.0, k + 1
        flush()
    print(f"wrote levels_sgx.txt ({len(levels)}), sgx_channels.txt ({len(A) + len(B)}), sgx_alphas.txt, {nch} chunks in tasks/")
