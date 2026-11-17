#!/usr/bin/env python3
"""Extension of the self-gravity system to every scattering permutation of the (n,8,8) levels (requested 2026-11-16).

The (n,8,8) levels, n = 9..45, are already in the system (Table I), with the channels of the paper's list and every
rule-allowed channel that also contains an (n,4,4) level (n <= 20).  This adds every other rule-allowed channel
(rules of rate_sve/gen_custom_input.py) among the levels with n <= 20 that contains an (n,8,8) level.

    python3 plan_n88.py            # report
    python3 plan_n88.py --write    # append the channels to sgx_channels.txt (backup kept) and write task chunks
                                   # tasks/chunk_n88a_*.txt (escaping fourth quantum) and chunk_n88b_*.txt (absorbed)
"""
import os, sys, shutil, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from states_lib import *
HERE = os.path.dirname(os.path.abspath(__file__))
ALPHAS = [float(x) for x in open(os.path.join(HERE, "sgx_alphas.txt")).read().split(",")]
st = lambda t: tuple(int(x) for x in t.split("-")) if "-" in t else tuple(int(c) for c in t)
key = lambda c: (frozenset([c[0], c[1]]), c[2])
LIST = os.path.join(HERE, "sgx_channels.txt")
have = [(st(p[0]), st(p[1]), st(p[2]), p[3]) for p in (l.split() for l in open(LIST)) if len(p) >= 4]
hk = set(key(c) for c in have)
S = sorted(set(expand(T1_10)) | set(SG21) | set(N44))
S20 = [s for s in S if s[0] <= 20]
new = [c for c in channels(S20) if any(s[1:] == (8, 8) for s in c[:3]) and key(c) not in hk]
cost = lambda c: sum((6.0 + (4.0 / x if c[3] == "BH" else 0.0)) for x in ALPHAS) + 5.0 * sum(1 for s in c[:3] if s[0] <= 8)
print(f"existing channels {len(have)}; new (n,8,8) channels {len(new)} (BH {sum(c[3] == 'BH' for c in new)}, Inf {sum(c[3] == 'Inf' for c in new)}); "
      f"est. {sum(cost(c) for c in new) / 3600:.0f} core-h")
if "--write" in sys.argv:
    bak = LIST + ".before_n88"
    if not os.path.exists(bak): shutil.copy(LIST, bak)
    with open(LIST, "a") as f:
        for c in new: f.write("    ".join((fmt(c[0]), fmt(c[1]), fmt(c[2]), c[3])) + " \n")
    td = os.path.join(HERE, "tasks"); n = 0
    for tag, sel in (("n88a", [c for c in new if c[3] == "Inf"]), ("n88b", [c for c in new if c[3] == "BH"])):
        cur, cc, k = [], 0.0, 0
        for c in sorted(sel, key=lambda c: (c[0], c[1], c[2])) + [None]:
            if c is not None:
                cur.append(c); cc += cost(c)
            if cur and (c is None or cc >= 5400.0):
                with open(os.path.join(td, f"chunk_{tag}_{k:05d}.txt"), "w") as f:
                    for x in cur: f.write(f"{fmt(x[0])} {fmt(x[1])} {fmt(x[2])} {x[3]}\n")
                cur, cc, k = [], 0.0, k + 1; n += 1
    print(f"appended {len(new)} channels to sgx_channels.txt (backup {os.path.basename(bak)}), wrote {n} chunks")
