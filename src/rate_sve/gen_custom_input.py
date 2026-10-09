#!/usr/bin/env python3
"""Write a load_rate_input file (S1 S2 S3 BH|Inf per line) for a chosen set of states.

Uses the same selection rules as Print_all_levels.py:
  * energies at alpha = 0.1, erg(n) = 1 - alpha^2/(2n^2) - alpha^4/(8n^4) (in units of mu);
  * erg(n1) + erg(n2) - erg(n3) > 1  ->  the fourth quantum escapes:      S1 S2 S3 Inf
  * otherwise                         ->  absorbed by the BH, allowed only if
                                          m1 + m2 = m3, l1 + l2 = l3 and l1 + l2 + l3 even:  S1 S2 S3 BH
  * (S1, S2, S3) and (S2, S1, S3) are the same channel; the ordering kept is the first one met when looping
    over states sorted by (n, l, m), i.e. the ordering used for the rate_sve/*_LvrHc_.dat file names.

usage:  python3 gen_custom_input.py 211 322 433 544 766 [-o load_rate_input_min_211_322_433_544_766.txt]
States are given as 'nlm' (single digits) or 'n-l-m'.
"""
import argparse
import os

ALPHA = 0.1


def parse_state(s):
    parts = s.split("-") if "-" in s else list(s)
    n, l, m = (int(p) for p in parts)
    if not (1 <= m <= l < n):
        raise ValueError(f"need 1 <= m <= l < n, got {s}")
    return n, l, m


def fmt(n, l, m):
    return f"{n}{l}{m}" if (n < 10 and l < 10 and m < 10) else f"{n}-{l}-{m}"


def erg(n, alph=ALPHA):
    return 1.0 - alph**2 / (2 * n**2) - alph**4 / (8 * n**4)


def channels(states):
    states = sorted(states)
    out, seen = [], set()
    for s1 in states:
        for s2 in states:
            for s3 in states:
                (n1, l1, m1), (n2, l2, m2), (n3, l3, m3) = s1, s2, s3
                key = (frozenset([s1, s2]) if s1 != s2 else frozenset([s1]), s3)
                if key in seen:
                    continue
                if erg(n1) + erg(n2) - erg(n3) > 1.0:
                    dest = "Inf"
                elif (m1 + m2 == m3) and (l1 + l2 == l3) and ((l1 + l2 + l3) % 2 == 0):
                    dest = "BH"
                else:
                    continue
                seen.add(key)
                out.append((fmt(*s1), fmt(*s2), fmt(*s3), dest))
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("states", nargs="+")
    ap.add_argument("-o", "--output", default=None)
    args = ap.parse_args()
    states = [parse_state(s) for s in args.states]
    chans = channels(states)
    fout = args.output or "load_rate_input_min_" + "_".join(fmt(*s) for s in sorted(states)) + ".txt"
    fout = os.path.join(os.path.dirname(os.path.abspath(__file__)), fout) if not os.path.isabs(fout) else fout
    with open(fout, "w") as f:
        for c in chans:
            f.write("    ".join(c) + " \n")
    print(f"{len(chans)} channels -> {fout}")
    for c in chans:
        print("  " + " ".join(c))
