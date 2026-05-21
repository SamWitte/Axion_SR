#!/usr/bin/env python3
"""
Generate a filtered version of a load_rate_input file keeping only lines
where every state has l == m, with the exception that "765" is always kept
(even though it has l=6, m=5).

Usage:
  python gen_lm_input.py                          # default: Nmax_15
  python gen_lm_input.py load_rate_input_Nmax_8.txt
"""

import os
import sys

HERE      = os.path.dirname(os.path.abspath(__file__))
EXCEPTION = {"765"}   # kept even though l != m


def parse_lm(token):
    """Return (l, m) from 'nlm' (3-char) or 'n-l-m' (hyphenated) format."""
    if "-" in token:
        parts = token.split("-")
        return int(parts[1]), int(parts[2])
    else:
        return int(token[1]), int(token[2])


def keep_line(s1, s2, s3):
    """Keep if every state has l==m (treating EXCEPTION tokens as l==m)."""
    for tok in (s1, s2, s3):
        if tok in EXCEPTION:
            continue
        l, m = parse_lm(tok)
        if l != m:
            return False
    return True


def filter_file(src):
    base, ext = os.path.splitext(src)
    dst = base + "_lm" + ext

    with open(src) as fh:
        lines = [l for l in fh if l.strip()]

    kept = []
    for line in lines:
        parts = line.split()
        if len(parts) < 3:
            continue
        if keep_line(parts[0], parts[1], parts[2]):
            kept.append(line)

    with open(dst, "w") as fh:
        fh.writelines(kept)

    print(f"Input  : {src}  ({len(lines)} lines)")
    print(f"Output : {dst}  ({len(kept)} lines kept, {len(lines)-len(kept)} removed)")
    return dst


if __name__ == "__main__":
    if len(sys.argv) > 1:
        src = sys.argv[1]
        if not os.path.isabs(src):
            src = os.path.join(HERE, src)
    else:
        src = os.path.join(HERE, "load_rate_input_Nmax_15.txt")

    filter_file(src)
