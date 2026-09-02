"""Extend load_rate_input_Nmax_15.txt with:

  1. the full l=m shell at n=9, i.e. (9,l,l) for l=1..8;
  2. the l=m=8 ladder extended to n=12..15, i.e. (12,8,8)..(15,8,8);
  3. all permutations in which at least one leg is one of those states,
     with the partner states drawn from the *complete* n<=8 block (every
     l, every m) rather than the (l>=3, m>=3) corner used previously;
  4. the truncation ("cap") processes that stop those ladders growing,
     (n,l,l) x (n,l,l) -> (2l+1, 2l, 2l) + BH, plus the inverse.

Selection rules are identical to Print_all_levels.py: emission to infinity
whenever the energy budget allows it, otherwise a BH-absorption channel if
l1+l2 = l3, m1+m2 = m3 and l1+l2+l3 is even.

Output: rate_sve/load_rate_input_Nmax_15_ext.txt
"""

import os

Nmax = 15
BASE_FILE = "rate_sve/load_rate_input_Nmax_15.txt"
OUT_FILE = "rate_sve/load_rate_input_Nmax_15_ext.txt"
# the added rows only (everything in OUT_FILE that is not already in BASE_FILE)
DIFF_FILE = "rate_sve/load_rate_input_Nmax_15_ext_only.txt"

# The l=m=10 ladder is tracked to n=15 in the Nmax_15 file but never capped
# (the |21,20,20> truncation only appears in the Nmax_18 file). Set False to
# reproduce that omission.
INCLUDE_L10_CAP = True

HERE = os.path.dirname(os.path.abspath(__file__))


def erg_shift_1(n, alph=0.1):
    return 1.0 * (1.0 - alph**2 / (2 * n**2) - alph**4 / (8 * n**4))


def format_state(n, l, m):
    if n < 10 and l < 10 and m < 10:
        return "{}{}{}".format(n, l, m)
    else:
        return "{}-{}-{}".format(n, l, m)


def parse_state(tok):
    if "-" in tok:
        return tuple(int(x) for x in tok.split("-"))
    return tuple(int(c) for c in tok)


# ---------------------------------------------------------------- state pools

# complete n<=8 block
block = [(n, l, m) for n in range(2, 9) for l in range(1, n) for m in range(1, l + 1)]

# high-n l=m ladders already present in the Nmax_15 file
ladders_old = (
    [(n, 6, 6) for n in range(9, 16)]
    + [(n, 8, 8) for n in range(9, 12)]
    + [(n, 10, 10) for n in range(11, 16)]
)

# newly tracked states
new_states = (
    [(9, l, l) for l in range(1, 9)]
    + [(n, 8, 8) for n in range(12, 16)]
)

pool = sorted(set(block) | set(ladders_old) | set(new_states))
NEW = set(new_states)


# ------------------------------------------------------------------ selection

def classify(s1, s2, s3):
    """Return 'Inf', 'BH' or None for the process s1 x s2 -> s3 x (Inf|BH)."""
    n1, l1, m1 = s1
    n2, l2, m2 = s2
    n3, l3, m3 = s3

    erg_diff = erg_shift_1(n1) + erg_shift_1(n2) - erg_shift_1(n3)

    if erg_diff > 1.0:
        return "Inf"

    if (m1 + m2 - m3) == 0 and (l1 + l2 - l3) == 0 and (l1 + l2 + l3) % 2 == 0:
        return "BH"

    return None


rows = []          # ordered list of (s1, s2, s3, tag)
seen = set()       # canonical keys, initial pair unordered


def try_add(s1, s2, s3, tag):
    key = (min(s1, s2), max(s1, s2), s3, tag)
    if key in seen:
        return False
    seen.add(key)
    rows.append((s1, s2, s3, tag))
    return True


# ------------------------------------------------ 1. carry over the base file

n_base = 0
with open(os.path.join(HERE, BASE_FILE)) as fh:
    for line in fh:
        p = line.split()
        if len(p) < 4:
            continue
        if try_add(parse_state(p[0]), parse_state(p[1]), parse_state(p[2]), p[3]):
            n_base += 1

# ---------------------------------------------- 2. permutations with new legs

n_perm = 0
for s1 in pool:
    for s2 in pool:
        if s1 > s2:
            continue
        for s3 in pool:
            if not (s1 in NEW or s2 in NEW or s3 in NEW):
                continue
            tag = classify(s1, s2, s3)
            if tag and try_add(s1, s2, s3, tag):
                n_perm += 1

# ------------------------------------------------------- 3. truncation / caps

# every l=m ladder that gained members needs its growth capped; the cap state
# |2l+1, 2l, 2l> is only reachable through these dedicated processes.
cap_sources = sorted(NEW | {(n, 8, 8) for n in range(9, 12)})
if INCLUDE_L10_CAP:
    cap_sources += [(n, 10, 10) for n in range(11, 16)]

n_cap = 0
for (n, l, m) in cap_sources:
    cap = (2 * l + 1, 2 * l, 2 * m)

    tag = classify((n, l, m), (n, l, m), cap)
    if tag and try_add((n, l, m), (n, l, m), cap, tag):
        n_cap += 1

    # inverse: two quanta in the cap state scatter back down
    tag = classify(cap, cap, (n, l, m))
    if tag == "Inf" and try_add(cap, cap, (n, l, m), tag):
        n_cap += 1

# ------------------------------------------------------------------- 4. write

def write_rows(path, rws):
    with open(os.path.join(HERE, path), "w") as fh:
        for (s1, s2, s3, tag) in rws:
            fh.write(
                format_state(*s1) + "    " + format_state(*s2) + "    "
                + format_state(*s3) + "    " + tag + " \n"
            )


write_rows(OUT_FILE, rows)
# rows[:n_base] are the carried-over base rows; everything after is new
write_rows(DIFF_FILE, rows[n_base:])

states_used = set()
for (s1, s2, s3, tag) in rows:
    states_used |= {s1, s2, s3}

print("base rows carried over : {}".format(n_base))
print("new permutation rows   : {}".format(n_perm))
print("new truncation rows    : {}".format(n_cap))
print("total rows written     : {}  -> {}".format(len(rows), OUT_FILE))
print("added rows only        : {}  -> {}".format(len(rows) - n_base, DIFF_FILE))
print("distinct states        : {}  (max n = {})".format(
    len(states_used), max(s[0] for s in states_used)))
