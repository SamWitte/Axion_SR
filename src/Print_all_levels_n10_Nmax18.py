import numpy as np


def erg_shift_1(n, alph=0.1):
    return 1.0 * (1.0 - alph**2 / (2 * n**2) - alph**4 / (8 * n**4))


def format_state(n, l, m):
    if n < 10 and l < 10 and m < 10:
        return "{}{}{}".format(n, l, m)
    else:
        return "{}-{}-{}".format(n, l, m)


# New states — at least one must appear in each recorded scattering
new_states = {(16, 10, 10), (17, 10, 10), (18, 10, 10)}

all_states = [
    (5, 4, 4),
    (6, 4, 4),
    (6, 5, 5),
    (7, 4, 4),
    (7, 5, 5),
    (7, 6, 6),
    (8, 5, 5),
    (8, 6, 6),
    (8, 7, 7),
    (9, 8, 8),
    (11, 8, 8),
    (11, 10, 10),
    (12, 10, 10),
    (13, 10, 10),
    (13, 12, 12),
    (14, 10, 10),
    (15, 10, 10),
    (16, 10, 10),
    (17, 10, 10),
    (18, 10, 10),
]

file_output = "rate_sve/load_rate_input_Nmax_18.dat"

hold_out_1 = []
hold_out_2 = []
hold_out_3 = []
hold_out_4 = []
final_keep_sve = []


def try_add(tag1, tag2, tag3, tag4, outtag):
    for i in range(len(hold_out_1)):
        if (
            (hold_out_1[i] == tag1 and hold_out_2[i] == tag2 and hold_out_3[i] == tag3)
            or (hold_out_1[i] == tag2 and hold_out_2[i] == tag1 and hold_out_3[i] == tag3)
        ):
            return
    hold_out_1.append(tag1)
    hold_out_2.append(tag2)
    hold_out_3.append(tag3)
    hold_out_4.append(tag4)
    final_keep_sve.append(tag1[1:-1] + "    " + tag2[1:-1] + "    " + tag3[1:-1] + "    " + outtag)


def process_triple(n1, l1, m1, n2, l2, m2, n3, l3, m3):
    tag1 = "|{}>".format(format_state(n1, l1, m1))
    tag2 = "|{}>".format(format_state(n2, l2, m2))
    tag3 = "|{}>".format(format_state(n3, l3, m3))

    erg_diff = erg_shift_1(n1) + erg_shift_1(n2) - erg_shift_1(n3)

    if erg_diff > 1.0:
        lf = l1 + l2 - l3
        mf = m1 + m2 - m3
        tag4 = "Inf ({}-{})".format(lf, mf)
        try_add(tag1, tag2, tag3, tag4, "Inf")

    elif erg_diff > 0.999:
        test_m = ((m1 + m2 - m3) == 0)
        test_l = ((l1 + l2 - l3) == 0)
        test1 = ((l1 + l2 + l3) % 2 == 0)
        if test1 and test_m and test_l:
            try_add(tag1, tag2, tag3, "BH", "BH")

    else:
        test_m = ((m1 + m2 - m3) == 0)
        test_l = ((l1 + l2 - l3) == 0)
        test1 = ((l1 + l2 + l3) % 2 == 0)
        if test1 and test_m and test_l:
            try_add(tag1, tag2, tag3, "BH", "BH")


# Main loop: all triples from the pool where at least one is a new state
for (n1, l1, m1) in all_states:
    for (n2, l2, m2) in all_states:
        for (n3, l3, m3) in all_states:
            if not ((n1, l1, m1) in new_states or (n2, l2, m2) in new_states or (n3, l3, m3) in new_states):
                continue
            process_triple(n1, l1, m1, n2, l2, m2, n3, l3, m3)


# Truncation: two |n,10,10> states (at least one new) scattering to |21,20,20>
# (outside the normal pool, analogous to n66 -> |13,12,12> in Print_all_levels_n6.py)
n1010_states = [(n, 10, 10) for n in range(11, 19)]

trunc_extra = []
for (n1, l1, m1) in n1010_states:
    for (n2, l2, m2) in n1010_states:
        if not ((n1, l1, m1) in new_states or (n2, l2, m2) in new_states):
            continue
        m3 = m1 + m2  # = 20
        l3 = m3       # = 20
        n3 = l3 + 1   # = 21
        trunc_extra.append((n1, l1, m1, n2, l2, m2, n3, l3, m3))
        process_triple(n1, l1, m1, n2, l2, m2, n3, l3, m3)

# Inverse truncation: |21,20,20> x |21,20,20> -> |n,10,10>
for (n1, l1, m1, n2, l2, m2, n3, l3, m3) in trunc_extra:
    na, la, ma = n3, l3, m3  # the truncation state |21,20,20>
    nb, lb, mb = n3, l3, m3
    nc, lc, mc = n1, l1, m1  # back to original n1010 state

    erg_diff = erg_shift_1(na) + erg_shift_1(nb) - erg_shift_1(nc)
    tag1 = "|{}>".format(format_state(na, la, ma))
    tag2 = "|{}>".format(format_state(nb, lb, mb))
    tag3 = "|{}>".format(format_state(nc, lc, mc))

    if erg_diff > 1.0:
        lf = la + lb - lc
        mf = ma + mb - mc
        tag4 = "Inf ({}-{})".format(lf, mf)
        try_add(tag1, tag2, tag3, tag4, "Inf")


with open(file_output, 'w') as f:
    for line in final_keep_sve:
        f.write(line + " \n")

print("Wrote {} scatterings to {}".format(len(final_keep_sve), file_output))
