"""State lists for the extended self-gravity system (2026-11-12).
T1_10 / T1_1E4: Table I of 'Toward the minimal set of states in axion superradiance' (draft DESY-26-156),
columns M = 10 and 1e4 Msun: (l, m) -> n values. SG21: the 21 levels of the first SG campaign.
"""
def rng(a, b): return list(range(a, b + 1))
T1_10 = {(1,1):[2,4], (4,1):[5], (2,2):[3,4,6,7], (4,2):[5], (3,3):rng(4,8), (4,3):[5], (6,3):[7], (7,3):[8],
         (4,4):rng(5,9), (5,4):[6,7], (6,4):[7], (7,4):[8], (5,5):rng(6,8), (7,5):[8],
         (6,6):rng(7,11)+rng(13,45), (7,6):[8], (8,6):[9], (12,6):[13], (7,7):[8,9], (8,8):rng(9,45),
         (10,8):[11], (12,8):[13], (10,10):rng(11,45), (12,12):[13], (14,14):[15], (16,16):rng(17,45),
         (18,18):[19], (20,20):rng(21,45)}
T1_1E4 = {(1,1):[2,4], (2,2):[3,4], (4,2):[5], (3,3):rng(4,8), (6,3):[7], (7,3):[8], (4,4):rng(5,9), (6,4):[7],
          (7,4):[8], (5,5):rng(6,8), (7,5):[8], (6,6):rng(7,11)+rng(14,45), (7,6):[8], (7,7):[8,9],
          (8,8):rng(9,45), (10,10):rng(11,45), (12,12):[13], (16,16):rng(17,45), (20,20):rng(22,45)}
def expand(t): return sorted((n, l, m) for (l, m), ns in t.items() for n in ns)
SG21 = [(2,1,1),(3,2,2),(4,1,1),(4,2,2),(4,3,3),(5,3,3),(5,4,4),(6,3,3),(6,4,4),(7,3,3),(7,4,4),(7,6,6),(8,3,3),
        (8,6,6),(9,6,6)] + [(n,6,6) for n in range(10,16)]
N44 = [(n,4,4) for n in range(5,21)]
# Fig. 2 (M = 10) and Fig. 23 (M = 1e4) of the draft, columns alpha = 0.05 and 0.2 (any f_a 1e12-1e18)
FIG_A02 = [(2,1,1),(3,2,2),(4,1,1),(4,2,2),(4,3,3),(5,4,1),(5,4,2),(5,4,4),(6,2,2),(6,4,4),(6,5,5),(7,2,2),(7,4,4),
           (7,5,5),(7,6,3),(7,6,4),(7,6,6),(8,7,3),(8,7,4),(8,7,7),(9,4,4),(9,6,6),(9,8,8),(10,6,6),(10,8,8),(11,6,6),
           (11,10,10),(13,12,6),(13,12,12),(14,6,6),(15,6,6)] + [(n,6,6) for n in range(16,46)] + [(n,8,8) for n in range(16,21)]
def fmt(s):
    n, l, m = s
    return f"{n}{l}{m}" if (n < 10 and l < 10 and m < 10) else f"{n}-{l}-{m}"
def erg(n, alph=0.1): return 1.0 - alph**2 / (2 * n**2) - alph**4 / (8 * n**4)
def channels(states):
    """gen_custom_input.py rules (same ordering convention)."""
    states = sorted(states); out, seen = [], set()
    for s1 in states:
        for s2 in states:
            for s3 in states:
                key = (frozenset([s1, s2]), s3)
                if key in seen: continue
                (n1,l1,m1),(n2,l2,m2),(n3,l3,m3) = s1, s2, s3
                if erg(n1) + erg(n2) - erg(n3) > 1.0: dest = "Inf"
                elif (m1+m2 == m3) and (l1+l2 == l3) and ((l1+l2+l3) % 2 == 0): dest = "BH"
                else: continue
                seen.add(key); out.append((s1, s2, s3, dest))
    return out
