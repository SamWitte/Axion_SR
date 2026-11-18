#!/usr/bin/env python3
"""Final spin, final mass and peak occupations of every run of a set (runs/<set>/gw_*/M_*/<variant>/alpha_*), cached
in runsum_<set>.pkl.   summary(set) -> {(gw, M, variant, alpha): dict(t, a, M, pk)}
    python3 runsum.py <set> [<set> ...]        (re)build the caches
    python3 runsum.py --table <ref> <set> [<set> ...]   one line per run: final spin of <ref> and the change in each set
"""
import glob, os, sys, pickle, warnings
import numpy as np
from multiprocessing import Pool
warnings.simplefilter("ignore", DeprecationWarning)
H = os.path.dirname(os.path.abspath(__file__)); X = os.path.dirname(H)
lab = lambda m: f"{m[0]}{m[1]}{m[2]}" if max(m) < 10 else f"{m[0]}-{m[1]}-{m[2]}"


def load(d):
    g = lambda k: (glob.glob(os.path.join(d, k + "_*.dat")) or [None])[0]
    if g("Time") is None or g("States") is None:
        return None
    t = np.loadtxt(g("Time")); spin = np.loadtxt(g("Spin")); mbh = np.loadtxt(g("MassBH"))
    modes = [tuple(int(float(x)) for x in l.split()[:3]) for l in open(g("Modes"))]
    pk = {}
    with open(g("States")) as f:
        for i, line in enumerate(f):
            p = np.fromstring(line, sep="\t").max()
            if p > 1e-20:
                pk[modes[i]] = p
    return dict(t=float(t[-1]), a=float(spin[-1]), M=float(mbh[-1] / mbh[0]), pk=pk)


def summary(name, refresh=False):
    cache = os.path.join(H, f"runsum_{name}.pkl")
    old = pickle.load(open(cache, "rb")) if os.path.exists(cache) and not refresh else {}
    dirs = sorted(glob.glob(f"{X}/runs/{name}/gw_*/M_*/*/alpha_*"))
    keyof = lambda d: (d.split("/")[-4][3:], d.split("/")[-3][2:], d.split("/")[-2], d.split("/")[-1][6:])
    todo = [d for d in dirs if keyof(d) not in old or old[keyof(d)] is None or old[keyof(d)]["t"] < 0.99e9]
    if todo:
        with Pool(min(12, len(todo))) as p:
            res = p.map(load, todo)
        for d, r in zip(todo, res):
            old[keyof(d)] = r
        pickle.dump(old, open(cache, "wb"))
    return old


if __name__ == "__main__":
    if "--table" in sys.argv:
        names = [a for a in sys.argv[1:] if not a.startswith("--")]
        S = {n: summary(n) for n in names}
        ref = S[names[0]]
        print("# final spin (reference set %s) and its change in the sets %s; peak changes > 20%% among the reference's top 8 levels" % (names[0], ", ".join(names[1:])))
        for k in sorted(ref, key=lambda k: (k[0], float(k[1]), k[2], float(k[3]))):
            r = ref[k]
            if r is None: continue
            out = f"{k[0]:6s} M={k[1]:<6s} {k[2]:14s} alpha {k[3]:<5s} a_f {r['a']:.4f}"
            for n in names[1:]:
                e = S[n].get(k)
                if e is None or e["t"] < 0.99e9:
                    out += f" | {n}: -"; continue
                top = sorted(r["pk"], key=lambda s: -r["pk"][s])[:8]
                ch = [f"{lab(s)} x{e['pk'].get(s, 0) / r['pk'][s]:.2g}" for s in top if abs(np.log(max(e["pk"].get(s, 0) / r["pk"][s], 1e-300))) > np.log(1.2)]
                out += f" | {n}: {e['a'] - r['a']:+.4f} M {e['M'] - r['M']:+.4f}" + (" [" + ", ".join(ch) + "]" if ch else "")
            print(out)
        for n in names[1:]:
            d = [(abs(S[n][k]["a"] - ref[k]["a"]), k) for k in ref if ref[k] and S[n].get(k) and S[n][k]["t"] >= 0.99e9]
            if d:
                d.sort(reverse=True)
                print(f"# {n}: {len(d)} runs, largest |change of final spin| {d[0][0]:.4f} ({' '.join(d[0][1])}); runs above 0.001: {sum(x[0] > 1e-3 for x in d)}; above 0.0001: {sum(x[0] > 1e-4 for x in d)}")
    else:
        for n in sys.argv[1:]:
            s = summary(n)
            print(n, len(s), "runs,", sum(1 for v in s.values() if v and v["t"] >= 0.99e9), "complete")
