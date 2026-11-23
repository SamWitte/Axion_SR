#!/usr/bin/env python3
"""Complete first post-Newtonian correction of the self-gravity amplitude of a channel N1 x N2 -> N3 x {BH, Inf},
term by term, relative to the amplitude of the rate tables A_NR (Newtonian exchange plus Newtonian recoil of the hole,
Klein-Gordon mode functions, Poisson equation in the Boyer-Lindquist radius).  Assembled from the pieces computed by
pieces_job.sh.  Derivation: Section 3.7 of ongoing_tests/sg_relativistic/note_source/sg_source_note.pdf.

Terms (amplitudes per multipole l4 of mode 4):
  gravitomagnetic   a_N * delta_GM                         (sg_1pn_gen.jl)
  hole potential    D(wrec) - a_N                          direct exchange in the harmonic radius x = r_BL - GM with the
                                                           weights 4 V_B (where the field acts) and 6 V_B (at the source)
  binding energy    D(wall) - D(wrec)                      frequency factors of the direct exchange
  retardation       a_N * delta_ret                        (sg_1pn_gen.jl)
  recoil            R(wall) - R(n0) + a_N * delta_ind      1PN correction of the recoil term: radial coordinate, weights
                                                           5/2 V_B and 6 V_B, retardation and gravitomagnetic induction in
                                                           the acceleration of the hole, frequency factors
  s-channel         12 pi G M_pl^2 a_SI / (8 alpha^2)      the pair term (what load_rates.jl adds)
  total             the sum
with D = a_N the direct (fixed-hole) amplitude of a set and R = a_NR - a_N its recoil term.

usage: complete_1pn.py <workdir> [--rel <dir>] [--tasks]
   prints, per channel and alpha, the terms divided by A_NR in units of alpha^2 (leading multipole) and the change of
   the rate, sum over multipoles of |A_NR + total|^2 / sum |A_NR|^2 - 1.
   --rel <dir>: result files of scan_gen.jl without newton_only (same names): adds (A_GR - A_NR)/A_NR of the
                relativistic t/u-channel exchange, to be compared with the total without the s-channel term.
   --tasks: print the task lines of sg_1pn_gen.jl ("pn ...") and sg_1pn_recoil.jl ("ind ...") from <workdir>/n0/res
"""
import glob, math, os, sys
import numpy as np

GMPL2 = 7484169213.942707 / 1.1219176708765e+66 * (1.22e19 * 1e9) ** 2          # G M_pl^2 with the constants of the code
NAMES = ["gm", "pot", "en", "ret", "rec", "pair"]
z = lambda n, i: complex(float(n[i]), float(n[i + 1]))
st = lambda s: tuple(int(x) for x in s.split("-")) if "-" in s else tuple(int(c) for c in s)
dipole = lambda S1, S2, S3: any(abs(st(s)[2] - st(S3)[2]) == 1 for s in (S1, S2))


def rows(path):
    return {int(r[10]): r for r in (l.split() for l in open(path) if l.strip())} if os.path.exists(path) else None


def table(path):
    return {int(r[0]): r for r in np.loadtxt(path, ndmin=2)} if os.path.exists(path) else None


def terms(W, f):
    """terms of one point: {l4: dict}, or None if a piece is missing"""
    b = os.path.basename(f); n0 = rows(f); q = next(iter(n0.values()))
    S1, S2, S3, als, dest = q[0], q[1], q[2], q[3], q[8]; a_ = float(als)
    wr, wa = rows(os.path.join(W, "wrec", "res", b)), rows(os.path.join(W, "wall", "res", b))
    pn = table(os.path.join(W, "pn", f"{S1}_{S2}_{S3}_{dest}_{als}.dat")); dip = dipole(S1, S2, S3)
    ind = table(os.path.join(W, "ind", f"{S1}_{S2}_{S3}_{dest}_{als}.dat")) if dip else {}
    if None in (wr, wa, pn, ind): return None
    per = {}
    for l4, q in n0.items():
        if l4 not in wr or l4 not in wa or l4 not in pn or (dip and l4 not in ind) or not all(math.isfinite(x) for x in pn[l4]): continue
        aSI, aN, aNR = z(q, 11), z(q, 13), z(q, 15)
        D = lambda x: z(x[l4], 13); R = lambda x: z(x[l4], 15) - z(x[l4], 13)
        t = dict(aNR=aNR, aN=aN, gm=aN * complex(pn[l4][3], pn[l4][4]), ret=aN * complex(pn[l4][9], pn[l4][10]), pot=D(wr) - aN, en=D(wa) - D(wr),
                 rec=(R(wa) - (aNR - aN)) + (aN * complex(ind[l4][3], ind[l4][4]) if dip else 0.0), pair=12 * math.pi * GMPL2 * aSI / (8 * a_ ** 2))
        t["total"] = sum(t[k] for k in NAMES); per[l4] = t
    return per or None


def main():
    W = sys.argv[1]; rel = sys.argv[sys.argv.index("--rel") + 1] if "--rel" in sys.argv else None
    key = lambda f: (os.path.basename(f).rsplit("_", 1)[0], float(os.path.basename(f)[:-4].rsplit("_", 1)[1]))
    base = sorted(glob.glob(os.path.join(W, "n0", "res", "*.txt")), key=key)
    if "--tasks" in sys.argv:
        for f in base:
            rr = [l.split() for l in open(f) if l.strip()]; q = rr[0]
            t = f"{q[0]} {q[1]} {q[2]} {q[8]} {q[3]} {','.join(r[10] for r in rr)}"
            print("pn", t)
            if dipole(q[0], q[1], q[2]): print("ind", t)
        return
    last = None
    for f in base:
        per = terms(W, f)
        if per is None: continue
        q = open(f).readline().split(); a_ = float(q[3]); t = per[min(per)]; r = lambda k: (t[k] / t["aNR"]).real / a_ ** 2
        rate = sum(abs(x["aNR"] + x["total"]) ** 2 for x in per.values()) / sum(abs(x["aNR"]) ** 2 for x in per.values()) - 1
        line = f"  {a_:6.3f} | {r('gm'):+9.4f} {r('pot'):+9.4f} {r('en'):+9.4f} {r('ret'):+9.4f} {r('rec'):+9.4f} {r('pair'):+9.4f} | {r('total'):+9.4f} | {100 * rate:+9.3f}"
        if rel:
            g = rows(os.path.join(rel, os.path.basename(f)))
            if g is not None and min(per) in g and math.isfinite(float(g[min(per)][17])):
                line += f" | {r('total') - r('pair'):+9.4f} {((z(g[min(per)], 17) - t['aNR']) / t['aNR']).real / a_ ** 2:+9.4f}"
        if key(f)[0] != last:
            last = key(f)[0]
            print(f"\n{q[0]} x {q[1]} -> {q[2]} x {q[8]}   leading multipole l4 = {min(per)},  A_NR/A_N = {(t['aNR'] / t['aN']).real:+.4f}")
            print("   alpha | term / A_NR / alpha^2: gravitomagnetic, hole potential, binding energy, retardation, recoil, s-channel | total | rate change [%]"
                  + (" | total without s-channel, relativistic exchange" if rel else ""))
        print(line)


if __name__ == "__main__":
    main()
