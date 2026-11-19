#!/usr/bin/env python3
"""Complete first post-Newtonian correction of the self-gravity amplitude of a channel N1 x N1 -> N3 x {BH, Inf},
relative to the Newtonian amplitude of the rate tables (Klein-Gordon mode functions, Poisson equation solved in the
Boyer-Lindquist radius), assembled from the pieces computed by pieces_job.sh.  Derivation: Section 3.7 of
ongoing_tests/sg_relativistic/note_source/sg_source_note.pdf.

Printed per channel and alpha, as coefficients of alpha^2 of (A - A_NR)/A_N for the leading multipole of mode 4:
  direct exchange   radial coordinate (harmonic x = r_BL - GM), potential of the hole where the field acts (4 V_B) and at
                    the source (6 V_B), frequency factors E = 2 (w_a w_4 + w_b w_3 - 2 mu^2)/mu^2, gravitomagnetism,
                    retardation
  recoil            (only for |m1 - m3| = 1) radial coordinate, 5/2 V_B where the field acts, 6 V_B at the source,
                    retardation and gravitomagnetic induction in the acceleration of the hole, frequency factors
  pair term         the s-channel graviton exchange, 12 pi G M_pl^2 a_SI / (8 alpha^2 a_N)  (what load_rates.jl adds)
and their sum.  With --rel <dir> (result files of scan_gen.jl without newton_only, same names) the relativistic value
(A_GR - A_NR)/A_N of the t/u-channel exchange is printed next to the sum without the pair term.

Frequency factors are taken at leading order (hydrogenic energies): -4/n1^2 for the direct term and
-(5/(2 n1^2) + 1/(2 n3^2)) times the Newtonian recoil term for the recoil.  Incoming levels must be identical: for
N1 != N2 the two exchanges carry different frequency factors and are not separated in the stored amplitudes.

usage: complete_1pn.py <workdir> [--rel <dir>] [--tasks]
   <workdir>/{n0,shift,w40,w46,wrec}/res   Newtonian amplitudes (scan_gen.jl newton_only) without and with the switches
   <workdir>/pn                            1PN pieces of sg_1pn_gen.jl (gravitomagnetism and retardation are used)
   <workdir>/ind                           recoil pieces of sg_1pn_recoil.jl (induction)
   --tasks: print the task lines for sg_1pn_gen.jl and sg_1pn_recoil.jl from <workdir>/n0/res (used by pieces_job.sh)
"""
import glob, math, os, sys
import numpy as np

GMPL2 = 7484169213.942707 / 1.1219176708765e+66 * (1.22e19 * 1e9) ** 2          # G M_pl^2 with the constants of the code
z = lambda n, i: complex(float(n[i]), float(n[i + 1]))
st = lambda s: tuple(int(x) for x in s.split("-")) if "-" in s else tuple(int(c) for c in s)


def first(path):
    return open(path).readline().split() if os.path.exists(path) else None


def main():
    W = sys.argv[1]; rel = sys.argv[sys.argv.index("--rel") + 1] if "--rel" in sys.argv else None
    base = sorted(glob.glob(os.path.join(W, "n0", "res", "*.txt")), key=lambda f: (os.path.basename(f).rsplit("_", 1)[0], float(os.path.basename(f)[:-4].rsplit("_", 1)[1])))
    if "--tasks" in sys.argv:
        for f in base:
            rows = [l.split() for l in open(f) if l.strip()]; q = rows[0]
            n1, n3 = st(q[0]), st(q[2])
            print("pn", q[0], q[1], q[2], q[8], q[3], ",".join(r[10] for r in rows))
            if q[0] == q[1] and abs(n1[2] - n3[2]) == 1: print("ind", q[0], q[2], q[8], q[3])
        return
    last = None
    for f in base:
        q = first(f); S1, S2, S3, als, dest, l4 = q[0], q[1], q[2], q[3], q[8], int(q[10]); a_ = float(als)
        if S1 != S2: continue
        n1, n3 = st(S1), st(S3); b = os.path.basename(f)
        v = {k: first(os.path.join(W, k, "res", b)) for k in ("shift", "w40", "w46", "wrec")}
        pn = os.path.join(W, "pn", f"{S1}_{S2}_{S3}_{dest}_{als}.dat")
        if any(x is None for x in v.values()) or not os.path.exists(pn): continue
        p = {int(r[0]): r for r in np.loadtxt(pn, ndmin=2)}[l4]
        gm, ret = p[3] / a_ ** 2, p[9] / a_ ** 2
        a0 = z(q, 13); N = {k: z(x, 13) for k, x in v.items()}
        coord = ((N["shift"] / a0).real - 1) / a_ ** 2
        wx = ((N["w40"] - N["shift"]) / a0).real / a_ ** 2; wy = ((N["w46"] - N["w40"]) / a0).real / a_ ** 2
        E = -4.0 / n1[0] ** 2
        direct = coord + wx + wy + E + gm + ret
        pair = (12 * math.pi * GMPL2 * z(q, 11) / (8 * a_ ** 2 * a0)).real / a_ ** 2
        line = f"  {a_:6.3f} | {coord:+.4f} {wx:+.4f} {wy:+.4f} {E:+.4f} {gm:+.4f} {ret:+.4f} = {direct:+.4f} |"
        rec = 0.0; dip = abs(n1[2] - n3[2]) == 1
        if dip:
            ind = os.path.join(W, "ind", f"{S1}_{S1}_{S3}_{dest}_{als}.dat")
            if not os.path.exists(ind): continue
            R = lambda x: (z(x, 15) - z(x, 13)) / a0                              # recoil term / A_N
            R0 = R(q).real; rc = (R(v["shift"]).real - R0) / a_ ** 2
            rx = 2.5 * (R(v["w40"]) - R(v["shift"])).real / 4 / a_ ** 2; ry = (R(v["w46"]) - R(v["w40"])).real / a_ ** 2
            rr = (R(v["wrec"]) - R(v["shift"])).real / a_ ** 2 - rx - ry        # retardation in the acceleration of the hole
            Er = -(5.0 / (2 * n1[0] ** 2) + 1.0 / (2 * n3[0] ** 2)) * R0; di = np.loadtxt(ind)[3] / a_ ** 2
            rec = rc + rx + ry + rr + Er + di
            line += f" {rc:+.4f} {rx:+.4f} {ry:+.4f} {rr:+.4f} {Er:+.4f} {di:+.4f} = {rec:+.4f} |"
        line += f" {pair:+.4f} | {direct + rec + pair:+.4f}"
        if rel:
            g = first(os.path.join(rel, b))
            if g is not None and math.isfinite(float(g[17])):
                line += f" | {direct + rec:+.4f}  {((z(g, 17) - z(g, 15)) / z(g, 13)).real / a_ ** 2:+.4f}"
        if (S1, S3) != last:
            last = (S1, S3)
            print(f"\n{S1} x {S1} -> {S3} x {dest}   leading multipole l4 = {l4}" + (f", Newtonian recoil term / A_N = {R0:+.4f}" if dip else ", no recoil"))
            print("   alpha | direct: coordinate, 4 V_B(x), 6 V_B(y), frequencies, gravitomagnetism, retardation = sum |"
                  + (" recoil: coordinate, 5/2 V_B(x), 6 V_B(y), retardation, frequencies, induction = sum |" if dip else "")
                  + " pair term | complete 1PN" + (" | without pair term, relativistic exchange" if rel else ""))
        print(line)


if __name__ == "__main__":
    main()
