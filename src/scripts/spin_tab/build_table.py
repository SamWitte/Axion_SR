#!/usr/bin/env python3
"""Spin table of the scattering rates from the campaign output.   python3 build_table.py [--fast DIR] [--pts DIR] [--out FILE]
Reads the rows of spin_fast.jl (fast_out/fast_*.txt; fix998_*.txt for the a = 0.998 node of rows with a broken cached radial
function, see wf_bad.txt; fixlev_*.txt for the nodes a = 0, 0.3 of rows with a level energy that sat on a wrong root, see
lev_bad.txt; fast7_*.txt, the run with the n = 8 levels hydrogenic, only as a cross-check) and writes the table
the loader reads (NL_SPIN_TABLE):
    # spins 0 0.3 0.5 0.7 0.8 0.9 0.95 0.99 0.998
    S1 S2 S3 S4 alpha flag  lnR_SI[9]  lnR_SG[9]
    R = [rate(alpha, a) / rate(alpha, 0.95)] / rho,   rho = r+(a) / r+(0.95) for BH-labelled channels, 1 otherwise
    (the tabulated rate: without kH for horizon channels with m4 != 0).  R = 1 is today's behaviour.
flag: 2 horizon channel, m4 = 0;  3 the same with the a = 0.998 node recomputed (the cached relativistic radial function of
      one of its levels is broken at 0.998; the function of a = 0.99 is used there, wf_check.jl / wf_bad.txt);
      4 horizon channel, m4 != 0 (5: with the recomputed 0.998 node).   Channels not in the file: R = 1 (all channels to
      infinity, and whatever is listed in the report).
The report (<out>.report.txt) lists every special case, and compares with the full-path validation points in pts/."""
import glob, math, collections, argparse
X = "/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext"
ap = argparse.ArgumentParser()
ap.add_argument("--fast", default=X + "/spin_tab/fast_out"); ap.add_argument("--pts", default=X + "/spin_tab/pts")
ap.add_argument("--out", default=X + "/spin_tab/spin_table_v1.txt"); ap.add_argument("--channels", default=X + "/sgx_channels.txt")
A = ap.parse_args()
SP = [0.0, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99, 0.998]; HI = [a for a in SP if a > 0.95]
AL = [float(x) for x in open(X + "/sgx_alphas.txt").read().strip().split(",")]
rp = lambda a: 1 + math.sqrt(1 - a * a)
st = lambda t: tuple(int(x) for x in t.split("-")) if "-" in t else tuple(int(c) for c in t)
chans = [tuple(l.split()[:4]) for l in open(A.channels) if l.strip()]

def read_fast(prefix):
    d = collections.defaultdict(dict)                    # (channel, alpha) -> spin -> (status, si, sg)
    for f in glob.glob(A.fast + "/" + prefix + "_*.txt"):
        spins = None
        for ln in open(f):
            p = ln.split()
            if ln.startswith("# spins"):
                spins = [float(x) for x in p[2:p.index("hmax")]]; continue
            if spins is None or len(p) != 6 + 4 * len(spins): continue
            for j, a in enumerate(spins):
                q = p[6 + 4 * j: 10 + 4 * j]; d[(tuple(p[:4]), float(p[4]))][a] = (int(q[0]), float(q[1]), float(q[2]))
    return d
fast = read_fast("fast"); fast7 = read_fast("fast7"); fix = read_fast("fix998")
fixlev = read_fast("fixlev"); levbad = set(); nlev = 0   # level energies on a wrong root (lev_bad.txt): nodes 0, 0.3 recomputed
for l in open(X + "/spin_tab/lev_bad.txt"):
    if not l.startswith("#"):
        q = l.split(); levbad.add(((int(q[0]), int(q[1]), int(q[2])), float(q[4])))
for k, d in fixlev.items():
    if k in fast: fast[k].update(d); nlev += 1
bad = set()                                              # broken cached radial functions: (level, alpha), all at a = 0.998
for l in open(X + "/spin_tab/wf_bad.txt"):
    if not l.startswith("#"):
        q = l.split(); assert float(q[3]) == 0.998; bad.add(((int(q[0]), int(q[1]), int(q[2])), float(q[4])))
full = collections.defaultdict(dict)                     # validation points: (channel, alpha) -> spin -> (si, sg, bound4)
for f in glob.glob(A.pts + "/pts_*.txt"):
    for ln in open(f):
        p = ln.split()
        if len(p) == 20 and float(p[6]) == 1: full[(tuple(p[:4]), float(p[5]))][float(p[4])] = (float(p[8]), float(p[9]), int(p[14]))

def lnS(d, spins):
    """ln[rate(a)/rate(0.95)] (SI, SG) at the given spins; None unless every node is computed (status 1) and positive"""
    if any(a not in d or d[a][0] != 1 or min(d[a][1:]) <= 0 for a in list(spins) + [0.95]): return None
    return [[math.log(d[a][i] / d[0.95][i]) for a in spins] for i in (1, 2)]

rep = collections.defaultdict(list); nflag = collections.Counter(); rows = []; val = collections.defaultdict(list)
xchk = []; smooth = []; big = []; I99 = SP.index(0.99); I998 = SP.index(0.998); I95 = SP.index(0.95)
for c in chans:
    s = [st(x) for x in c[:3]]; m4 = s[0][2] + s[1][2] - s[2][2]; ns = [x[0] for x in s]
    lrho = [math.log(rp(a) / rp(0.95)) if c[3] == "BH" else 0.0 for a in SP]
    kind = "m4!=0" if (c[3] == "BH" and m4 != 0) else ("m4=0" if c[3] == "BH" else "Inf")
    for al in AL:
        k = (c, al); F = fast.get(k); flag = 0
        if F and k not in fixlev and any((x, al) in levbad for x in s):
            rep["level energy on a wrong root at spin 0 or 0.3 and no recomputed row: no table row (R = 1)"].append((c, al)); F = None
        ln = lnS(F, SP) if F else None
        if F and ln is None:
            stt = sorted({v[0] for v in F.values()})
            if stt != [2]:
                rep["row unusable: fourth quantum bound at some spins only (status 2), eigenvalue search failed (0) or rate <= 0"].append((c, al, "status", stt))
        if F is None and c[3] == "BH" and abs(m4) < 20: rep["horizon channel without a row (R = 1)"].append((c, al))
        if ln is not None:
            flag = 4 if kind == "m4!=0" else 2
            badl = [x for x in s if (x, al) in bad]
            if badl:
                x = fix.get(k, {}).get(0.998)
                if x and x[0] == 1 and min(x[1:]) > 0:
                    old = max(abs(ln[i][I998] - math.log(x[1 + i] / F[0.95][1 + i])) for i in (0, 1))
                    if old > 0.1: nflag["0.998 node changed by > 10 % by the recomputation"] += 1
                    for i in (0, 1): ln[i][I998] = math.log(x[1 + i] / F[0.95][1 + i])
                    h = lnS(fast7.get(k, {}), HI)
                    if h is not None and all(y[0] == 8 for y in badl) and not any((y, al) in bad for y in s if y[0] != 8):
                        xchk.append((max(abs((ln[i][I998] - ln[i][I99]) - (h[i][1] - h[i][0])) for i in (0, 1)), c, al))
                else:
                    rep["broken radial function at 0.998 but no recomputed value: 0.998 node set equal to the 0.99 node"].append((c, al))
                    for i in (0, 1): ln[i][I998] = ln[i][I99]
                flag += 1
            e = max(abs((ln[i][I998] - ln[i][I99]) - 0.45 * (ln[i][I99] - ln[i][I95])) for i in (0, 1))
            smooth.append((e, c, al))
            if c[3] != "BH": rep["Inf-labelled channel with a bound fourth quantum: gets a table row"].append((c, al))
            v = [[ln[i][j] - lrho[j] for j in range(len(SP))] for i in (0, 1)]
            big.append((max(abs(x) for x in v[0] + v[1]), c, al))
            rows.append(" ".join(c) + " %g %d " % (al, flag) + " ".join("%.5f" % x for x in v[0] + v[1]))
        nflag[(kind, flag)] += 1
        # validation against the full path (spin ratios, and absolute values at each spin)
        P = full.get(k)
        if P and F and 0.95 in P and min(P[0.95][:2]) > 0:
            for a, p in P.items():
                if a in F and F[a][0] == 1 and p[2] == 1 and min(p[:2]) > 0 and min(F[a][1:]) > 0:
                    ab = [abs(math.log(F[a][1 + i] / p[i])) for i in (0, 1)]
                    ok95 = F.get(0.95, (0,))[0] == 1 and a != 0.95
                    ra = [abs(math.log((F[a][1 + i] / F[0.95][1 + i]) / (p[i] / P[0.95][i]))) for i in (0, 1)] if ok95 else None
                    cls = kind if kind != "m4=0" else ("m4=0, all n<=8" if max(ns) <= 8 else ("m4=0, n>8 with a level n<=8" if min(ns) <= 8 else "m4=0, all n>8"))
                    val[(cls, al)].append((ab, ra, c, a, p[0] / p[1]))
with open(A.out, "w") as o:
    o.write("# spin table of the scattering rates, built by spin_tab/build_table.py; R = [rate(a)/rate(0.95)]/rho, natural logs; channels not listed: R = 1\n")
    o.write("# spins " + " ".join("%g" % a for a in SP) + "\n# S1 S2 S3 S4 alpha flag lnR_SI[%d] lnR_SG[%d]\n" % (len(SP), len(SP)))
    o.write("\n".join(rows) + "\n")
with open(A.out + ".report.txt", "w") as o:
    o.write("rows written: %d\n(class, flag) -> number of (channel, alpha):\n" % len(rows))
    for k, v in sorted((k, v) for k, v in nflag.items() if isinstance(k, tuple)): o.write("  %-6s flag %d : %d\n" % (k + (v,)))
    for k, v in nflag.items():
        if not isinstance(k, tuple): o.write("%s: %d\n" % (k, v))
    o.write("rows with the nodes a = 0, 0.3 recomputed after the correction of level energies (lev_bad.txt): %d\n" % nlev)
    q = lambda z, f: z[int(f * (len(z) - 1))] if z else float("nan")
    o.write("\nfast vs full path (validation points): |ln fast/full| of the rate, and of the spin ratio rate(a)/rate(0.95); SG | SI\n")
    o.write("  (SI outliers are channels in which the SI amplitude nearly cancels: SI/SG of the worst point is given)\n")
    for (cls, al), v in sorted(val.items()):
        o.write("  %-28s alpha %-5g n %4d" % (cls, al, len(v)))
        for i, nm in ((1, "SG"), (0, "SI")):
            ab = sorted(x[0][i] for x in v); ra = sorted(x[1][i] for x in v if x[1] is not None); w = max(v, key=lambda x: x[0][i])
            o.write("  | %s rate: med %.4f 90%% %.4f max %.4f (%s a=%g, SI/SG %.0e)  ratio: med %.4f 90%% %.4f max %.4f" %
                    (nm, q(ab, 0.5), q(ab, 0.9), ab[-1], " ".join(w[2][:3]), w[3], w[4], q(ra, 0.5), q(ra, 0.9), q(ra, 1.0)))
        o.write("\n")
    z = sorted(x[0] for x in xchk)
    o.write("\nrecomputed 0.998 node against the run with the n = 8 levels hydrogenic (rows whose only broken function is an n = 8 one): "
            "%d rows, |difference of ln rate(0.998)/rate(0.99)| median %.4f 90%% %.4f 99%% %.4f max %.4f\n" % (len(z), q(z, .5), q(z, .9), q(z, .99), q(z, 1.0)))
    for x in sorted(xchk, reverse=True)[:8]: o.write("  %.4f  %s  alpha %g\n" % (x[0], " ".join(x[1]), x[2]))
    z = sorted(x[0] for x in smooth)
    o.write("\nsmoothness at the top nodes, e = [lnR(0.998) - lnR(0.99)] - 0.45 [lnR(0.99) - lnR(0.95)] (0.45: ratio of the intervals in "
            "sqrt(1 - a^2)): %d rows, |e| median %.4f 90%% %.4f 99%% %.4f max %.3f, %d rows > 0.1\n" % (len(z), q(z, .5), q(z, .9), q(z, .99), q(z, 1.0), sum(x > 0.1 for x in z)))
    for x in sorted(smooth, reverse=True)[:25]: o.write("  %.3f  %s  alpha %g\n" % (x[0], " ".join(x[1]), x[2]))
    z = sorted(x[0] for x in big)
    o.write("\nsize of the correction, max |ln R| over the nodes of a row: median %.3f 90%% %.3f 99%% %.3f max %.2f; rows > 1: %d\n" %
            (q(z, .5), q(z, .9), q(z, .99), q(z, 1.0), sum(x > 1 for x in z)))
    for x in sorted(big, reverse=True)[:25]: o.write("  %.2f  %s  alpha %g\n" % (x[0], " ".join(x[1]), x[2]))
    for name, L in rep.items():
        o.write("\n%s: %d\n" % (name, len(L)))
        for x in L[:40]: o.write("  " + " ".join(x[0]) + "  alpha %g  %s\n" % (x[1], " ".join(str(y) for y in x[2:])))
        if len(L) > 40: o.write("  ... (%d more)\n" % (len(L) - 40))
print(open(A.out + ".report.txt").read()[:4000])
