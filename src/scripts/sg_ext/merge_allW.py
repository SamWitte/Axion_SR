#!/usr/bin/env python3
"""Rate tables of the self-gravity campaign in their corrected form (19 Nov 2026), as an overlay directory:
  (1) points in the "to infinity" regime recomputed with all outgoing multipoles and the corrected outgoing condition
      and outer radius (allW/rates_out, gf_radial inf_all_l + inf_wkb);
  (2) every row converted to the Klein-Gordon normalisation of the three cloud modes, rate -> rate / (Q_1 Q_2 Q_3), with
      the charges of ../kg_norm/kg_levels.txt (kg_levels.jl; the points themselves were computed with int |R|^2 r^2 dr = 1).
Channels with new points: table built as merge_rates.py does, from the stored points (../rates_out) with every alpha
that has a valid new point replaced by it.  Channels without new points (absorption by the hole): the production table
as it was before these tables were installed (rate_sve/backup_pre_v1120, i.e. still in the flat normalisation) times
the factor.  Do not point PROD at rate_sve itself after the installation (22 Nov 2026): the factor would be applied twice.
Rows off the alpha grid (rows of a first-campaign table above the grid) are copied
unchanged: they keep the old normalisation, single multipole and outgoing condition.
    <ch>_LvrHc_SG.dat   alpha, G_SI, G_SG, G_x, kappa   (also under the other label for regime-split channels)
    <ch>_LvrHc_.dat     alpha, G_SI   for the channels whose SI table was written by merge_rates.py (sgx_si_written.txt)
Tables go to allW/tables (NL_SG_OVERLAY=<dir> uses them without touching rate_sve).

    python3 merge_allW.py            report only
    python3 merge_allW.py --apply    write allW/tables and allW/merge_report.txt
"""
import os, sys, math, collections
HERE = os.path.dirname(os.path.abspath(__file__)); X = os.path.dirname(HERE)
RS = "/data/dust/user/samwitte/Axion_SR_sg/src/rate_sve/"; BAK = RS + "backup_pre_sgx/"
PROD = RS + "backup_pre_v1120/"   # production tables before the installation of 22 Nov 2026 (flat normalisation)
OLD = os.path.join(X, "rates_out"); NEW = os.path.join(HERE, "rates_out"); OUT = os.path.join(HERE, "tables")
apply = "--apply" in sys.argv
chans = [tuple(l.split()[:4]) for l in open(RS + "sgx_channels.txt") if l.strip()]
split = set()
for l in open(RS + "regime_split.txt"):
    p = l.split()
    if len(p) >= 4 and not p[0].startswith("#"): split.add(tuple(p[:3]))
mine = set(l.strip() for l in open(os.path.join(X, "sgx_si_written.txt")))
Q = {}
for l in open(os.path.join(X, "kg_norm", "kg_levels.txt")):
    p = l.split(); Q[(p[0], round(float(p[1]), 10))] = float(p[2])
if apply: os.makedirs(OUT, exist_ok=True)


def kgfac(c, al):
    try: return 1.0 / (Q[(c[0], round(al, 10))] * Q[(c[1], round(al, 10))] * Q[(c[2], round(al, 10))])
    except KeyError: return None


def best_rows(path):
    best = {}
    if os.path.isfile(path):
        for l in open(path):
            p = l.split()
            if len(p) < 9: continue
            try: r = [float(x) for x in p[:9]]
            except ValueError: continue
            al = round(r[0], 10)
            if al not in best or r[7] > best[al][7]: best[al] = r
    return best


def write(path, rows, ncol):
    with open(path + ".tmp", "w") as f:
        for r in rows: f.write("\t".join(repr(x) for x in r[:ncol]) + "\n")
    os.replace(path + ".tmp", path)


def scaled(c, r, ncol):
    f = kgfac(c, r[0])
    if f is None: nrow["rows off the alpha grid (copied unchanged)"] += 1; return list(r[:ncol])
    nrow["rows converted to the Klein-Gordon normalisation"] += 1; fmin[0] = min(fmin[0], f)
    return [r[0]] + [x * f for x in r[1:min(ncol, 4)]] + list(r[4:ncol])


valid = lambda r: all(math.isfinite(x) for x in r[:5]) and r[1] > 0 and r[2] > 0 and abs(r[4]) <= 1.05
stat = collections.Counter(); nrow = collections.Counter(); kap = collections.Counter(); big = []; fmin = [1.0]
hist = {"SI": collections.Counter(), "SG": collections.Counter()}
bins = [(0, 1e-2, "< 0.01"), (1e-2, 0.1, "0.01 - 0.1"), (0.1, 0.5, "0.1 - 0.5"), (0.5, 0.9, "0.5 - 0.9"), (0.9, 1.1, "0.9 - 1.1"), (1.1, 2, "1.1 - 2"), (2, 10, "2 - 10"), (10, 1e99, "> 10")]
for c in chans:
    nm = "_".join(c); new = best_rows(os.path.join(NEW, nm + ".pts"))
    labels = [c[3]] + ([("Inf" if c[3] == "BH" else "BH")] if c[:3] in split else [])
    if not new:
        # no recomputed points: production table times the normalisation factor
        done = False
        for lab in labels:
            src = PROD + "_".join(c[:3] + (lab,)) + "_LvrHc_SG.dat"
            if not os.path.isfile(src): continue
            rows = [[float(x) for x in l.split()[:5]] for l in open(src) if len(l.split()) >= 5]
            out_rows = [scaled(c, r, 5) for r in rows]; done = True
            if apply: write(os.path.join(OUT, "_".join(c[:3] + (lab,)) + "_LvrHc_SG.dat"), out_rows, 5)
        stat["no new points: production table x normalisation" if done else "no new points and no production table (skipped)"] += 1
        if done and nm in mine and os.path.isfile(PROD + nm + "_LvrHc_.dat"):
            rows = [[float(x) for x in l.split()[:2]] for l in open(PROD + nm + "_LvrHc_.dat") if len(l.split()) >= 2]
            if apply: write(os.path.join(OUT, nm + "_LvrHc_.dat"), [scaled(c, r, 2) for r in rows], 2)
            stat["SI table rebuilt (channels without an archive SI table)"] += 1
        continue
    old = best_rows(os.path.join(OLD, nm + ".pts"))
    rows = dict(old)
    for al, r in new.items():
        if valid(r): rows[al] = r; nrow["rows replaced by a new point"] += 1
        else: nrow["new point invalid (stored point kept)"] += 1
    nrow["rows of these channels without a new point"] += sum(1 for al in rows if al not in new)
    ok = [rows[k] for k in sorted(rows) if valid(rows[k])]
    if len(ok) < 2: stat["too few valid points (skipped)"] += 1; continue
    for r in ok:
        if round(r[0], 10) in new:
            a = abs(r[4]); kap["> 0.99" if a > 0.99 else "0.9 - 0.99" if a > 0.9 else "0.5 - 0.9" if a > 0.5 else "< 0.5"] += 1
            o = old.get(round(r[0], 10))
            if o is not None and valid(o):
                for j, lab in ((1, "SI"), (2, "SG")):
                    q = r[j] / o[j]
                    hist[lab][next(b[2] for b in bins if b[0] <= q < b[1])] += 1
                    if q > 2 or q < 0.5: big.append((q, lab, nm, r[0]))
    old_hi = []
    src = BAK + nm + "_LvrHc_SG.dat" if os.path.isfile(BAK + nm + "_LvrHc_SG.dat") else None
    if src:
        for l in open(src):
            p = [float(x) for x in l.split()[:5]]
            if len(p) == 5 and p[0] > 1.02 * ok[-1][0]: old_hi.append(p)
    stat["table rebuilt from points (new points in the regime to infinity)"] += 1
    out_rows = [scaled(c, r, 5) for r in ok] + [list(p) for p in old_hi]; nrow["rows off the alpha grid (copied unchanged)"] += len(old_hi)
    if apply:
        for lab in labels: write(os.path.join(OUT, "_".join(c[:3] + (lab,)) + "_LvrHc_SG.dat"), out_rows, 5)
        if nm in mine:
            write(os.path.join(OUT, nm + "_LvrHc_.dat"), [r[:2] for r in out_rows[:len(ok)]], 2); stat["SI table rebuilt (channels without an archive SI table)"] += 1
out = [f"{v:8d}  {k}" for k, v in sorted(stat.items())] + [f"{v:8d}  {k}" for k, v in sorted(nrow.items())]
out.append(f"smallest normalisation factor 1 / (Q1 Q2 Q3): {fmin[0]:.4f}")
out.append("kappa of the new points: " + ", ".join(f"|kappa| {k}: {v}" for k, v in sorted(kap.items())))
for lab in ("SI", "SG"):
    tot = sum(hist[lab].values())
    out.append(f"{lab} rate, new point / stored point (before the normalisation factor), {tot} points: " + ", ".join(f"{b[2]}: {hist[lab][b[2]]} ({100 * hist[lab][b[2]] / max(tot, 1):.1f} %)" for b in bins))
big.sort(key=lambda t: -abs(math.log(t[0])))
out.append(f"points with a rate changed by more than 2x: {len(big)}; the largest:")
out += [f"      {lab} x{q:.4g}  {nm}  alpha {al:g}" for q, lab, nm, al in big[:30]]
print("\n".join(out))
if apply: open(os.path.join(HERE, "merge_report.txt"), "w").write("\n".join(out) + "\n")
