#!/usr/bin/env python3
"""Install the computed points (rates_out/<channel>.pts, from Compute_sg_points.jl) as rate tables in the dust
worktree's rate_sve (/data/dust/user/samwitte/Axion_SR_sg/src/rate_sve):

  <ch>_LvrHc_SG.dat   alpha, G_SI, G_SG, G_x, kappa   for every channel of sgx_channels.txt with computed points.
                      An SG table of the first campaign keeps its rows above the new grid (alpha > 1.02 max new alpha);
                      its rows inside the new grid are replaced (original saved in rate_sve/backup_pre_sgx/).
                      Regime-split channels (regime_split.txt) get the same table under both the BH and Inf label.
  <ch>_LvrHc_.dat     alpha, G_SI                     only for channels that have no SI table at all (the new ones).
  sgx_channels_ready.txt   the channels of sgx_channels.txt that now have an SG table (the run list).

    python3 merge_rates.py            # dry run: report only
    python3 merge_rates.py --apply
For a channel computed with several settings the row with the finest step is used at each alpha.
"""
import os, sys, math, shutil, collections
HERE = os.path.dirname(os.path.abspath(__file__))
RS = "/data/dust/user/samwitte/Axion_SR_sg/src/rate_sve/"
OUT = os.path.join(HERE, "rates_out")
BAK = RS + "backup_pre_sgx/"
apply = "--apply" in sys.argv
chans = [tuple(ln.split()[:4]) for ln in open(os.path.join(HERE, "sgx_channels.txt")) if ln.strip()]
split = set()
for ln in open(RS + "regime_split.txt"):
    p = ln.split()
    if len(p) >= 4 and not p[0].startswith("#"): split.add(tuple(p[:3]))
if apply: os.makedirs(BAK, exist_ok=True)
# SI tables written by an earlier merge (channels that had none): rewritten whenever more points exist
MAN = os.path.join(HERE, "sgx_si_written.txt")
mine = set(ln.strip() for ln in open(MAN)) if os.path.isfile(MAN) else set()
stat = collections.Counter(); bad = []; ready = []; mism = 0; npt = 0; kneg = 0
for c in chans:
    nm = "_".join(c); fsi = RS + nm + "_LvrHc_.dat"; fsg = RS + nm + "_LvrHc_SG.dat"; fp = os.path.join(OUT, nm + ".pts")
    have_si, have_sg = os.path.isfile(fsi) and nm not in mine, os.path.isfile(fsg)
    # a channel with an SI table from the paper archive stays in the run list even without an SG table (SI only)
    si_archive = os.path.isfile(fsi) and nm not in mine
    if not os.path.isfile(fp):
        stat["not computed (old SG table kept)" if have_sg else
             "not computed (archive SI table only: kept, SI only)" if si_archive else "not computed (dropped from the run list)"] += 1
        if have_sg or si_archive: ready.append(c)
        continue
    best = {}
    for ln in open(fp):
        p = ln.split()
        if len(p) < 9: continue
        al = round(float(p[0]), 10); hf = float(p[7])
        if al not in best or hf > best[al][7]: best[al] = [float(x) for x in p[:9]]
    rows = [best[k] for k in sorted(best)]
    ok = [r for r in rows if all(math.isfinite(x) for x in r[:5]) and r[1] > 0 and r[2] > 0 and abs(r[4]) <= 1.05]
    if len(ok) < len(rows): bad.append((nm, len(rows) - len(ok), len(rows)))
    mism += sum(1 for r in ok if (r[8] == 1) != (c[3] == "Inf")); npt += len(ok); kneg += sum(1 for r in ok if r[4] < 0)
    if len(ok) < 2:
        stat["too few valid points" + (" (old SG table kept)" if have_sg else " (archive SI table only: kept, SI only)" if si_archive else " (dropped)")] += 1
        if have_sg or si_archive: ready.append(c)
        continue
    old_hi = []
    if have_sg:
        src = BAK + nm + "_LvrHc_SG.dat" if os.path.isfile(BAK + nm + "_LvrHc_SG.dat") else fsg   # rerun: start from the original
        for ln in open(src):
            p = [float(x) for x in ln.split()[:5]]
            if len(p) == 5 and p[0] > 1.02 * ok[-1][0]: old_hi.append(p)
    stat["SG table merged into the first-campaign table" if have_sg else ("SG table new" if have_si else "SI + SG tables new")] += 1
    if apply:
        if have_sg and not os.path.isfile(BAK + nm + "_LvrHc_SG.dat"): shutil.copy2(fsg, BAK + nm + "_LvrHc_SG.dat")
        labels = [c[3]] + ([("Inf" if c[3] == "BH" else "BH")] if c[:3] in split else [])
        for lab in labels:        # write to a temporary name and rename, so an interrupted merge never leaves a cut table
            ft = RS + "_".join(c[:3] + (lab,)) + "_LvrHc_SG.dat"
            with open(ft + ".tmp", "w") as f:
                for r in ok + old_hi: f.write("\t".join(repr(x) for x in r[:5]) + "\n")
            os.replace(ft + ".tmp", ft)
        if not have_si:
            with open(fsi + ".tmp", "w") as f:
                for r in ok: f.write(f"{r[0]!r}\t{r[1]!r}\n")
            os.replace(fsi + ".tmp", fsi)
            mine.add(nm)
    ready.append(c)
for k, v in sorted(stat.items()): print(f"{v:7d}  {k}")
print(f"valid points {npt} (kappa < 0 at {kneg}); points whose energy-based regime differs from the channel label: {mism}")
print(f"channels with rejected points: {len(bad)}", bad[:12])
if apply:
    with open(MAN, "w") as f: f.write("".join(x + "\n" for x in sorted(mine)))
    with open(os.path.join(HERE, "sgx_channels_ready.txt"), "w") as f:
        for c in ready: f.write("    ".join(c) + " \n")
    print(f"wrote tables; sgx_channels_ready.txt has {len(ready)} of {len(chans)} channels")
