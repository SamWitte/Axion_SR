#!/usr/bin/env python3
"""Evolution runs with the fully relativistic self-gravity tables: the SG runs of submit_sgx_runs.py (set "full")
repeated with NL_SG_OVERLAY=<tables dir>, everything else identical.
Output: ../runs/<name>/gw_<model>/M_<mass>/<variant>/alpha_<a>/   (default name: exact)

    python3 submit_exact_runs.py [--name exact] [--tables tables] [--mass ...] [--alpha ...] [--variants ...] [--gw ...] [--rate_input <list>] [--no-submit]
"""
import glob, os, sys, subprocess
HERE = os.path.dirname(os.path.abspath(__file__)); BASE = os.path.dirname(HERE)
SRC = "/data/dust/user/samwitte/Axion_SR_sg/src"; RUN_DIR = os.path.join(SRC, "Nlevels_Runs")
JULIA = os.path.expanduser("~/julia-1.8.0/bin/julia")
SPIN_BH, TAU_MAX, FA_OFF, N_TIMES = 0.95, 1e9, 1e30, 20000
VARIANTS = {"SG_only": FA_OFF, "SG_SI_fa_1e16": 1e16, "SG_SI_fa_5e16": 5e16, "SG_SI_fa_1e17": 1e17}
# self-interaction only (--self_grav false; with NL_SI_FROM_SG=1 Gamma_SI is still read from the SG tables / overlay);
# not in the default list, select with --variants
SI_ONLY = {"SI_only_fa_1e16": 1e16, "SI_only_fa_5e16": 5e16, "SI_only_fa_1e17": 1e17}

def opt(name, default, n=1):
    if name not in sys.argv: return default
    i = sys.argv.index(name); vals = []
    for a in sys.argv[i + 1:]:
        if a.startswith("--"): break
        vals.append(a)
    return vals if n != 1 else vals[0]
fmt = lambda x: f"{x:g}"
st = lambda tok: tuple(int(x) for x in tok.split("-")) if "-" in tok else tuple(int(c) for c in tok)
def status(outdir):
    hits = glob.glob(os.path.join(outdir, "Time_*.dat"))
    if not hits: return "new"
    with open(hits[0], "rb") as f:
        f.seek(0, 2); f.seek(max(0, f.tell() - 200)); tail = f.read().split()
    return "done" if (float(tail[-1]) if tail else 0.0) >= 0.99 * TAU_MAX else "resume"

name = opt("--name", "exact"); tables = os.path.abspath(os.path.join(HERE, opt("--tables", "tables")))
masses = [float(m) for m in opt("--mass", ["10", "1e4"], n=0)]
alphas = [float(a) for a in opt("--alpha", ["0.05", "0.075", "0.1", "0.15", "0.2", "0.25"], n=0)]
variants = opt("--variants", list(VARIANTS), n=0); VARIANTS.update(SI_ONLY); gws = opt("--gw", ["nonrel", "rel"], n=0)
rin = opt("--rate_input", "sgx_channels.txt")          # channel list: file name in rate_sve/ or an absolute path
chans = [ln.split()[:4] for ln in open(rin if os.path.isabs(rin) else os.path.join(SRC, "rate_sve", rin)) if ln.strip()]
levels = sorted(set(st(s) for c in chans for s in c[:3]))
keep_levels = ",".join(f"{n}-{l}-{m}" for n, l, m in levels)
log_dir = os.path.join(BASE, "logs", "runs", name); os.makedirs(log_dir, exist_ok=True)
cmds = []
for gw in gws:
    for mass in masses:
        for alpha in alphas:
            for v in variants:
                outdir = os.path.join(BASE, "runs", name, f"gw_{gw}", f"M_{fmt(mass)}", v, f"alpha_{fmt(alpha)}")
                s = status(outdir)
                if s == "done": continue
                os.makedirs(outdir, exist_ok=True)
                log = os.path.join(log_dir, f"gw_{gw}_M{fmt(mass)}_{v}_alpha_{fmt(alpha)}.log")
                cmds.append((mass, f"{JULIA} run_Nlevels.jl --MassBH {fmt(mass)} --SpinBH {SPIN_BH} --f_a {fmt(VARIANTS[v])} --alpha {fmt(alpha)} "
                             f"--Nmax 45 --tau_max {fmt(TAU_MAX)} --gw_model {gw} --self_grav {'false' if v in SI_ONLY else 'true'} "
                             f"--rate_input {rin} --keep_levels \"{keep_levels}\" --gw_table_nmax sgx --regime_split true "
                             f"--solver radau --analytic_jac --timescale_saves false --maxiters 1e8 --n_times {N_TIMES} "
                             f"--outdir {outdir}" + (" --resume" if s == "resume" else "") + f" >> {log} 2>&1"))
cmds.sort(key=lambda t: t[0])
print(f"{len(levels)} levels, {len(chans)} channels, tables {tables} ({len(glob.glob(tables + '/*_SG.dat'))} files), {len(cmds)} runs to do")
if not cmds: sys.exit(0)
tf = os.path.join(HERE, f"run_tasks_{name}.sh")
open(tf, "w").write("\n".join(c for _, c in cmds) + "\n")
script = os.path.join(HERE, f"runs_{name}.sh")
open(script, "w").write(f"""#!/bin/bash
#SBATCH --partition=maxcpu --constraint="7402|75F3|7H12" --nodes=1 --time=2-00:00:00
#SBATCH --job-name=sgx_{name}
#SBATCH --output={log_dir}/runs_{name}_%j.out
source ~/.bashrc
cd {RUN_DIR}
export JULIA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export NL_MAX_WALL_SEC={int(1.9 * 24 * 3600)} NL_PROGRESS_SEC=1800 NL_SI_FROM_SG=1 NL_SG_OVERLAY={tables}
echo "node $(hostname), {len(cmds)} runs, start $(date)"
pids=(); while IFS= read -r line; do [ -z "$line" ] && continue; bash -c "$line" & pids+=($!); sleep 2; done < {tf}
rc=0; for p in "${{pids[@]}}"; do wait "$p" || {{ echo "run pid $p exit $?"; rc=1; }}; done
echo "end $(date) rc=$rc"; exit $rc
""")
os.chmod(script, 0o755)
if "--no-submit" not in sys.argv:
    print(subprocess.run(f"sbatch --parsable {script}", shell=True, capture_output=True, text=True).stdout.strip())
