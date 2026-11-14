#!/usr/bin/env python3
"""Evolution runs of the extended self-gravity system (requested 2026-11-12).

M_BH = 10 and 1e4 Msun, spin 0.95, alpha = 0.05 0.075 0.1 0.15 0.2 0.25, tau_max = 1e9 yr, five physics variants
(noSG_noSI, SG_only, SG + SI at f_a = 1e16, 5e16, 1e17 GeV; "no SI" = f_a 1e30 as in the earlier campaigns), each
with non-relativistic (gw_nonrel) and relativistic (gw_rel, Kerr/Teukolsky tables where they exist) GW emission.

Code: /data/dust/user/samwitte/Axion_SR_sg (main bac01b2 merged with relativistic-sg), run_Nlevels.jl with
  --keep_levels <the levels of the set>   only these levels are evolved (no level without scattering channels)
  --rate_input sgx_channels[_core].txt    the channel list (written by merge_rates.py / this script into rate_sve/)
  --gw_table_nmax sgx                     GW table for exactly these levels (build_gw_nr_custom.jl)
  --solver radau --analytic_jac --timescale_saves false --maxiters 1e8      (as the paper's tower runs)
  NL_SI_FROM_SG=1                         Gamma_SI from the SG tables (same alpha grid and step as Gamma_SG, kappa)
Sets:  full = every level of levels_sgx.txt (Table I + first-campaign levels + (n,4,4) to n = 20), Nmax 45
       p1   = core + the levels populated at alpha <= 0.2 in Fig. 2 with n <= 20 (60 levels), Nmax 20
       core = first-campaign levels + (n,4,4) to n = 20 only (34 levels, the p0 channels), Nmax 20

    python3 submit_sgx_runs.py [--set full|core] [--mass 10 1e4] [--alpha ...] [--variants ...] [--gw nonrel rel]
                               [--per-node 60] [--threads 1] [--no-submit] [--dependency JOBID]
Rerunning skips finished runs (Time_ file reaches 0.99 tau_max) and resumes stopped ones (--resume).
Output: runs/<set>/gw_<model>/M_<mass>/<variant>/alpha_<a>/ ; logs: logs/runs/<set>/
"""
import glob, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from states_lib import SG21, N44, FIG_A02, fmt as sfmt

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = "/data/dust/user/samwitte/Axion_SR_sg/src"
RUN_DIR = os.path.join(SRC, "Nlevels_Runs")
JULIA = os.path.expanduser("~/julia-1.8.0/bin/julia")
SPIN_BH, TAU_MAX, FA_OFF = 0.95, 1e9, 1e30
MASSES = [10.0, 1e4]
ALPHAS = [0.05, 0.075, 0.1, 0.15, 0.2, 0.25]
VARIANTS = {"noSG_noSI": (FA_OFF, False), "SG_only": (FA_OFF, True), "SG_SI_fa_1e16": (1e16, True),
            "SG_SI_fa_5e16": (5e16, True), "SG_SI_fa_1e17": (1e17, True)}
GW_MODELS = ["nonrel", "rel"]
N_TIMES = 20000
WALL_SEC_CAP = int(6.5 * 24 * 3600)


def opt(name, default, n=1):
    if name not in sys.argv:
        return default
    i = sys.argv.index(name); vals = []
    for a in sys.argv[i + 1:]:
        if a.startswith("--"): break
        vals.append(a)
    return vals if n != 1 else vals[0]


def fmt(x): return f"{x:g}"


def st(tok): return tuple(int(x) for x in tok.split("-")) if "-" in tok else tuple(int(c) for c in tok)


def status(outdir):
    hits = glob.glob(os.path.join(outdir, "Time_*.dat"))
    if not hits: return "new"
    with open(hits[0], "rb") as f:
        f.seek(0, 2); f.seek(max(0, f.tell() - 200)); tail = f.read().split()
    return "done" if (float(tail[-1]) if tail else 0.0) >= 0.99 * TAU_MAX else "resume"


sset = opt("--set", "full")
masses = [float(m) for m in opt("--mass", [fmt(m) for m in MASSES], n=0)]
alphas = [float(a) for a in opt("--alpha", [fmt(a) for a in ALPHAS], n=0)]
variants = opt("--variants", list(VARIANTS), n=0)
gws = opt("--gw", GW_MODELS, n=0)
per_node = int(opt("--per-node", "60")); threads = int(opt("--threads", "1"))
dependency = opt("--dependency", None); submit = "--no-submit" not in sys.argv

# ---- level set and channel list of this set
ready = os.path.join(HERE, "sgx_channels_ready.txt")
chans = [ln.split()[:4] for ln in open(ready) if ln.strip()]
if sset in ("core", "p1"):
    keep = set(SG21) | set(N44)
    if sset == "p1":      # + the levels populated at alpha <= 0.2 in Fig. 2 of the draft with n <= 20 (priority-1 channels)
        keep |= set(s for s in FIG_A02 if s[0] <= 20)
    chans = [c for c in chans if all(st(s) in keep for s in c[:3])]
    nmax, rin = 20, f"sgx_channels_{sset}.txt"
else:
    nmax, rin = 45, "sgx_channels.txt"
levels = sorted(set(st(s) for c in chans for s in c[:3]))
with open(os.path.join(SRC, "rate_sve", rin), "w") as f:
    for c in chans: f.write("    ".join(c) + " \n")
keep_levels = ",".join(f"{n}-{l}-{m}" for n, l, m in levels)
print(f"set {sset}: {len(levels)} levels, {len(chans)} channels -> rate_sve/{rin}, Nmax {nmax}")

log_dir = os.path.join(HERE, "logs", "runs", sset); os.makedirs(log_dir, exist_ok=True)
cmds = []
for gw in gws:
    for mass in masses:
        for alpha in alphas:
            for v in variants:
                fa, sg = VARIANTS[v]
                outdir = os.path.join(HERE, "runs", sset, f"gw_{gw}", f"M_{fmt(mass)}", v, f"alpha_{fmt(alpha)}")
                s = status(outdir)
                if s == "done":
                    print(f"done, skipping   gw={gw:6s} M={fmt(mass):5s} {v:14s} alpha={fmt(alpha)}"); continue
                os.makedirs(outdir, exist_ok=True)
                log = os.path.join(log_dir, f"gw_{gw}_M{fmt(mass)}_{v}_alpha_{fmt(alpha)}.log")
                cmds.append((mass, f"{JULIA} run_Nlevels.jl --MassBH {fmt(mass)} --SpinBH {SPIN_BH} --f_a {fmt(fa)} --alpha {fmt(alpha)} "
                             f"--Nmax {nmax} --tau_max {fmt(TAU_MAX)} --gw_model {gw} --self_grav {'true' if sg else 'false'} "
                             f"--rate_input {rin} --keep_levels \"{keep_levels}\" --gw_table_nmax sgx --regime_split true "
                             f"--solver radau --analytic_jac --timescale_saves false --maxiters 1e8 --n_times {N_TIMES} "
                             f"--outdir {outdir}" + (" --resume" if s == "resume" else "") + f" >> {log} 2>&1"))
# heavy (10 Msun) runs first within each node, spread evenly
cmds.sort(key=lambda t: t[0])
ntask = max(1, -(-len(cmds) // per_node))
task_dir = os.path.join(HERE, "run_tasks", sset); os.makedirs(task_dir, exist_ok=True)
tfs = []
for k in range(ntask):
    part = [c for i, (_, c) in enumerate(cmds) if i % ntask == k]
    if not part: continue
    tf = os.path.join(task_dir, f"task_{k:02d}.sh")
    with open(tf, "w") as f: f.write("\n".join(part) + "\n")
    tfs.append(tf)
print(f"{len(cmds)} runs to do in {len(tfs)} node task(s) of <= {per_node}")
if not tfs:
    sys.exit(0)
script = os.path.join(HERE, f"sgx_runs_{sset}.sh")
with open(script, "w") as f:
    f.write(f"""#!/bin/bash
#SBATCH --partition=maxcpu --constraint="7402|75F3|7H12" --nodes=1 --time=7-00:00:00
#SBATCH --job-name=sgx_{sset} --array=0-{len(tfs) - 1}
#SBATCH --output={log_dir}/sgx_{sset}_%A_%a.out
source ~/.bashrc
cd {RUN_DIR}
export JULIA_NUM_THREADS=1 OMP_NUM_THREADS={threads} OPENBLAS_NUM_THREADS={threads}
export NL_MAX_WALL_SEC={WALL_SEC_CAP} NL_PROGRESS_SEC=1800 NL_SI_FROM_SG=1
tasks=({" ".join(tfs)})
task=${{tasks[$SLURM_ARRAY_TASK_ID]}}
echo "node $(hostname), task $task, $(wc -l < "$task") runs, start $(date)"
pids=(); while IFS= read -r line; do [ -z "$line" ] && continue; bash -c "$line" & pids+=($!); sleep 2; done < "$task"
rc=0; for p in "${{pids[@]}}"; do wait "$p" || {{ echo "run pid $p exit $?"; rc=1; }}; done
echo "end $(date) rc=$rc"; exit $rc
""")
os.chmod(script, 0o755)
print(f"-> {script}")
if submit:
    import subprocess
    dep = f"--dependency=afterany:{dependency} " if dependency else ""
    jid = subprocess.run(f"sbatch --parsable {dep}{script}", shell=True, capture_output=True, text=True).stdout.strip().split(";")[0]
    print(f"submitted run array {jid}")
    if jid:      # summary table + figures once the runs are over
        os.system(f"sbatch --dependency=afterany:{jid} --export=ALL,SGX_SET={sset} {os.path.join(HERE, 'post.sh')}")
