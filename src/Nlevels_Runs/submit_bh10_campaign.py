#!/usr/bin/env python3
"""Submit full BH_10 production grid with master rates (Fresh_cpy).

Grid: fa × alpha × Nmax = 4 × 8 × 8 = 256 jobs (one Nmax per job).
Skips if NL_output_Nmax_N.tar.gz already exists and is non-tiny.

Wall-clock checkpoints: Julia stops via NL_MAX_WALL_SEC (default 4.5d) so it
can write Time_/States_/… before the 5d Slurm limit. Incomplete runs skip
finish/prune and leave loose .dat for a later --resume hop.

Usage on cluster:
  cd .../Fresh_cpy/Axion_SR/src/Nlevels_Runs
  python3 submit_bh10_campaign.py --no-submit   # dry create
  python3 submit_bh10_campaign.py               # cold start + sbatch
  python3 submit_bh10_campaign.py --resume \\
      --fa 1e16 --alpha 1.3 --nmax 15 18      # continue from .dat
"""
from __future__ import annotations

import argparse
import itertools
import os
import subprocess
import sys

BASE_DIR = os.environ.get(
    "NL_BASE_DIR",
    os.path.dirname(os.path.abspath(__file__)),
)
RUN_TAG = os.environ.get("RUN_TAG", "run_nodrag").strip()
OUTPUT_DIR = (
    os.path.join(BASE_DIR, "output", RUN_TAG)
    if RUN_TAG
    else os.path.join(BASE_DIR, "output")
)
SCRIPT_DIR = os.path.join(BASE_DIR, "scripts_bh10")
LOG_DIR = os.path.join(BASE_DIR, "logs")
JULIA_BIN = os.environ.get(
    "JULIA_BIN", "/groups/astro/spieksma/julia-1.8.0/bin/julia"
)

MassBH = 10.0
SpinBH = 0.95
tau_max = 5e7
FA_VALS = [1e18, 1e16, 1e14, 1e12]
ALPHA_VALS = [0.05, 0.2, 0.4, 0.6, 0.8, 1.0, 1.3, 1.6]
NMAX_LIST = [3, 4, 5, 6, 7, 8, 15, 18]

PARTITION = os.environ.get("NL_PARTITION", "astro3_long")
WALLCLOCK = os.environ.get("NL_WALLCLOCK", "5-00:00:00")
# Shell timeout must exceed NL_MAX_WALL_SEC so Julia can write after ODE stop.
TIMEOUT = os.environ.get("NL_TIMEOUT", "118h")
# Default 4.5 days — shorter than astro3_long's 5d Slurm wall.
NL_MAX_WALL_SEC = os.environ.get("NL_MAX_WALL_SEC", str(int(4.5 * 24 * 3600)))


def fa_tag(fa: float) -> str:
    import math

    return f"fa_1e{int(round(math.log10(fa)))}"


def alpha_tag(a: float) -> str:
    return f"alpha_{a:g}"


def resources(nmax: int) -> tuple[int, int]:
    if nmax == 18:
        return 10, 96
    if nmax == 15:
        return 10, 80
    if nmax >= 8:
        return 10, 64
    return 10, 48


def already_done(outdir: str, nmax: int) -> bool:
    path = os.path.join(outdir, f"NL_output_Nmax_{nmax}.tar.gz")
    marker = os.path.join(outdir, f"NL_output_Nmax_{nmax}.pruned")
    return (
        os.path.isfile(path)
        and os.path.getsize(path) >= 1024
        and os.path.isfile(marker)
    )


def write_script(fa: float, alpha: float, nmax: int, *, resume: bool) -> str:
    outdir = os.path.join(OUTPUT_DIR, "BH_10", fa_tag(fa), alpha_tag(alpha))
    tag = f"BH_10_{fa_tag(fa)}_{alpha_tag(alpha)}_Nmax_{nmax}"
    if resume:
        tag = f"{tag}_resume"
    script = os.path.join(SCRIPT_DIR, f"NL_{tag}.sh")
    cpus, mem = resources(nmax)
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(LOG_DIR, exist_ok=True)
    resume_arg = " \\\n  --resume" if resume else ""
    if resume:
        prep = f"""# Resume: keep existing Time_/States_/… checkpoint for Nmax=$NMAX
TIMEFILE=$(ls "$OUTDIR"/Time_*Nmax_${{NMAX}}.dat 2>/dev/null | head -1 || true)
if [[ -z "${{TIMEFILE}}" ]]; then
  echo "ERROR: --resume but no Time_*Nmax_${{NMAX}}.dat in $OUTDIR"
  exit 1
fi
echo "Resuming from: $TIMEFILE"
# Drop any prior archive/marker so a completed hop can re-finish
rm -f "$ARCHIVE" "$MARKER" || true
"""
    else:
        prep = """# Cold start: drop any stale loose outputs / partial archive for this Nmax
rm -f "$OUTDIR"/{Time,Spin,States,Modes,MassBH}_*Nmax_${NMAX}.dat \\
      "$OUTDIR"/{Time,Spin,States,Modes,MassBH}_*Nmax_${NMAX}.dat.pruned \\
      "$ARCHIVE" "$MARKER" || true
"""

    body = f"""#!/bin/bash
#SBATCH --job-name=NL_M10
#SBATCH --time={WALLCLOCK}
#SBATCH --partition={PARTITION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --threads-per-core=1
#SBATCH --mem={mem}G
#SBATCH --output={LOG_DIR}/NL_{tag}_%j.out
#SBATCH --error={LOG_DIR}/NL_{tag}_%j.err

set -euo pipefail
unset LD_PRELOAD
unset LD_LIBRARY_PATH
cd {BASE_DIR}
mkdir -p logs
export NL_BASE_DIR={BASE_DIR}
export JULIA_NUM_THREADS={cpus}
export OMP_NUM_THREADS={cpus}
export PYTHONPATH={BASE_DIR}:${{PYTHONPATH:-}}
export PYTHON_BIN=/lustre/hpc/astro/spieksma/miniforge3/bin/python3
export NL_MAX_WALL_SEC={NL_MAX_WALL_SEC}

OUTDIR={outdir}
NMAX={nmax}
TAU_MAX={tau_max:g}
ARCHIVE="$OUTDIR/NL_output_Nmax_${{NMAX}}.tar.gz"
MARKER="$OUTDIR/NL_output_Nmax_${{NMAX}}.pruned"
if [[ -f "$ARCHIVE" && -s "$ARCHIVE" && -f "$MARKER" ]]; then
  echo "Already complete (pruned+compressed): $ARCHIVE"
  exit 0
fi

{prep}
echo "=== BH_10 fa={fa:g} alpha={alpha:g} Nmax=$NMAX resume={str(resume).lower()} NL_MAX_WALL_SEC=$NL_MAX_WALL_SEC ==="
set +e
timeout {TIMEOUT} srun --exclusive {JULIA_BIN} run_Nlevels.jl \\
  --MassBH 10 --SpinBH {SpinBH} \\
  --f_a {fa:g} --alpha {alpha:g} \\
  --Nmax $NMAX --tau_max $TAU_MAX \\
  --outdir "$OUTDIR"{resume_arg}
status=$?
set -e
if [ "$status" -eq 124 ]; then echo "TIMEOUT Nmax=$NMAX"; exit 124; fi
if [ "$status" -ne 0 ]; then echo "Julia failed status=$status"; exit "$status"; fi
wait

TIMEFILE=$(ls "$OUTDIR"/Time_*Nmax_${{NMAX}}.dat 2>/dev/null | head -1 || true)
if [[ -z "${{TIMEFILE}}" ]]; then
  echo "ERROR: no Time_*Nmax_${{NMAX}}.dat after Julia"
  exit 1
fi

# If wall-clock stopped early, leave loose .dat for a later --resume hop.
incomplete=$("$PYTHON_BIN" -c "
import numpy as np
t = np.atleast_1d(np.loadtxt('${{TIMEFILE}}'))
print(1 if float(t[-1]) < 0.99 * float('${{TAU_MAX}}') else 0)
")
if [[ "$incomplete" == "1" ]]; then
  echo "INCOMPLETE checkpoint (t < 0.99*tau_max); skipping finish/prune"
  echo "Re-submit with: python3 submit_bh10_campaign.py --resume --fa {fa:g} --alpha {alpha:g} --nmax $NMAX"
  ls -lh "$OUTDIR"/{{Time,Spin,States,Modes,MassBH}}_*Nmax_${{NMAX}}.dat
  echo CHECKPOINT_SAVED
  exit 0
fi

echo "=== prune + compress Nmax=$NMAX ==="
bash {BASE_DIR}/finish_nmax_output.sh "$OUTDIR" $NMAX
if [[ ! -f "$ARCHIVE" || ! -s "$ARCHIVE" ]]; then
  echo "ERROR: missing archive after finish_nmax_output: $ARCHIVE"
  exit 1
fi
if [[ ! -f "$MARKER" ]]; then
  echo "ERROR: missing prune marker: $MARKER"
  exit 1
fi
ls -lh "$ARCHIVE" "$MARKER"
echo DONE
"""
    with open(script, "w") as f:
        f.write(body)
    os.chmod(script, 0o755)
    return script


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-submit", action="store_true")
    ap.add_argument("--resume", action="store_true",
                    help="Continue from existing Time_/States_ checkpoint (.dat)")
    ap.add_argument("--nmax", type=int, nargs="*", default=None)
    ap.add_argument("--fa", type=float, nargs="*", default=None)
    ap.add_argument("--alpha", type=float, nargs="*", default=None)
    args = ap.parse_args()

    fas = args.fa or FA_VALS
    alphas = args.alpha or ALPHA_VALS
    nmaxs = args.nmax or NMAX_LIST

    os.makedirs(SCRIPT_DIR, exist_ok=True)
    scripts = []
    skipped = 0
    for fa, alpha, nmax in itertools.product(fas, alphas, nmaxs):
        outdir = os.path.join(OUTPUT_DIR, "BH_10", fa_tag(fa), alpha_tag(alpha))
        if already_done(outdir, nmax):
            skipped += 1
            continue
        if args.resume:
            # Need a Time_ checkpoint to resume; skip holes with nothing.
            import glob as _glob

            hits = _glob.glob(os.path.join(outdir, f"Time_*Nmax_{nmax}.dat"))
            if not hits:
                print(f"skip resume (no checkpoint): {fa_tag(fa)}/{alpha_tag(alpha)} Nmax={nmax}")
                skipped += 1
                continue
        scripts.append(write_script(fa, alpha, nmax, resume=args.resume))

    print(f"OUTPUT_DIR={OUTPUT_DIR}")
    print(f"NL_MAX_WALL_SEC={NL_MAX_WALL_SEC} TIMEOUT={TIMEOUT} WALLCLOCK={WALLCLOCK}")
    print(f"resume={args.resume} scripts_created={len(scripts)} skipped={skipped}")
    if args.no_submit:
        print("(--no-submit) not submitting")
        return

    job_ids = []
    for script in scripts:
        out = subprocess.check_output(["sbatch", "--parsable", script], text=True)
        job_ids.append(out.strip())
    print(f"submitted={len(job_ids)}")
    if job_ids:
        print("first_ids:", ", ".join(job_ids[:5]))
        print("last_ids:", ", ".join(job_ids[-3:]))


if __name__ == "__main__":
    main()
