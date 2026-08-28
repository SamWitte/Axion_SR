#!/usr/bin/env python3
"""
Generate and submit a Slurm job array to compute rates for all entries in
load_rate_input_Nmax_18.dat.

Each array task covers one (S1, S2, S3, S4) combination and uses NCPUS cores.
Jobs whose output .dat already exists are skipped at submission time (to keep
the array compact) and also at runtime as a safety net.
"""

import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_DIR    = os.path.dirname(SCRIPT_DIR)          # .../src
_input_arg = next((sys.argv[i+1] for i, a in enumerate(sys.argv) if a == "--input" and i+1 < len(sys.argv)), None)
INPUT_FILE = os.path.abspath(_input_arg) if _input_arg else os.path.join(SCRIPT_DIR, "load_rate_input_Nmax_18.dat")
RATE_DIR   = os.path.join(SRC_DIR, "rate_sve")
LOG_DIR    = os.path.join(SCRIPT_DIR, "rate_logs")

NCPUS      = 2
MEM_GB     = 16
WALLCLOCK  = "7-00:00:00"
PARTITION  = "long"

FTAG       = "_LvrHc_"
ALPHA_MIN  = 0.05
ALPHA_PTS  = 14
RUN_LEAVER = True
CHECK_ERR  = False
SUBMIT     = "--no-submit" not in sys.argv

# ---------------------------------------------------------------------------
# Read input and filter already-computed entries
# ---------------------------------------------------------------------------
with open(INPUT_FILE) as fh:
    all_entries = [l.split() for l in fh if l.strip()]

pending = [
    (s1, s2, s3, s4) for s1, s2, s3, s4 in all_entries
    if not os.path.isfile(os.path.join(RATE_DIR, f"{s1}_{s2}_{s3}_{s4}{FTAG}.dat"))
]

n_total   = len(all_entries)
n_done    = n_total - len(pending)
n_pending = len(pending)

print(f"Input entries : {n_total}")
print(f"Already done  : {n_done}")
print(f"To compute    : {n_pending}")

if n_pending == 0:
    print("Nothing to do.")
    raise SystemExit(0)

# ---------------------------------------------------------------------------
# Write a compact input file containing only the pending entries
# (keeps the array indices small and avoids empty tasks)
# ---------------------------------------------------------------------------
os.makedirs(LOG_DIR, exist_ok=True)

pending_file = os.path.join(SCRIPT_DIR, "load_rate_pending.dat")
with open(pending_file, "w") as fh:
    for row in pending:
        fh.write("  ".join(row) + "\n")

print(f"Pending list  : {pending_file}  ({n_pending} entries)")

# ---------------------------------------------------------------------------
# Write the array job script
# ---------------------------------------------------------------------------
array_script = os.path.join(SCRIPT_DIR, "Script_Rate_array.sh")

rl_str  = "true"  if RUN_LEAVER else "false"
ce_str  = "true"  if CHECK_ERR  else "false"

with open(array_script, "w") as f:
    f.write("#!/bin/bash\n")
    f.write(f"#SBATCH --job-name=Rates\n")
    f.write(f"#SBATCH --partition={PARTITION}\n")
    f.write(f"#SBATCH --ntasks=1\n")
    f.write(f"#SBATCH --cpus-per-task={NCPUS}\n")
    f.write(f"#SBATCH --mem={MEM_GB}G\n")
    f.write(f"#SBATCH --time={WALLCLOCK}\n")
    f.write(f"#SBATCH --array=0-{n_pending - 1}\n")
    f.write(f"#SBATCH --output={LOG_DIR}/rate_%a.log\n")
    f.write("\n")
    f.write("source /home/phys2564/.bashrc\n")
    f.write(f"cd {SRC_DIR}\n")
    f.write("\n")
    f.write(f"# Pick this task's entry from the pending list (0-indexed)\n")
    f.write(f"read S1 S2 S3 S4 < <(sed -n \"$((SLURM_ARRAY_TASK_ID + 1))p\" {pending_file})\n")
    f.write("\n")
    f.write("# Safety-net skip if output appeared since submission\n")
    f.write(f'outfile="{RATE_DIR}/${{S1}}_${{S2}}_${{S3}}_${{S4}}{FTAG}.dat"\n')
    f.write('if [ -f "$outfile" ]; then\n')
    f.write('    echo "Already exists, skipping: $outfile"\n')
    f.write('    exit 0\n')
    f.write('fi\n')
    f.write("\n")
    f.write(f"export OMP_NUM_THREADS={NCPUS}\n")
    f.write('echo "Running: S1=$S1  S2=$S2  S3=$S3  S4=$S4"\n')
    f.write("julia Compute_all_rates.jl \\\n")
    f.write(f"    --alpha_min {ALPHA_MIN} \\\n")
    f.write(f"    --alpha_pts {ALPHA_PTS} \\\n")
    f.write( "    --S1 $S1 --S2 $S2 --S3 $S3 --S4 $S4 \\\n")
    f.write(f'    --ftag "{FTAG}" \\\n')
    f.write(f"    --run_leaver {rl_str} \\\n")
    f.write(f"    --check_err {ce_str}\n")

os.chmod(array_script, 0o755)
print(f"Written       : {array_script}")

# ---------------------------------------------------------------------------
# Submit
# ---------------------------------------------------------------------------
if SUBMIT:
    ret = os.system(f"sbatch {array_script}")
    if ret != 0:
        print(f"WARNING: sbatch returned exit code {ret}")
else:
    print(f"Skipped submission (--no-submit).  To submit manually:")
    print(f"  sbatch {array_script}")
