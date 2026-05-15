"""
submit_scan_jobs.py
-------------------
Discover every (BH, fa, alpha) directory under ./output/ and submit one SLURM
job per directory to run scan_active_states.py --alpha-dir <dir>.

Each job writes a JSON result to ./scan_results/<bh>_<fa>_<alpha>.json.
After all scan jobs finish, a final aggregation job runs
scan_active_states.py --aggregate ./scan_results/ to print the combined table.

Usage:
  python submit_scan_jobs.py            # generate scripts AND submit
  python submit_scan_jobs.py --no-submit  # generate scripts only

Scripts  -> ./scan_scripts/
Logs     -> ./scan_logs/
Results  -> ./scan_results/
"""

import os
import sys
import re

SUBMIT    = "--no-submit" not in sys.argv
BASE_DIR  = os.path.dirname(os.path.abspath(__file__))
OUT_ROOT  = os.path.join(BASE_DIR, "output")
RES_DIR   = os.path.join(BASE_DIR, "scan_results")
SCRPT_DIR = os.path.join(BASE_DIR, "scan_scripts")
LOG_DIR   = os.path.join(BASE_DIR, "scan_logs")

# SLURM resources per scan job
PARTITION       = "short"
MEM_GB          = 8
WALLCLOCK       = "2:00:00"
AGG_MEM_GB      = 2
AGG_WALLCLOCK   = "00:30:00"

SERVER_BASE = "/data/phys-axion-ns/phys2564/Axion_SR/src/Nlevels_Runs"
BASHRC      = "/home/phys2564/.bashrc"

BH_MASSES = ["BH_1e4", "BH_1e8"]

os.makedirs(SCRPT_DIR, exist_ok=True)
os.makedirs(LOG_DIR,   exist_ok=True)
os.makedirs(RES_DIR,   exist_ok=True)


def alpha_sort_key(name):
    m = re.search(r"alpha_(\d+\.?\d*)", name)
    return float(m.group(1)) if m else 0.0

def fa_sort_key(name):
    m = re.search(r"fa_1e(\d+)", name)
    return int(m.group(1)) if m else 0


def find_alpha_dirs(out_root, bh_masses):
    """Return sorted list of (bh, fa, alpha, alpha_dir) tuples."""
    tasks = []
    for bh in bh_masses:
        bh_dir = os.path.join(out_root, bh)
        if not os.path.isdir(bh_dir):
            continue
        fa_list = sorted(
            [d for d in os.listdir(bh_dir) if os.path.isdir(os.path.join(bh_dir, d))],
            key=fa_sort_key,
        )
        for fa in fa_list:
            fa_path = os.path.join(bh_dir, fa)
            alpha_list = sorted(
                [d for d in os.listdir(fa_path)
                 if d.startswith("alpha_") and os.path.isdir(os.path.join(fa_path, d))],
                key=alpha_sort_key,
            )
            for alpha in alpha_list:
                alpha_dir = os.path.join(fa_path, alpha)
                tasks.append((bh, fa, alpha, alpha_dir))
    return tasks


def server_path(local_path):
    rel = os.path.relpath(local_path, BASE_DIR)
    return os.path.join(SERVER_BASE, rel)


tasks   = find_alpha_dirs(OUT_ROOT, BH_MASSES)
scripts = []
skipped = 0

for bh, fa, alpha, alpha_dir in tasks:
    tag       = f"{bh}_{fa}_{alpha}"
    json_file = os.path.join(RES_DIR, f"{tag}.json")

    if os.path.exists(json_file):
        skipped += 1
        continue

    job_name = f"Scan_{tag}"
    script   = os.path.join(SCRPT_DIR, f"scan_{tag}.sh")

    srv_alpha_dir = server_path(alpha_dir)
    srv_json      = os.path.join(SERVER_BASE, "scan_results", f"{tag}.json")
    srv_log       = os.path.join(SERVER_BASE, "scan_logs", f"{tag}_%j.log")

    lines = [
        "#!/bin/bash\n",
        f"#SBATCH --job-name={job_name}\n",
        f"#SBATCH --partition={PARTITION}\n",
        "#SBATCH --ntasks=1\n",
        "#SBATCH --cpus-per-task=1\n",
        f"#SBATCH --mem={MEM_GB}G\n",
        f"#SBATCH --time={WALLCLOCK}\n",
        f"#SBATCH --output={srv_log}\n",
        f"#SBATCH --chdir={SERVER_BASE}\n",
        "\n",
        f"source {BASHRC}\n",
        "\n",
        f"python3 scan_active_states.py"
        f" --alpha-dir {srv_alpha_dir}"
        f" --output-json {srv_json}\n",
    ]

    with open(script, "w") as fh:
        fh.writelines(lines)
    os.system(f"chmod u+x {script}")
    scripts.append(script)

print(f"Generated {len(scripts)} scan scripts  ({skipped} already done, skipped)")

# Aggregation script (always regenerate)
agg_script = os.path.join(SCRPT_DIR, "scan_aggregate.sh")
srv_res_dir = os.path.join(SERVER_BASE, "scan_results")
srv_agg_log = os.path.join(SERVER_BASE, "scan_logs", "aggregate_%j.log")
agg_lines = [
    "#!/bin/bash\n",
    "#SBATCH --job-name=ScanAggregate\n",
    f"#SBATCH --partition={PARTITION}\n",
    "#SBATCH --ntasks=1\n",
    "#SBATCH --cpus-per-task=1\n",
    f"#SBATCH --mem={AGG_MEM_GB}G\n",
    f"#SBATCH --time={AGG_WALLCLOCK}\n",
    f"#SBATCH --output={srv_agg_log}\n",
    f"#SBATCH --chdir={SERVER_BASE}\n",
    "\n",
    f"source {BASHRC}\n",
    "\n",
    f"python3 scan_active_states.py --aggregate {srv_res_dir}\n",
]
with open(agg_script, "w") as fh:
    fh.writelines(agg_lines)
os.system(f"chmod u+x {agg_script}")

if not scripts and skipped > 0:
    print("All scan results already exist — running aggregation only.")
    if SUBMIT:
        os.system(f"sbatch {agg_script}")
    sys.exit(0)

if not scripts:
    print("No output directories found. Nothing to submit.")
    sys.exit(0)

if SUBMIT:
    job_ids = []
    for sp in scripts:
        result = os.popen(f"sbatch --parsable {sp}").read().strip()
        if result.isdigit():
            job_ids.append(result)
        else:
            print(f"  WARNING: unexpected sbatch output for {os.path.basename(sp)}: {result!r}")

    print(f"Submitted {len(job_ids)}/{len(scripts)} scan jobs.")

    if job_ids:
        dep = "afterok:" + ":".join(job_ids)
        ret = os.system(f"sbatch --dependency={dep} {agg_script}")
        if ret == 0:
            print(f"Submitted aggregation job (runs after all scan jobs finish).")
        else:
            print(f"WARNING: failed to submit aggregation job. Run manually: sbatch {agg_script}")
    else:
        print("No job IDs collected — skipping aggregation submission.")
        print(f"Run aggregation manually when done: sbatch {agg_script}")
else:
    print("Skipped submission (--no-submit).  To submit manually:")
    print(f"  for f in {SCRPT_DIR}/scan_BH_*.sh; do sbatch $f; done")
    print(f"  # then after all finish:")
    print(f"  sbatch {agg_script}")
