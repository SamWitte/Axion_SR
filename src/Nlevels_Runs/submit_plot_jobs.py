"""
submit_plot_jobs.py
-------------------
Discover every leaf directory under ./output/ and submit one SLURM job per
directory to run plot_runs.py --directory <dir>.

Usage:
  python submit_plot_jobs.py            # generate scripts AND submit
  python submit_plot_jobs.py --no-submit  # generate scripts only

Scripts  -> ./plot_scripts/
Logs     -> ./plot_logs/
Figures  -> ./figures/  (written by plot_runs.py)

A job is skipped if the corresponding output figure already exists.
"""

import os
import sys

SUBMIT    = "--no-submit" not in sys.argv
BASE_DIR  = os.path.dirname(os.path.abspath(__file__))
OUT_ROOT  = os.path.join(BASE_DIR, "output")
FIG_ROOT  = os.path.join(BASE_DIR, "figures")
SCRPT_DIR = os.path.join(BASE_DIR, "plot_scripts")
LOG_DIR   = os.path.join(BASE_DIR, "plot_logs")

# SLURM resources per job (one directory = one job)
PARTITION  = "short"
MEM_GB     = 20
WALLCLOCK  = "8:00:00"

SERVER_BASE = "/data/phys-axion-ns/phys2564/Axion_SR/src/Nlevels_Runs"
BASHRC      = "/home/phys2564/.bashrc"

os.makedirs(SCRPT_DIR, exist_ok=True)
os.makedirs(LOG_DIR,   exist_ok=True)


def find_leaf_dirs(root):
    leaves = []
    for dirpath, _, filenames in os.walk(root):
        if any(f.endswith(".dat") for f in filenames):
            leaves.append(dirpath)
    return sorted(leaves)


def figure_exists(leaf_dir):
    """Return True if the corresponding figure already exists."""
    rel   = os.path.relpath(leaf_dir, OUT_ROOT)
    parts = rel.split(os.sep)
    if len(parts) >= 3:
        bh, fa, alpha = parts[0], parts[1], parts[2]
    elif len(parts) == 2:
        bh, fa, alpha = parts[0], parts[1], "plot"
    else:
        bh, fa, alpha = parts[0], "unknown_fa", "plot"
    fig_path = os.path.join(FIG_ROOT, bh, fa, alpha + ".jpg")
    return os.path.exists(fig_path)


def server_path(local_path):
    """Translate a local absolute path to the equivalent server path."""
    rel = os.path.relpath(local_path, BASE_DIR)
    return os.path.join(SERVER_BASE, rel)


leaves    = find_leaf_dirs(OUT_ROOT)
scripts   = []
skipped   = 0

for leaf in leaves:
    rel       = os.path.relpath(leaf, OUT_ROOT)
    tag       = rel.replace(os.sep, "_")
    job_name  = f"Plot_{tag}"
    log_file  = os.path.join(LOG_DIR, f"{tag}_%j.log")
    script    = os.path.join(SCRPT_DIR, f"plot_{tag}.sh")

    # Translate paths for the server
    srv_dir     = server_path(leaf)
    srv_log     = os.path.join(SERVER_BASE, "plot_logs", f"{tag}_%j.log")
    srv_script  = os.path.join(SERVER_BASE, "plot_scripts", f"plot_{tag}.sh")

    if figure_exists(leaf):
        skipped += 1
        continue

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
        f"python3 plot_runs.py --directory {srv_dir} --no-zip\n",
    ]

    with open(script, "w") as fh:
        fh.writelines(lines)
    os.system(f"chmod u+x {script}")
    scripts.append(script)

print(f"Generated {len(scripts)} plot scripts  ({skipped} already done, skipped)")

if not scripts:
    print("Nothing to submit.")
    sys.exit(0)

if SUBMIT:
    submitted = 0
    for sp in scripts:
        cmd = f"sbatch {sp}"
        ret = os.system(cmd)
        if ret == 0:
            submitted += 1
        else:
            print(f"  WARNING: sbatch returned {ret} for {os.path.basename(sp)}")
    print(f"Submitted {submitted}/{len(scripts)} plot jobs.")
else:
    print("Skipped submission (--no-submit).  To submit manually:")
    print(f"  for f in {SCRPT_DIR}/plot_*.sh; do sbatch $f; done")
