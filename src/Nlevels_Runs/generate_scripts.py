"""
generate_scripts.py
-------------------
Generates and (optionally) submits queue scripts for the Nlevels superradiance
parameter sweep.

Grid:
  MassBH  : [1e8, 1e4, 10]  solar masses
  SpinBH  : 0.99             (fixed)
  f_a     : [1e18, 1e16, 1e14, 1e12]  eV
  alpha   : [0.05, 0.1, 0.25, 0.4, 0.6, 0.9, 1.2, 1.5]
  Nmax    : [3, 4, 5, 6, 7, 8, 15]    (run sequentially in each script)

Grouping: one shell script per (MassBH, f_a, alpha) -- 96 scripts total.
Each script runs all 7 Nmax values in sequence for that parameter combination.

Output layout:
  <BASE_DIR>/output/
    BH_1e8/
      fa_1e18/
        alpha_0.05/
          Time_FullRel_fa_...Nmax_3.dat
          Spin_FullRel_fa_...Nmax_3.dat
          ...
    BH_1e4/
      ...
    BH_10/
      ...

Usage:
  python generate_scripts.py           # create scripts AND submit
  python generate_scripts.py --no-submit  # create scripts only
"""

import os
import re
import sys
import itertools

# ---------------------------------------------------------------------------
# Control
# ---------------------------------------------------------------------------
SUBMIT = "--no-submit" not in sys.argv   # pass --no-submit to skip addqueue

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR   = os.path.dirname(os.path.abspath(__file__))   # .../Nlevels_Runs
SCRIPT_DIR = os.path.join(BASE_DIR, "scripts")
OUTPUT_DIR = os.path.join(BASE_DIR, "output")

os.makedirs(SCRIPT_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Parameter grids
# ---------------------------------------------------------------------------
# MassBH_vals = [1e8, 1e4, 10.0]
MassBH_vals = [1e4]
SpinBH      = 0.95
fa_vals     = [1e18, 1e16, 1e14, 1e12]
# alpha_vals  = [0.05, 0.2, 0.4, 0.6]
alpha_vals  = [0.05, 0.2, 0.4, 0.6, 0.8, 1.0, 1.3, 1.6]
Nmax_vals        = [3, 4, 5, 6, 7]   # run sequentially in one script
Nmax_solo_vals   = [8, 15]            # each gets its own dedicated script

tau_max_map = {
    1e8:  1e9,
    1e4:  1e9,
    10.0: 5e7,
}

# Queue resource request
MEM_GB    = 5
WALLCLOCK = "1 week"

# ---------------------------------------------------------------------------
# Helper formatters for directory names
# ---------------------------------------------------------------------------
def bh_tag(mass):
    """Return a clean directory-friendly tag for a BH mass."""
    if mass >= 1e7:
        exp = int(round(__import__("math").log10(mass)))
        return f"BH_1e{exp}"
    elif mass >= 1e3:
        exp = int(round(__import__("math").log10(mass)))
        return f"BH_1e{exp}"
    else:
        return f"BH_{int(mass)}"

def fa_tag(fa):
    """Return a clean directory-friendly tag for f_a."""
    import math
    exp = int(round(math.log10(fa)))
    return f"fa_1e{exp}"

def alpha_tag(a):
    """Return a clean directory-friendly tag for alpha."""
    return f"alpha_{a:g}"

def fmt(x):
    """Format a float in clean scientific notation (no leading zeros in exponent)."""
    s = f"{x:.6g}"
    s = re.sub(r'e\+0*(\d+)', r'e\1', s)
    s = re.sub(r'e-0*(\d+)',  r'e-\1', s)
    return s

# ---------------------------------------------------------------------------
# Generate scripts
# ---------------------------------------------------------------------------
script_paths = []   # (script_path, bh_label) pairs for submission

job_count = 0
for MassBH, fa, alpha in itertools.product(MassBH_vals, fa_vals, alpha_vals):

    tau_max = tau_max_map[MassBH]

    # Directory that Julia will write output files into
    outdir = os.path.join(
        OUTPUT_DIR,
        bh_tag(MassBH),
        fa_tag(fa),
        alpha_tag(alpha),
    )

    base_name = (
        f"Script_NL_{bh_tag(MassBH)}_"
        f"{fa_tag(fa)}_"
        f"{alpha_tag(alpha)}"
    )

    def make_script(name, nmax_list):
        lines = []
        lines.append("#!/bin/bash\n")
        lines.append(f"source /mnt/users/switte/.bashrc\n")
        lines.append(f"cd {BASE_DIR}\n")
        lines.append("\n")
        lines.append(
            f"# MassBH={fmt(MassBH)} M_sun | f_a={fmt(fa)} eV | alpha={alpha:g}\n"
        )
        lines.append("\n")
        for Nmax in nmax_list:
            lines.append(f"# --- Nmax = {Nmax} ---\n")
            lines.append(
                f"srun --exclusive julia run_Nlevels.jl"
                f" --MassBH {fmt(MassBH)}"
                f" --SpinBH {SpinBH}"
                f" --f_a {fmt(fa)}"
                f" --alpha {alpha:g}"
                f" --Nmax {Nmax}"
                f" --tau_max {fmt(tau_max)}"
                f" --outdir {outdir}\n"
            )
            lines.append("wait\n")
            lines.append("\n")
        path = os.path.join(SCRIPT_DIR, name + ".sh")
        with open(path, "w") as fh:
            fh.writelines(lines)
        os.system(f"chmod u+x {path}")
        return path

    # Small Nmax values run sequentially in one script
    script_paths.append(make_script(base_name, Nmax_vals))
    job_count += 1

    # Nmax=8 and Nmax=15 each get their own dedicated script
    for Nmax in Nmax_solo_vals:
        script_paths.append(make_script(f"{base_name}_Nmax_{Nmax}", [Nmax]))
        job_count += 1

print(f"Created {job_count} scripts in {SCRIPT_DIR}")

# ---------------------------------------------------------------------------
# Submit
# ---------------------------------------------------------------------------
if SUBMIT:
    # addqueue expects to be called from the directory containing the script
    os.chdir(SCRIPT_DIR)
    submitted = 0
    for sp in script_paths:
        sname = os.path.basename(sp)
        cmd = f'addqueue -s -n 1 -c "{WALLCLOCK}" -m {MEM_GB} ./{sname}'
        ret = os.system(cmd)
        if ret == 0:
            submitted += 1
        else:
            print(f"  WARNING: addqueue returned {ret} for {sname}")
    print(f"Submitted {submitted}/{job_count} jobs.")
else:
    print("Skipped submission (--no-submit).  Run with:")
    print(f"  cd {SCRIPT_DIR}")
    print(f'  addqueue -s -n 1 -c "{WALLCLOCK}" -m {MEM_GB} ./Script_NL_*.sh')
