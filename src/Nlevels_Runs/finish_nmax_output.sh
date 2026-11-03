#!/bin/bash
# Prune inactive modes on loose .dat, then compress to tarball (post-evolution).
# Usage: finish_nmax_output.sh <outdir> <Nmax>
# Called automatically at the end of each NL_M10 evolution job.
set -euo pipefail

# Prefer env with numpy; compute-node `python3` often has none.
PYTHON_BIN="${PYTHON_BIN:-/lustre/hpc/astro/spieksma/miniforge3/bin/python3}"

OUTDIR="${1:?outdir required}"
NMAX="${2:?Nmax required}"
NL="${NL_BASE_DIR:-$(cd "$(dirname "$0")" && pwd)}"
ARCHIVE="${OUTDIR}/NL_output_Nmax_${NMAX}.tar.gz"
MARKER="${OUTDIR}/NL_output_Nmax_${NMAX}.pruned"

if [[ ! -d "$OUTDIR" ]]; then
  echo "ERROR: outdir missing: $OUTDIR" >&2
  exit 1
fi
# run_Nlevels.jl writes this when the run stopped before tau_max (wall-clock limit):
# keep the loose .dat files so the run can be continued with --resume.
# Fallback for runs started before the marker existed: the job's own log.
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  JOBLOG="$(scontrol show job "$SLURM_JOB_ID" 2>/dev/null | sed -n 's/^ *StdOut=//p')"
  if [[ -n "$JOBLOG" && -f "$JOBLOG" ]] && grep -q "Terminating integration due to wall-clock limit" "$JOBLOG"; then
    echo "stopped at the wall-clock limit (see $JOBLOG); continue with --resume" > "${OUTDIR}/NL_incomplete_Nmax_${NMAX}"
  fi
fi
if [[ -f "${OUTDIR}/NL_incomplete_Nmax_${NMAX}" ]]; then
  echo "[finish] Nmax=$NMAX incomplete ($(cat "${OUTDIR}/NL_incomplete_Nmax_${NMAX}")); not pruning/archiving" >&2
  exit 3
fi
if [[ ! -f "${NL}/prune_states.py" ]]; then
  echo "ERROR: prune_states.py not found in $NL" >&2
  exit 1
fi
if [[ ! -x "${NL}/compress_nmax_output.sh" && ! -f "${NL}/compress_nmax_output.sh" ]]; then
  echo "ERROR: compress_nmax_output.sh not found in $NL" >&2
  exit 1
fi

cd "$OUTDIR"
rm -f "$ARCHIVE" "${ARCHIVE}.tmp" "$MARKER"

echo "[finish] pruning inactive modes in $OUTDIR Nmax=$NMAX"
export PYTHONPATH="${NL}:${PYTHONPATH:-}"
"${PYTHON_BIN}" - <<PY
import sys
sys.path.insert(0, "${NL}")
from prune_states import prune_states_modes
n_in, n_out = prune_states_modes("${OUTDIR}", ${NMAX})
print(f"[finish] prune done: input_modes={n_in} kept_modes={n_out}")
PY

echo "[finish] compressing Nmax=$NMAX"
FORCE=1 bash "${NL}/compress_nmax_output.sh" "${OUTDIR}" "${NMAX}"

if [[ ! -f "$ARCHIVE" || ! -s "$ARCHIVE" ]]; then
  echo "ERROR: compress did not produce $ARCHIVE" >&2
  exit 1
fi

touch "$MARKER"
echo "[finish] pruned + compressed -> $ARCHIVE ($(du -h "$ARCHIVE" | awk '{print $1}'))"
