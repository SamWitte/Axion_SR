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
