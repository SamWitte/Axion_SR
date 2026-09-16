#!/bin/bash
# Prune + plot one evolution leaf (directory with NL_output_Nmax_*.tar.gz).
#
# Usage:
#   bash postprocess_leaf.sh output/BH_10/fa_1e16/alpha_0.6
#   RUN_TAG=run_nodrag bash postprocess_leaf.sh output/run_nodrag/BH_10/fa_1e16/alpha_0.6
#   bash postprocess_leaf.sh --prune-only output/BH_10/fa_1e16/alpha_0.6
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
export NL_BASE_DIR="${NL_BASE_DIR:-$HERE}"

PRUNE_ONLY=0
LEAF=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --prune-only) PRUNE_ONLY=1; shift ;;
    -h|--help)
      sed -n '2,8p' "$0"
      exit 0
      ;;
    *)
      LEAF="$1"
      shift
      ;;
  esac
done

if [[ -z "$LEAF" ]]; then
  echo "Usage: bash postprocess_leaf.sh [--prune-only] <leaf_dir>" >&2
  exit 1
fi

# Resolve relative paths against NL_BASE_DIR when needed
if [[ "$LEAF" != /* ]]; then
  if [[ -d "$LEAF" ]]; then
    LEAF="$(cd "$LEAF" && pwd)"
  elif [[ -d "${NL_BASE_DIR}/${LEAF}" ]]; then
    LEAF="$(cd "${NL_BASE_DIR}/${LEAF}" && pwd)"
  else
    echo "Leaf not found: $LEAF" >&2
    exit 1
  fi
fi

cd "$HERE"
PY=(python3 -u plot_output_pipeline.py --nl-base "$NL_BASE_DIR" --directory "$LEAF")
if [[ -n "${RUN_TAG:-}" ]]; then
  PY+=(--run-tag "$RUN_TAG")
fi
if (( PRUNE_ONLY )); then
  PY+=(--prune-only)
fi

echo "[postprocess] leaf=$LEAF"
"${PY[@]}"
echo "[postprocess] done"
