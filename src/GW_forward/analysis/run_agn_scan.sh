#!/bin/bash
# AGN companion-resonance pipeline over a set of boson masses, one run directory
# per alpha (alpha = G M_ref mu, M_ref = 1e6 M_sun by default):
#
#   ./run_agn_scan.sh [options] alpha1 [alpha2 ...]
#
# e.g.  ./run_agn_scan.sh -j 48 -N 100 -n 20 -m 4 0.05 0.1 0.2 0.3 0.5
#
# For every alpha:
#   1. export   N accretion histories into OUT/alpha<alpha> (forward seeding with
#               the LISA-band cut, agn_julia_pipeline.py PARAMS). The same seed
#               gives the same histories for every alpha.
#   2. evolve   one Julia job per history (evolve_agn_history.jl); the histories
#               of all alphas share one pool of JOBS single-threaded processes.
#   3. analyze  n_inj companions per history -> OUT/alpha<alpha>/first_resonance_hist.png,
#               first_resonances.csv, crossings.csv, ...
# then agn_scan_summary.py combines the run directories into OUT/scan_summary.png
# and OUT/scan_{kinds,transitions}.csv.
#
# Options:
#   -j JOBS   parallel Julia jobs (default: all cores)
#   -N NHIST  histories per alpha (default 100)
#   -n NINJ   companions per history (default 20)
#   -m NMAX   Nmax of the cloud evolution (default 4)
#   -f FA     axion decay constant [GeV] (default 1e18)
#   -s SEED   history seed (default 7)
#   -r MREF   reference mass for alpha [M_sun] (default 1e6)
#   -o OUT    output root (default scan_runs)
#   -q        also analyze with quadrupole tides only (l* = 2; files tagged _quad)
#   -A        analysis only: skip export and evolution
#
# Re-running the same command resumes: existing histories and finished
# evolutions (evol_XXXX.csv) are kept, so a killed scan can simply be restarted.
# Environment: PYTHON (default python3), JULIA (default julia), NL_MAX_WALL_SEC
# (wall-clock cap per solver call in s, default 3600, see evolve_agn_history.jl).
set -euo pipefail
# one core per job: no multithreaded BLAS / Julia threads
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1} JULIA_NUM_THREADS=${JULIA_NUM_THREADS:-1}
export PYTHON=${PYTHON:-python3} JULIA=${JULIA:-julia}

ncores() { nproc 2>/dev/null || sysctl -n hw.ncpu; }
JOBS=$(ncores); NHIST=100; NINJ=20; NMAX=4; FA=1e18; SEED=7; MREF=1e6; OUT=scan_runs; QUAD=0; ANALYZE_ONLY=0
while getopts "j:N:n:m:f:s:r:o:qAh" opt; do
    case $opt in
        j) JOBS=$OPTARG ;; N) NHIST=$OPTARG ;; n) NINJ=$OPTARG ;; m) NMAX=$OPTARG ;;
        f) FA=$OPTARG ;; s) SEED=$OPTARG ;; r) MREF=$OPTARG ;; o) OUT=$OPTARG ;;
        q) QUAD=1 ;; A) ANALYZE_ONLY=1 ;;
        *) sed -n '2,35p' "$0"; exit 1 ;;
    esac
done
shift $((OPTIND - 1))
[ $# -ge 1 ] || { echo "usage: $0 [options] alpha1 [alpha2 ...]  (-h for help)"; exit 1; }
ALPHAS=("$@")
cd "$(dirname "$0")"
mkdir -p "$OUT"
export FA NMAX
log() { echo "[$(date '+%F %T')] $*"; }

rundir() { echo "$OUT/alpha$1"; }
n_files() { ls "$1"/$2 2>/dev/null | wc -l | tr -d ' '; }

if [ "$ANALYZE_ONLY" = 0 ]; then
    # 1. export (skipped if the run directory already holds N histories for this mu)
    for a in "${ALPHAS[@]}"; do
        d=$(rundir "$a")
        if [ -f "$d/params.json" ]; then
            "$PYTHON" - "$d/params.json" "$a" "$MREF" <<'EOF'
import json, sys
import agn_julia_pipeline as P
mu_run = json.load(open(sys.argv[1]))["mu_eV"]
mu_req = P.mu_from_alpha(float(sys.argv[2]), float(sys.argv[3]))
if abs(mu_run / mu_req - 1) > 1e-3:
    sys.exit(f"{sys.argv[1]}: mu = {mu_run:.6g} eV, but alpha = {sys.argv[2]} at M_ref = {sys.argv[3]} "
             f"needs {mu_req:.6g} eV. Use another -o or remove the directory.")
EOF
        fi
        if [ "$(n_files "$d" 'hist_*.json')" -lt "$NHIST" ]; then
            log "export alpha=$a -> $d"
            "$PYTHON" agn_julia_pipeline.py export --rundir "$d" --N "$NHIST" --seed "$SEED" \
                --alpha "$a" --M_ref "$MREF" > "$d.export.log" 2>&1 || { cat "$d.export.log"; exit 1; }
            mv "$d.export.log" "$d/export.log"
            grep acceptance "$d/export.log" || true
        fi
    done

    # 2. evolve every history without an evol_XXXX.csv, all alphas in one pool
    jobs_file="$OUT/jobs.txt"
    : > "$jobs_file"
    for a in "${ALPHAS[@]}"; do
        d=$(rundir "$a")
        mu=$("$PYTHON" -c "import json; print(repr(json.load(open('$d/params.json'))['mu_eV']))")
        for h in "$d"/hist_*.txt; do
            id=$(basename "$h" .txt); id=${id#hist_}
            [ -f "$d/evol_$id.csv" ] || echo "$d $((10#$id)) $mu" >> "$jobs_file"
        done
    done
    njobs=$(wc -l < "$jobs_file" | tr -d ' ')
    if [ "$njobs" -gt 0 ]; then
        log "precompiling the Julia code once before starting $njobs jobs on $JOBS cores"
        "$JULIA" -e 'using Suppressor; @suppress include(joinpath("..", "..", "super_rad.jl"))' > /dev/null 2>&1 || true
        log "evolving $njobs histories"
        # shellcheck disable=SC2016
        xargs -P "$JOBS" -n 3 bash -c '
            "$JULIA" evolve_agn_history.jl --rundir "$0" --mu "$2" --f_a "$FA" --Nmax "$NMAX" --ids "$1" \
                > "$0/evolve_$1.log" 2>&1 || echo "  julia failed: $0 id $1 (see $0/evolve_$1.log)"
            echo "  done: $0 id $1"' < "$jobs_file"
    fi
    "$PYTHON" - "${ALPHAS[@]/#/$OUT/alpha}" <<'EOF'
import glob, json, sys
for d in sys.argv[1:]:
    st = [json.load(open(f))["status"] for f in sorted(glob.glob(f"{d}/evol_*_meta.json"))]
    bad = [s for s in st if s != "ok"]
    print(f"  {d}: {len(st) - len(bad)} ok, {len(bad)} not ok" + (f"  e.g. {bad[0][:100]}" if bad else ""))
EOF
fi

# 3. analyze every alpha, then combine
for a in "${ALPHAS[@]}"; do
    d=$(rundir "$a")
    [ "$(n_files "$d" 'evol_*.csv')" -gt 0 ] || { log "no evolved histories in $d, skipped"; continue; }
    log "analyze $d"
    "$PYTHON" agn_julia_pipeline.py analyze --rundir "$d" --n_inj "$NINJ" > "$d/analyze.log" 2>&1 \
        || { tail -5 "$d/analyze.log"; continue; }
    sed -n '/companions in/,/no_SR_state/p' "$d/analyze.log"
    if [ "$QUAD" = 1 ]; then
        "$PYTHON" agn_julia_pipeline.py analyze --rundir "$d" --n_inj "$NINJ" --l_stars 2 --tag quad \
            > "$d/analyze_quad.log" 2>&1 || tail -5 "$d/analyze_quad.log"
    fi
done
dirs=(); for a in "${ALPHAS[@]}"; do dirs+=("$(rundir "$a")"); done
"$PYTHON" agn_scan_summary.py --out "$OUT" --M_ref "$MREF" "${dirs[@]}"
[ "$QUAD" = 1 ] && "$PYTHON" agn_scan_summary.py --out "$OUT" --M_ref "$MREF" --tag quad "${dirs[@]}"
log "finished; see $OUT/scan_summary.png and $OUT/alpha*/first_resonance_hist.png"
