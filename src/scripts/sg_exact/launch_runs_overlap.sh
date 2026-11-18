#!/bin/bash
# Evolution runs with an overlay table set, started inside a running allocation (takes no job slot):
#   launch_runs_overlap.sh <jobid> <set name> <tables dir> [extra submit_exact_runs.py options, e.g. --alpha 0.15]
# Output: ../runs/<set name>/gw_<model>/M_<mass>/<variant>/alpha_<a>/ ; runs already complete are skipped.
JID=$1; NAME=$2; TAB=$3; shift 3
X=/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext/sg_exact
cd $X || exit 1
/software/python/3.12/bin/python3 submit_exact_runs.py --name $NAME --tables $TAB --no-submit "$@" || exit 1
T=$X/run_tasks_$NAME.sh
[ -s $T ] || { echo "nothing to do"; exit 0; }
cp $T $X/run_tasks_${NAME}_$JID.sh; T=$X/run_tasks_${NAME}_$JID.sh
n=$(wc -l < $T)
L=/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext/logs/runs/$NAME
nohup srun --jobid=$JID --overlap -n1 -c$(( n > 90 ? 90 : n )) bash -c "source ~/.bashrc; cd /data/dust/user/samwitte/Axion_SR_sg/src/Nlevels_Runs
  export JULIA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NL_PROGRESS_SEC=1800 NL_SI_FROM_SG=1 NL_SG_OVERLAY=$X/$TAB
  pids=(); while IFS= read -r line; do [ -z \"\$line\" ] && continue; bash -c \"\$line\" & pids+=(\$!); sleep 1; done < $T
  for p in \"\${pids[@]}\"; do wait \$p; done; echo finished \$(date)" > $L/overlap_$JID.out 2>&1 &
echo "started $n runs of set $NAME (tables $TAB) inside job $JID"
