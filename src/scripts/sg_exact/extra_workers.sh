#!/bin/bash
# Extra exact-exchange workers inside the allocation of a running job (takes no job slot):
#   extra_workers.sh <jobid> <n workers> <threads> <outdir> <job files...>
# One call = one job step: if the kernel kills one worker for memory, Slurm ends the whole step, so start big workers
# one per call.  PASSES=<n> (environment): pass over the job files n times (a pass ends at the end of each list, so
# points released later need another pass).
JID=$1; NW=$2; NT=$3; OUT=$4; shift 4; JOBS="$*"; TAG=$(date +%H%M%S%N | cut -c1-10)
X=/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext/sg_exact
nohup srun --jobid=$JID --overlap -n1 -c$(( NW * NT )) bash -c "
  source ~/.bashrc; cd $X
  export JULIA_NUM_THREADS=$NT OPENBLAS_NUM_THREADS=$NT SGX_WALL_SEC=80000 SLURM_JOB_ID=$JID
  for w in \$(seq 1 $NW); do
    ( for pass in \$(seq 1 ${PASSES:-1}); do for jf in $JOBS; do ~/julia-1.8.0/bin/julia /data/dust/user/samwitte/Axion_SR_sg/src/ongoing_tests/sg_relativistic/scan_gen.jl \$jf $OUT; done; [ ${PASSES:-1} -gt 1 ] && sleep 120; done ) > $OUT/logs/x_${JID}_${TAG}_\$w.log 2>&1 &
    sleep 3
  done
  wait" > $OUT/logs/extra_${JID}_$TAG.out 2>&1 &
echo "launched $NW extra workers ($NT threads) in job $JID"
