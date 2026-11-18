#!/bin/bash
# Single-threaded workers for the Kerr-state Newtonian (+ recoil) and contact amplitudes (scan_gen.jl ... newton_only)
# inside the allocation of a running job:   newton_workers.sh <jobid> <n workers> [outdir=prodN] [jobfile=jobs_newton.txt]
JID=$1; NW=$2; OUT=${3:-prodN}; JF=${4:-jobs_newton.txt}
X=/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext/sg_exact
mkdir -p $X/$OUT/logs
nohup srun --jobid=$JID --overlap -n1 -c$NW bash -c "
  source ~/.bashrc; cd $X
  export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 SGX_WALL_SEC=80000 SLURM_JOB_ID=$JID
  for w in \$(seq 1 $NW); do
    ~/julia-1.8.0/bin/julia /data/dust/user/samwitte/Axion_SR_sg/src/ongoing_tests/sg_relativistic/scan_gen.jl $JF $OUT newton_only > $OUT/logs/n_${JID}_\$w.log 2>&1 &
    sleep 2
  done
  wait" > $X/$OUT/logs/newton_$JID.out 2>&1 &
echo "launched $NW Newtonian workers in job $JID -> $OUT"
