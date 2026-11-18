#!/bin/bash
#SBATCH --partition=maxcpu
#SBATCH --constraint="7402|75F3|7H12"
#SBATCH --nodes=1
#SBATCH --time=14:00:00
#SBATCH --job-name=sgx_exact
#SBATCH --output=/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext/sg_exact/logs/exact_%j.out
# As exact_job.sh, but robust against the node's memory watchdog: each worker passes over the job files SGX_PASSES
# times (a worker that is terminated on a heavy point starts again with the next file), and the allocation is kept
# until the file HOLD_<jobid>.stop exists or SGX_HOLD_H hours have passed, so that evolution runs can be started in it
# with srun --overlap without depending on the exact workers.   End it with: touch HOLD_<jobid>.stop  (or scancel).
#   SGX_JOBS (default: jobs_todoA.txt jobs_todoB.txt), SGX_WORKERS (3), SGX_THREADS (cores / workers), SGX_OUT (prod2)
source ~/.bashrc
X=/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext/sg_exact
cd $X
NW=${SGX_WORKERS:-3}
NCPU=$(nproc)
export JULIA_NUM_THREADS=${SGX_THREADS:-$(( NCPU / NW ))}
export OPENBLAS_NUM_THREADS=$JULIA_NUM_THREADS
export SGX_WALL_SEC=${SGX_WALL_SEC:-40000}
OUT=${SGX_OUT:-prod2}
JOBS=${SGX_JOBS:-"jobs_todoA.txt jobs_todoB.txt"}
J=~/julia-1.8.0/bin/julia
S=/data/dust/user/samwitte/Axion_SR_sg/src/ongoing_tests/sg_relativistic/scan_gen.jl
mkdir -p $OUT/logs
echo "node $(hostname) cpus $NCPU workers $NW threads $JULIA_NUM_THREADS jobs $JOBS out $OUT start $(date)"
for w in $(seq 1 $NW); do
  ( for pass in $(seq 1 ${SGX_PASSES:-10}); do
      [ -e $X/HOLD_${SLURM_JOB_ID}.stop ] && break
      for jf in $JOBS; do $J $S $jf $OUT; done
      sleep 240
    done ) > $OUT/logs/w_${SLURM_JOB_ID}_$w.log 2>&1 &
  sleep 3
done
t0=$(date +%s)
while [ ! -e $X/HOLD_${SLURM_JOB_ID}.stop ] && [ $(( $(date +%s) - t0 )) -lt $(( ${SGX_HOLD_H:-12} * 3600 )) ]; do sleep 60; done
echo "hold released $(date)"
