#!/bin/bash
#SBATCH --partition=maxcpu
#SBATCH --constraint="7402|75F3|7H12"
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --job-name=sgx_exact
#SBATCH --output=/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext/sg_exact/logs/exact_%j.out
# One node of workers for the exact self-gravity exchange (scan_gen.jl).  Workers claim jobs from the job files
# (in the order given) with lock directories, so any number of these jobs can run on the same lists.
#   SGX_JOBS    job files, space separated (default: jobs_A.txt jobs_B.txt jobs_C.txt)
#   SGX_WORKERS workers on this node (default 12),  SGX_THREADS threads per worker (default cores / workers)
#   SGX_OUT     output directory (default prod)
source ~/.bashrc
cd /data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext/sg_exact
NW=${SGX_WORKERS:-12}
NCPU=$(nproc)
export JULIA_NUM_THREADS=${SGX_THREADS:-$(( NCPU / NW ))}
export OPENBLAS_NUM_THREADS=$JULIA_NUM_THREADS
export SGX_WALL_SEC=${SGX_WALL_SEC:-82000}
OUT=${SGX_OUT:-prod}
JOBS=${SGX_JOBS:-"jobs_A.txt jobs_B.txt jobs_C.txt"}
J=~/julia-1.8.0/bin/julia
S=/data/dust/user/samwitte/Axion_SR_sg/src/ongoing_tests/sg_relativistic/scan_gen.jl
mkdir -p $OUT/logs
echo "node $(hostname) cpus $NCPU workers $NW threads $JULIA_NUM_THREADS jobs $JOBS out $OUT"
for w in $(seq 1 $NW); do
  ( for jf in $JOBS; do $J $S $jf $OUT; done ) > $OUT/logs/w_${SLURM_JOB_ID}_$w.log 2>&1 &
  sleep 3
done
wait
