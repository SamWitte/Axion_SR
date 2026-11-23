#!/bin/bash
#SBATCH --partition=allcpu --constraint="75F3|7542|7402|7H12" --nodes=1 --time=12:00:00 --requeue
#SBATCH --job-name=spintab --output=/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext/spin_tab/logs/job_%A_%a.out
# One node of the spin-table campaign. Environment: MODE, LIST=<name in tasks/>, NPER=<workers on this node>,
# NARR=<number of array tasks sharing the list>. Restart-safe (preemptible partition): finished rows are skipped.
#   MODE=levels  level energies and n <= 8 radial functions at the 9 spins (spin_levels.jl)
#   MODE=fast    spin_fast.jl, all spins, levels with n <= 8 relativistic            -> fast_out/fast_*
#   MODE=fast7   spin_fast.jl, spins 0.95 0.99 0.998, the n = 8 levels hydrogenic    -> fast_out/fast7_*
#   MODE=fix     spin_fast.jl, spin 0.998 only, rows with a broken cached radial function (wf_bad.txt) -> fast_out/fix998_*
#   MODE=full    spin_points.jl (production gf_radial; validation)                   -> pts/
source ~/.bashrc
D=/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext/spin_tab
X=/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext
cd $D
export JULIA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
JL=~/julia-1.8.0/bin/julia
ALPHAS=$(cat $X/sgx_alphas.txt)
SPINS=0,0.3,0.5,0.7,0.8,0.9,0.95,0.99,0.998
NPER=${NPER:-64}; NARR=${NARR:-1}; NTOT=$(( NPER * NARR )); OFF=$(( ${SLURM_ARRAY_TASK_ID:-0} * NPER ))
echo "job ${SLURM_JOB_ID:-local} task ${SLURM_ARRAY_TASK_ID:-0} on $(hostname): MODE=$MODE LIST=$LIST workers $NPER of $NTOT, start $(date)"
for i in $(seq $NPER); do
  k=$(( OFF + i ))
  case "$MODE" in
    levels) $JL spin_levels.jl tasks/levels.txt $ALPHAS $SPINS cache levels_out $k $NTOT both >> logs/levels_$k.log 2>&1 & ;;
    fast)   $JL spin_fast.jl tasks/$LIST.txt $ALPHAS $SPINS levels_out cache fast_out $k $NTOT hmax=0.25 >> logs/fast_${LIST}_$k.log 2>&1 & ;;
    fast7)  $JL spin_fast.jl tasks/$LIST.txt $ALPHAS 0.95,0.99,0.998 levels_out cache fast_out $k $NTOT rel_nmax=7 tag=fast7 hmax=0.25 >> logs/fast7_${LIST}_$k.log 2>&1 & ;;
    fix)    $JL spin_fast.jl tasks/$LIST.txt $ALPHAS 0.998 levels_out cache fast_out $k $NTOT hmax=0.25 tag=fix998 wfsub=wf_bad.txt only_sub=1 >> logs/fix998s_${LIST}_$k.log 2>&1 & ;;
    full)   $JL spin_points.jl tasks/$LIST.txt pts cache $k $NTOT >> logs/${LIST}_$k.log 2>&1 & ;;
    *) echo "unknown MODE '$MODE'"; exit 2 ;;
  esac
  sleep 0.2
done
wait
echo "end $(date)"
