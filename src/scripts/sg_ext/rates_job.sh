#!/bin/bash
#SBATCH --partition=maxcpu --constraint="7402|75F3|7H12" --nodes=1 --time=7-00:00:00
#SBATCH --job-name=sgx_rates --output=/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext/logs/rates_%j.out
# Extended SG rate campaign: one node works through tasks/chunk_*.txt (file-name order = priority) with as many
# single-threaded Julia processes as 3/4 of its hardware threads. Any number of these jobs can run at once: a chunk is
# claimed with an atomic mkdir in locks/ and marked in done/ when its task list finished; every rate point is appended
# to rates_out/<channel>.pts as soon as it is computed, and points already there are skipped, so a killed job loses
# nothing. First, the level energies and relativistic wavefunctions are precomputed into cache/ (once; other jobs wait).
source ~/.bashrc
X=/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext
cd /data/dust/user/samwitte/Axion_SR_sg/src
export JULIA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
JL=~/julia-1.8.0/bin/julia
ALPHAS=$(cat $X/sgx_alphas.txt)
ME=${SLURM_JOB_ID:-local$$}
NCPU=${SLURM_CPUS_ON_NODE:-96}; NW=${SGX_WORKERS:-$(( NCPU * 3 / 4 ))}
echo "job $ME on $(hostname): $NCPU threads, $NW workers, start $(date)"

# ---- precompute (one job does it, the others wait for the marker)
if [ ! -e $X/cache/PRECOMPUTE_DONE ]; then
  if mkdir $X/locks/_precompute 2>/dev/null; then
    echo $ME > $X/locks/_precompute/job
    for i in $(seq $NW); do
      $JL Compute_sg_points.jl --precompute $X/levels_sgx.txt --alpha_list $ALPHAS --cache $X/cache --part $i/$NW > $X/logs/pre/pre_${ME}_$i.log 2>&1 &
    done
    wait
    n=$(grep -l PRECOMPUTE_DONE $X/logs/pre/pre_${ME}_*.log | wc -l)
    echo "precompute: $n of $NW parts finished, $(date)"
    [ "$n" -eq "$NW" ] && touch $X/cache/PRECOMPUTE_DONE
    rm -rf $X/locks/_precompute
  else
    echo "waiting for the precompute of job $(cat $X/locks/_precompute/job 2>/dev/null)"
    while [ ! -e $X/cache/PRECOMPUTE_DONE ] && [ -d $X/locks/_precompute ]; do sleep 30; done
  fi
fi

# ---- stale locks (owner job gone), older than 3 minutes
for l in $(find $X/locks -mindepth 1 -maxdepth 1 -type d -mmin +3 -name 'chunk_*'); do
  j=$(cat $l/job 2>/dev/null)
  if [ -z "$j" ] || ! squeue -h -j "$j" 2>/dev/null | grep -q .; then rm -rf "$l"; fi
done

worker() {
  while :; do
    got=""
    for f in $X/tasks/${SGX_PREFIX:-chunk}_*.txt; do      # SGX_PREFIX=chunk_p0 restricts a job to one priority tier
      b=$(basename $f .txt)
      [ -e $X/done/$b ] && continue
      [ -e $X/locks/$b ] && continue
      if mkdir $X/locks/$b 2>/dev/null; then
        echo $ME > $X/locks/$b/job
        if [ -e $X/done/$b ] || { [ -e $X/fail/$b ] && [ "$(wc -l < $X/fail/$b)" -ge 2 ]; }; then
          [ -e $X/done/$b ] || touch $X/locks/$b/skip      # failed twice: leave it locked for this job so nobody loops on it
          [ -e $X/done/$b ] && rm -rf $X/locks/$b
          continue
        fi
        $JL Compute_sg_points.jl --tasks $f --out $X/rates_out --alpha_list $ALPHAS --cache $X/cache >> $X/logs/chunks/$b.log 2>&1
        if tail -n 1 $X/logs/chunks/$b.log | grep -q TASKS_DONE; then touch $X/done/$b; else echo "$ME $(date)" >> $X/fail/$b; fi
        rm -rf $X/locks/$b
        got=1; break
      fi
    done
    [ -z "$got" ] && break
  done
}
for i in $(seq $NW); do worker & sleep 0.3; done
wait
echo "job $ME (${SGX_PREFIX:-chunk}): no chunks left, $(ls $X/done | wc -l) done of $(ls $X/tasks | wc -l), end $(date)"
