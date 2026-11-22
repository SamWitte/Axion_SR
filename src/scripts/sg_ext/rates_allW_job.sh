#!/bin/bash
#SBATCH --partition=allcpu --constraint="75F3|7542|7402|7H12" --nodes=1 --time=1-00:00:00 --requeue
#SBATCH --job-name=sgx_allW --output=/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext/allW/logs/rates_%j.out
# Recomputation of every rate point in the "to infinity" regime with all outgoing multipoles AND the corrected
# outgoing condition / outer radius (gf_radial inf_all_l, inf_wkb).  Supersedes allL.
# Same scheme as ../rates_job.sh with its own directories: a node works through tasks/chunk_*.txt, a chunk is claimed
# with an atomic mkdir in locks/ and marked in done/; points go to rates_out/<channel>.pts as soon as they are computed
# and points already there are skipped, so the job can be stopped and resubmitted. Any number of these jobs can run.
source ~/.bashrc
X=/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext; A=$X/allW
cd $A/src   # snapshot of the code taken when the campaign was set up (src/SNAPSHOT_INFO.txt), so that later edits of the worktree cannot change it
export JULIA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
JL=~/julia-1.8.0/bin/julia
ALPHAS=$(cat $X/sgx_alphas.txt)
ME=${SLURM_JOB_ID:-local$$}
NCPU=${SLURM_CPUS_ON_NODE:-96}; NW=${SGX_WORKERS:-$(( NCPU * 3 / 4 ))}
echo "job $ME on $(hostname): $NCPU threads, $NW workers, start $(date)"
unlock() { rm -f "$1/job" "$1/skip"; rmdir "$1" 2>/dev/null; }
# allcpu is preemptible (requeue): a restarted job first drops the locks of its own earlier incarnation
for l in $A/locks/chunk_*; do [ -d "$l" ] && [ "$(cat $l/job 2>/dev/null)" = "$ME" ] && unlock "$l"; done
# stale locks (owner job gone), older than 3 minutes
for l in $(find $A/locks -mindepth 1 -maxdepth 1 -type d -mmin +3 -name 'chunk_*'); do
  j=$(cat $l/job 2>/dev/null)
  if [ -z "$j" ] || ! squeue -h -j "$j" 2>/dev/null | grep -q .; then unlock "$l"; fi
done
worker() {
  while :; do
    got=""
    for f in $A/tasks/chunk_*.txt; do
      b=$(basename $f .txt)
      [ -e $A/done/$b ] && continue
      [ -e $A/locks/$b ] && continue
      if mkdir $A/locks/$b 2>/dev/null; then
        echo $ME > $A/locks/$b/job
        if [ -e $A/done/$b ]; then unlock $A/locks/$b; continue; fi
        if [ -e $A/fail/$b ] && [ "$(wc -l < $A/fail/$b)" -ge 2 ]; then touch $A/locks/$b/skip; continue; fi   # failed twice: stays locked for this job
        $JL Compute_sg_points.jl --tasks $f --out $A/rates_out --alpha_list $ALPHAS --cache $X/cache --inf_wkb true >> $A/logs/chunks/$b.log 2>&1
        if tail -n 1 $A/logs/chunks/$b.log | grep -q TASKS_DONE; then touch $A/done/$b; else echo "$ME $(date)" >> $A/fail/$b; fi
        unlock $A/locks/$b
        got=1; break
      fi
    done
    [ -z "$got" ] && break
  done
}
for i in $(seq $NW); do worker & sleep 0.3; done
wait
echo "job $ME: no chunks left, $(ls $A/done | wc -l) done of $(ls $A/tasks | wc -l), end $(date)"
