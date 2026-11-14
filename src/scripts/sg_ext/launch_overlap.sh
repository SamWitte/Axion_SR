#!/bin/bash
# Run the evolution runs of a set inside the allocations of running jobs (srun --overlap), so they take no extra job
# slot:   ./launch_overlap.sh <set> <jobid> [<jobid> ...]     (one node task per job id; tasks from submit_sgx_runs.py)
X=/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext
cd $X
SET=$1; shift; JOBS=("$@")
NRUN=$(( ${#JOBS[@]} ))
/software/python/3.12/bin/python3 -I submit_sgx_runs.py --set $SET --no-submit --per-node ${PER_NODE:-30} ${EXTRA_ARGS} | tail -3
NT=$(ls run_tasks/$SET/task_*.sh | wc -l)
k=0
for t in $(seq 0 $((NT - 1))); do
  j=${JOBS[$(( k % NRUN ))]}; k=$((k + 1))
  n=$(wc -l < run_tasks/$SET/task_$(printf %02d $t).sh)
  nohup srun --jobid=$j --overlap -n1 -c$n bash -c "SLURM_ARRAY_TASK_ID=$t bash $X/sgx_runs_$SET.sh" > logs/runs/$SET/overlap_task${t}_job$j.out 2>&1 &
  echo "task $t ($n runs) -> job $j"
done
