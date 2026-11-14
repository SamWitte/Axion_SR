#!/bin/bash
# Finish the rate campaign (633 channels were incomplete on 2026-11-14: the first rate job was killed by the OOM killer and
# its chunks stayed locked). Leftovers are split in two tiers:
#   p0  cheap: escaping fourth quantum, or every level at n <= 20 and l <= 8   -> one job, then chain (merge, runs, figures)
#   p1  heavy: bound fourth quanta of the high (n,10,10)/(n,16,16)/(n,20,20) towers, minutes and GBs per point
#                                                                              -> NHEAVY jobs with few workers, then chain again
#   ./finish_rates.sh [heavy_workers=20] [heavy_jobs=2]
X=/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext
cd $X
HW=${1:-20}; NH=${2:-2}
/software/python/3.12/bin/python3 -I plan_rates.py --only-missing --write --chunk-sec 1800 | tail -5
rm -rf locks/chunk_* fail/* 2>/dev/null
a=$(sbatch --parsable --export=ALL,SGX_PREFIX=chunk_p0,SGX_WORKERS=48 rates_job.sh)
c1=$(sbatch --parsable --dependency=afterany:$a --export=ALL,SGX_RUNTAG=v1_gap chain.sh)
echo "cheap tier: rates job $a, chain $c1"
dep=""
for i in $(seq $NH); do b=$(sbatch --parsable --export=ALL,SGX_PREFIX=chunk_p1,SGX_WORKERS=$HW rates_job.sh); dep="$dep:$b"; echo "heavy tier: rates job $b ($HW workers)"; done
c2=$(sbatch --parsable --dependency=afterany$dep:$c1 --export=ALL,SGX_RUNTAG=v2_noheavy chain.sh)
echo "final chain $c2"
