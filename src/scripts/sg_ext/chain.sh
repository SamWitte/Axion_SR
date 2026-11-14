#!/bin/bash
#SBATCH --partition=maxcpu --constraint="7402|75F3|7H12" --nodes=1 --time=0-06:00:00
#SBATCH --job-name=sgx_chain --output=/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext/logs/chain_%j.out
# Runs after the rate jobs: report what is missing, install the tables, submit the evolution runs (full set).
# If more than 2% of the channels are incomplete it stops instead (CHAIN_STOPPED.txt) so the gap can be looked at.
source ~/.bashrc
X=/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext
PY=/software/python/3.12/bin/python3
cd $X
echo "chain start $(date)"
$PY -I plan_rates.py --only-missing | tee logs/chain_missing.txt
miss=$(grep -oE "^tasks: [0-9]+" logs/chain_missing.txt | grep -oE "[0-9]+")
tot=$(wc -l < sgx_channels.txt)
echo "incomplete channels: $miss of $tot"
if [ -z "$miss" ] || [ $(( miss * 50 )) -gt "$tot" ]; then
  echo "more than 2% of the channels incomplete: not merging; rerun 'sbatch rates_job.sh' (it resumes) and then this script" | tee CHAIN_STOPPED.txt
  exit 1
fi
rm -f CHAIN_STOPPED.txt
$PY -I merge_rates.py --apply | tee logs/chain_merge.txt
# quota: the runs write ~15-20 GB in total, and only at the end of each run; do not start them on a nearly full fileset
used=$(mmlsquota -j user.it.samwitte dust 2>/dev/null | awk '$1=="dust"{print $3}')
echo "dust fileset: ${used:-?} KB used of 2147483648"
if [ -n "$used" ] && [ "$used" -gt 2084000000 ]; then
  echo "less than ~60 GB left on the dust fileset: runs not submitted. Free space, then: python3 submit_sgx_runs.py --set full" | tee CHAIN_STOPPED.txt
  exit 2
fi
# SGX_RUNTAG=<tag>: the existing runs of the full set were made with the previous tables; set them aside so the runs are redone
if [ -n "$SGX_RUNTAG" ] && [ -d runs/full ]; then
  mv runs/full runs/full_$SGX_RUNTAG
  [ -d figs/full ] && mv figs/full figs/full_$SGX_RUNTAG
  for f in summary_full.txt summary_full.csv; do [ -e $f ] && mv $f ${f/full/full_$SGX_RUNTAG}; done
  [ -d logs/runs/full ] && mv logs/runs/full logs/runs/full_$SGX_RUNTAG
fi
$PY -I submit_sgx_runs.py --set full | tee logs/chain_submit.txt
echo "chain end $(date)"
