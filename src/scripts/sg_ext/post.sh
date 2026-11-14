#!/bin/bash
#SBATCH --partition=maxcpu --constraint="7402|75F3|7H12" --nodes=1 --time=0-04:00:00
#SBATCH --job-name=sgx_post --output=/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext/logs/post_%j.out
# After the evolution runs of a set (SGX_SET): summary table and figures.
source ~/.bashrc
X=/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext; SET=${SGX_SET:-full}
cd $X
/software/python/3.12/bin/python3 -I summarize_runs.py --set $SET --csv summary_$SET.csv > summary_$SET.txt 2>&1
./make_figs.sh $SET nonrel; ./make_figs.sh $SET rel
echo "post done $(date)"
