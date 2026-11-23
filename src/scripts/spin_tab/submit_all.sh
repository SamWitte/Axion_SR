#!/bin/bash
# Submits the spin-table campaign. Run by hand, once, on Sam's go (given 20 Nov 2026 for 16:00). Nothing else calls it.
#   BEGIN=16:00 bash submit_all.sh     (BEGIN: Slurm start time of the first job; everything else depends on it)
# At most 7 nodes at a time (Sam: keep it under 12):
#   wave 1: level energies and n <= 8 radial functions at the 9 spins                 3 nodes, ~30 min
#   gate  : check_levels.py must find every energy                                    seconds
#   wave 2: spin_fast.jl on all horizon channels (5 nodes), its high-spin variant for channels with an n = 8 level
#           (1 node) and the full-path validation sample (1 node), started by the gate
set -e
D=/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext/spin_tab; cd $D
sub() { sbatch --parsable "$@"; }
L=$(sub ${BEGIN:+--begin=$BEGIN} --job-name=spt_lev --array=0-2 --export=ALL,MODE=levels,NPER=64,NARR=3 node_job.sh)
G=$(sub --dependency=afterok:$L check_job.sh)
echo "$(date -Is) levels $L (begin ${BEGIN:-now})   gate $G" | tee -a logs/submitted.txt
F=$(sub --job-name=spt_fast --dependency=afterok:$G --array=0-4 --export=ALL,MODE=fast,LIST=fast_bh,NPER=96,NARR=5 node_job.sh)
H=$(sub --job-name=spt_fast7 --dependency=afterok:$G --export=ALL,MODE=fast7,LIST=fast_n8,NPER=96,NARR=1 node_job.sh)
V=$(sub --job-name=spt_check --dependency=afterok:$G --export=ALL,MODE=full,LIST=full_check,NPER=64,NARR=1 node_job.sh)
echo "fast $F   fast7 $H   full_check $V" | tee -a logs/submitted.txt
