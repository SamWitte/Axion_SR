#!/bin/bash
#SBATCH --partition=maxcpu --constraint="7402|75F3|7H12" --nodes=1 --time=0-06:00:00
#SBATCH --job-name=sgx_frame --output=/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext/sg_exact/test/frame_check/job_%j.out
# Check of the frame / recoil treatment of the relativistic exchange (exchange_gen.jl, production path via scan_gen.jl):
#   free  : the normal calculation (struts cancelled, hole recoils freely) at low alpha for recoil-dominated channels
#   fixed : SGX_HOLD_FIXED=1, the hole held at rest by the struts  -> must reproduce the Newtonian amplitude WITHOUT recoil
source ~/.bashrc
D=/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext/sg_exact/test/frame_check; cd $D
S=/data/dust/user/samwitte/Axion_SR_sg/src/ongoing_tests/sg_relativistic/scan_gen.jl; J=~/julia-1.8.0/bin/julia
export JULIA_NUM_THREADS=12 OPENBLAS_NUM_THREADS=12
for w in 1 2 3; do ( $J $S jobs_free.txt free > free/logs/w$w.log 2>&1 ) & sleep 5; done
for w in 1 2 3 4; do ( SGX_HOLD_FIXED=1 $J $S jobs_fixed.txt fixed > fixed/logs/w$w.log 2>&1 ) & sleep 5; done
wait; echo "end $(date)"
