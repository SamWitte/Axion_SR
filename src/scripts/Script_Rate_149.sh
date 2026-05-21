#!/bin/bash 
source /mnt/users/switte/.bashrc 
cd .. 
alpha_min=0.05 
alpha_pts=15 
S1="21-20-20" 
S2="21-20-20" 
S3="13-10-10" 
S4="Inf" 
ftag="_LvrHc_" 
RL=true 
check_err=false 
srun --exclusive julia Compute_all_rates.jl --alpha_min $alpha_min --alpha_pts $alpha_pts --S1 $S1 --S2 $S2 --S3 $S3 --S4 $S4 --ftag $ftag --run_leaver $RL --check_err $check_err --use_heunc true 
wait 
S1="21-20-20" 
S2="21-20-20" 
S3="14-10-10" 
S4="Inf" 
ftag="_LvrHc_" 
RL=true 
check_err=false 
srun --exclusive julia Compute_all_rates.jl --alpha_min $alpha_min --alpha_pts $alpha_pts --S1 $S1 --S2 $S2 --S3 $S3 --S4 $S4 --ftag $ftag --run_leaver $RL --check_err $check_err --use_heunc true 
wait 
