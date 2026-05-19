#!/bin/bash 
source /mnt/users/switte/.bashrc 
cd .. 
alpha_min=0.05 
alpha_pts=15 
S1="11-8-8" 
S2="17-10-10" 
S3="988" 
S4="Inf" 
ftag="_LvrHc_" 
RL=true 
check_err=false 
srun --exclusive julia Compute_all_rates.jl --alpha_min $alpha_min --alpha_pts $alpha_pts --S1 $S1 --S2 $S2 --S3 $S3 --S4 $S4 --ftag $ftag --run_leaver $RL --check_err $check_err --use_heunc true 
wait 
S1="11-8-8" 
S2="18-10-10" 
S3="544" 
S4="Inf" 
ftag="_LvrHc_" 
RL=true 
check_err=false 
srun --exclusive julia Compute_all_rates.jl --alpha_min $alpha_min --alpha_pts $alpha_pts --S1 $S1 --S2 $S2 --S3 $S3 --S4 $S4 --ftag $ftag --run_leaver $RL --check_err $check_err --use_heunc true 
wait 
S1="11-8-8" 
S2="18-10-10" 
S3="644" 
S4="Inf" 
ftag="_LvrHc_" 
RL=true 
check_err=false 
srun --exclusive julia Compute_all_rates.jl --alpha_min $alpha_min --alpha_pts $alpha_pts --S1 $S1 --S2 $S2 --S3 $S3 --S4 $S4 --ftag $ftag --run_leaver $RL --check_err $check_err --use_heunc true 
wait 
S1="11-8-8" 
S2="18-10-10" 
S3="655" 
S4="Inf" 
ftag="_LvrHc_" 
RL=true 
check_err=false 
srun --exclusive julia Compute_all_rates.jl --alpha_min $alpha_min --alpha_pts $alpha_pts --S1 $S1 --S2 $S2 --S3 $S3 --S4 $S4 --ftag $ftag --run_leaver $RL --check_err $check_err --use_heunc true 
wait 
S1="11-8-8" 
S2="18-10-10" 
S3="744" 
S4="Inf" 
ftag="_LvrHc_" 
RL=true 
check_err=false 
srun --exclusive julia Compute_all_rates.jl --alpha_min $alpha_min --alpha_pts $alpha_pts --S1 $S1 --S2 $S2 --S3 $S3 --S4 $S4 --ftag $ftag --run_leaver $RL --check_err $check_err --use_heunc true 
wait 
S1="11-8-8" 
S2="18-10-10" 
S3="755" 
S4="Inf" 
ftag="_LvrHc_" 
RL=true 
check_err=false 
srun --exclusive julia Compute_all_rates.jl --alpha_min $alpha_min --alpha_pts $alpha_pts --S1 $S1 --S2 $S2 --S3 $S3 --S4 $S4 --ftag $ftag --run_leaver $RL --check_err $check_err --use_heunc true 
wait 
S1="11-8-8" 
S2="18-10-10" 
S3="766" 
S4="Inf" 
ftag="_LvrHc_" 
RL=true 
check_err=false 
srun --exclusive julia Compute_all_rates.jl --alpha_min $alpha_min --alpha_pts $alpha_pts --S1 $S1 --S2 $S2 --S3 $S3 --S4 $S4 --ftag $ftag --run_leaver $RL --check_err $check_err --use_heunc true 
wait 
S1="11-8-8" 
S2="18-10-10" 
S3="866" 
S4="Inf" 
ftag="_LvrHc_" 
RL=true 
check_err=false 
srun --exclusive julia Compute_all_rates.jl --alpha_min $alpha_min --alpha_pts $alpha_pts --S1 $S1 --S2 $S2 --S3 $S3 --S4 $S4 --ftag $ftag --run_leaver $RL --check_err $check_err --use_heunc true 
wait 
S1="11-8-8" 
S2="18-10-10" 
S3="877" 
S4="Inf" 
ftag="_LvrHc_" 
RL=true 
check_err=false 
srun --exclusive julia Compute_all_rates.jl --alpha_min $alpha_min --alpha_pts $alpha_pts --S1 $S1 --S2 $S2 --S3 $S3 --S4 $S4 --ftag $ftag --run_leaver $RL --check_err $check_err --use_heunc true 
wait 
S1="11-8-8" 
S2="18-10-10" 
S3="988" 
S4="Inf" 
ftag="_LvrHc_" 
RL=true 
check_err=false 
srun --exclusive julia Compute_all_rates.jl --alpha_min $alpha_min --alpha_pts $alpha_pts --S1 $S1 --S2 $S2 --S3 $S3 --S4 $S4 --ftag $ftag --run_leaver $RL --check_err $check_err --use_heunc true 
wait 
