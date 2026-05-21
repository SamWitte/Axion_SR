#!/bin/bash 
source /mnt/users/switte/.bashrc 
cd .. 
alpha_min=0.05 
alpha_pts=15 
S1="744" 
S2="17-10-10" 
S3="544" 
S4="Inf" 
ftag="_LvrHc_" 
RL=true 
check_err=false 
srun --exclusive julia Compute_all_rates.jl --alpha_min $alpha_min --alpha_pts $alpha_pts --S1 $S1 --S2 $S2 --S3 $S3 --S4 $S4 --ftag $ftag --run_leaver $RL --check_err $check_err --use_heunc true 
wait 
S1="744" 
S2="17-10-10" 
S3="644" 
S4="Inf" 
ftag="_LvrHc_" 
RL=true 
check_err=false 
srun --exclusive julia Compute_all_rates.jl --alpha_min $alpha_min --alpha_pts $alpha_pts --S1 $S1 --S2 $S2 --S3 $S3 --S4 $S4 --ftag $ftag --run_leaver $RL --check_err $check_err --use_heunc true 
wait 
S1="744" 
S2="17-10-10" 
S3="655" 
S4="Inf" 
ftag="_LvrHc_" 
RL=true 
check_err=false 
srun --exclusive julia Compute_all_rates.jl --alpha_min $alpha_min --alpha_pts $alpha_pts --S1 $S1 --S2 $S2 --S3 $S3 --S4 $S4 --ftag $ftag --run_leaver $RL --check_err $check_err --use_heunc true 
wait 
S1="744" 
S2="18-10-10" 
S3="544" 
S4="Inf" 
ftag="_LvrHc_" 
RL=true 
check_err=false 
srun --exclusive julia Compute_all_rates.jl --alpha_min $alpha_min --alpha_pts $alpha_pts --S1 $S1 --S2 $S2 --S3 $S3 --S4 $S4 --ftag $ftag --run_leaver $RL --check_err $check_err --use_heunc true 
wait 
S1="744" 
S2="18-10-10" 
S3="655" 
S4="Inf" 
ftag="_LvrHc_" 
RL=true 
check_err=false 
srun --exclusive julia Compute_all_rates.jl --alpha_min $alpha_min --alpha_pts $alpha_pts --S1 $S1 --S2 $S2 --S3 $S3 --S4 $S4 --ftag $ftag --run_leaver $RL --check_err $check_err --use_heunc true 
wait 
S1="755" 
S2="755" 
S3="16-10-10" 
S4="BH" 
ftag="_LvrHc_" 
RL=true 
check_err=false 
srun --exclusive julia Compute_all_rates.jl --alpha_min $alpha_min --alpha_pts $alpha_pts --S1 $S1 --S2 $S2 --S3 $S3 --S4 $S4 --ftag $ftag --run_leaver $RL --check_err $check_err --use_heunc true 
wait 
S1="755" 
S2="755" 
S3="17-10-10" 
S4="BH" 
ftag="_LvrHc_" 
RL=true 
check_err=false 
srun --exclusive julia Compute_all_rates.jl --alpha_min $alpha_min --alpha_pts $alpha_pts --S1 $S1 --S2 $S2 --S3 $S3 --S4 $S4 --ftag $ftag --run_leaver $RL --check_err $check_err --use_heunc true 
wait 
S1="755" 
S2="755" 
S3="18-10-10" 
S4="BH" 
ftag="_LvrHc_" 
RL=true 
check_err=false 
srun --exclusive julia Compute_all_rates.jl --alpha_min $alpha_min --alpha_pts $alpha_pts --S1 $S1 --S2 $S2 --S3 $S3 --S4 $S4 --ftag $ftag --run_leaver $RL --check_err $check_err --use_heunc true 
wait 
S1="755" 
S2="855" 
S3="16-10-10" 
S4="BH" 
ftag="_LvrHc_" 
RL=true 
check_err=false 
srun --exclusive julia Compute_all_rates.jl --alpha_min $alpha_min --alpha_pts $alpha_pts --S1 $S1 --S2 $S2 --S3 $S3 --S4 $S4 --ftag $ftag --run_leaver $RL --check_err $check_err --use_heunc true 
wait 
S1="755" 
S2="855" 
S3="17-10-10" 
S4="BH" 
ftag="_LvrHc_" 
RL=true 
check_err=false 
srun --exclusive julia Compute_all_rates.jl --alpha_min $alpha_min --alpha_pts $alpha_pts --S1 $S1 --S2 $S2 --S3 $S3 --S4 $S4 --ftag $ftag --run_leaver $RL --check_err $check_err --use_heunc true 
wait 
