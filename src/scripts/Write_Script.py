import numpy as np
import os

num_scripts = 60 # num scripts
nstart = 700
vstart = 700
script_dir = "scripts/"
script_tag = "_CygC_rerun_N5_"

arr_text = np.empty(num_scripts, dtype=object)
    
for i in range(num_scripts):
    arr_text[i] = "#!/bin/bash \n"
    arr_text[i] += "source /mnt/users/switte/.bashrc \n"
    arr_text[i] += "cd .. \n"
    arr_text[i] += "alpha_max_cut=10.0 \n"
    arr_text[i] += "alpha_min_cut=0.001 \n"
    arr_text[i] += "tau_max_override=5.0e7 \n"
    arr_text[i] += "Ftag=\"_v{}_\" \n".format(i + vstart)
    arr_text[i] += "input_data=\"Me\" \n" ### Check this
    arr_text[i] += "one_BH=\"CygX1_Conserve\" \n" ### Check this
    arr_text[i] += "numwalkers=4 \n" ## Check this
    arr_text[i] += "numsamples_perwalker=500 \n"
    arr_text[i] += "solve_n4=true \n" ### Check this
    arr_text[i] += "solve_n5=true \n" ### Check this
    arr_text[i] += "burnin=300 \n"
    arr_text[i] += "stop_on_a=0.7 \n" ### Check this
    arr_text[i] += "non_rel=false \n" ### Check this
    arr_text[i] += "srun --exclusive julia MCMC.jl --alpha_max_cut $alpha_max_cut --alpha_min_cut $alpha_min_cut --tau_max_override $tau_max_override --Ftag $Ftag --input_data $input_data --one_BH $one_BH --numwalkers $numwalkers --numsamples_perwalker $numsamples_perwalker --solve_n4 $solve_n4 --stop_on_a $stop_on_a --burnin $burnin --non_rel $non_rel --solve_n5 $solve_n5 \n"
   


for i in range(num_scripts):
    text_file = open("Script_Run_"+script_tag+"RT_{:.0f}.sh".format(i + nstart), "w")
    text_file.write(arr_text[i])
    text_file.close()


    os.system("chmod u+x Script_Run_"+script_tag+"RT_{:.0f}.sh".format(i + nstart))
#    os.system("cd "+script_dir)
    os.system("addqueue -s -n 1 -c \"1 week\" -m 30 ./Script_Run_"+script_tag+"RT_{:.0f}.sh".format(i + nstart))
