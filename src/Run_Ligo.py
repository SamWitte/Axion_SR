import numpy as np
import os

massL = np.logspace(-14, -12, 10)
chainsN = 5
arr_text = np.empty(len(massL), dtype=object)

for i in range(len(massL)):
    arr_text[i] = "#!/bin/bash \n"
    arr_text[i] += "source /mnt/users/switte/.bashrc \n"
    arr_text[i] += "cd .. \n"
    arr_text[i] += "dataname=\"LIGO_test\" \n"
    arr_text[i] += "ax_mass={:.3e} \n".format(massL[i])
    arr_text[i] += "tau_max_override=4.5e7 \n"
    arr_text[i] += "delt_M=0.05 \n"
    arr_text[i] += "Nmax=3 \n"
    arr_text[i] += "numsamples_perwalker=1000 \n"
    arr_text[i] += "burnin=500 \n"
    arr_text[i] += "Ftag=\"_\" \n"
    
    arr_text[i] += "srun --exclusive julia --threads " + str(chainsN) + " Ligo_1dL.jl --dataname $dataname --ax_mass $ax_mass --tau_max_override $tau_max_override --delt_M $delt_M --Nmax $Nmax --numsamples_perwalker $numsamples_perwalker --burnin $burnin --Ftag $Ftag \n"
   

for i in range(len(massL)):
   
    text_file = open("Script_Ligo_{:.0f}.sh".format(i), "w")
    text_file.write(arr_text[i])
    text_file.close()


    os.system("chmod u+x Script_Ligo_{:.0f}.sh".format(i))
    os.system("addqueue -s -n " + str(chainsN) + " -c \"1 week\" -m 30 ./Script_Ligo_{:.0f}.sh".format(i))
