import numpy as np
import os

massL = np.logspace(np.log10(3e-14), -12, 40)
chainsN = 2
arr_text = np.empty(len(massL), dtype=object)
tg="_v0"
for i in range(len(massL)):
    arr_text[i] = "#!/bin/bash \n"
    arr_text[i] += "source /mnt/users/switte/.bashrc \n"
    # arr_text[i] += "cd .. \n"
    arr_text[i] += "dataname=\"GW231123\" \n" # GW231123_NRSur
    arr_text[i] += "ax_mass={:.3e} \n".format(massL[i])
    arr_text[i] += "tau_max_override=1e5 \n"
    arr_text[i] += "delt_M=0.05 \n"
    arr_text[i] += "Nmax=5 \n"
    arr_text[i] += "numsamples_perwalker=500 \n"
    arr_text[i] += "burnin=200 \n"
    arr_text[i] += "Ftag=\"_4dKDE_nr_\" \n"
    arr_text[i] += "oneBH=false \n"
    arr_text[i] += "cuthigh=false \n"
    
    arr_text[i] += "srun --exclusive julia --threads " + str(chainsN) + " Ligo_1dL.jl --dataname $dataname --ax_mass $ax_mass --tau_max_override $tau_max_override --delt_M $delt_M --Nmax $Nmax --numsamples_perwalker $numsamples_perwalker --burnin $burnin --Ftag $Ftag --one_BH $oneBH --cut_high_spin $cuthigh \n"
   

for i in range(len(massL)):
   
    text_file = open("Script_Ligo" + tg + "_{:.0f}.sh".format(i), "w")
    text_file.write(arr_text[i])
    text_file.close()


    os.system("chmod u+x Script_Ligo" + tg + "_{:.0f}.sh".format(i))
    os.system("addqueue -s -n " + str(chainsN) + " -c \"1 week\" -m 6 ./Script_Ligo" + tg + "_{:.0f}.sh".format(i))
