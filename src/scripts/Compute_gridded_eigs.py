import numpy as np
import os
import time


runL=False
zeroS=True
gridD=False

if runL and zeroS:
    fnOut= "../rate_sve/Imag_zero_"
elif not runL and zeroS:
    fnOut= "../rate_sve/Imag_zeroC_"
elif runL and gridD:
    fnOut= "../rate_sve/Imag_erg_neg_"
elif not runL and gridD:
    fnOut= "../rate_sve/Imag_ergC_neg_"

    
N_start = 1

nmax=8
listnlm = []
for nn in range(2, nmax+1):
    for ll in range(1, nn):
        for mm in range(1, ll+1):
            listnlm.append(str(nn)+str(ll)+str(mm))
            
arr_text = np.empty(len(listnlm), dtype=object)
for i in range(len(listnlm)):
    arr_text[i] = "#!/bin/bash \n"
    arr_text[i] += "source /mnt/users/switte/.bashrc \n"
    arr_text[i] += "cd .. \n"
    if runL:
        arr_text[i] += "RunLeaver=true \n"
    else:
        arr_text[i] += "RunLeaver=false \n"
    if zeroS:
        arr_text[i] += "zerosS=true \n"
    else:
        arr_text[i] += "zerosS=false \n"

    if gridD:
        arr_text[i] += "grided=true \n"
    else:
        arr_text[i] += "grided=false \n"
        
    arr_text[i] += "nlm=\""+listnlm[i]+"\" \n"
    
    arr_text[i] += "srun --exclusive julia grid_sr_growth.jl --run_leaver $RunLeaver --solve_for_zeros $zerosS --solve_gridded $grided --nlmIn $nlm  \n"
   
if zeroS:
    snme = "Script_Grid_"
elif gridD:
    snme = "Script_FLGrid_"
    
for i in range(len(listnlm)):
    filecheck = fnOut +listnlm[i]+".dat"
    if os.path.isfile(filecheck):
        continue
    print(filecheck, "\t", os.path.isfile(filecheck))
    text_file = open(snme + "{:.0f}.sh".format(i + N_start), "w")
    text_file.write(arr_text[i])
    text_file.close()


    os.system("chmod u+x "+snme + "{:.0f}.sh".format(i + N_start))
#    os.system("cd "+script_dir)
    os.system("addqueue -s -n 1 -c \"1 week\" -m 40 ./"+snme+"{:.0f}.sh".format(i + N_start))
    time.sleep(3)

        
