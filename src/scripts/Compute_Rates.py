import numpy as np
import os

nstart = 1

S1 = []
S2 = []
S3 = []
S4 = []

file1 = open("load_compute_rates_new.dat", 'r')
Lines = file1.readlines()

count = 0
for line in Lines:
    line_parts = line.split()
    S1.append(line_parts[0])
    S2.append(line_parts[1])
    S3.append(line_parts[2])
    S4.append(line_parts[3])
     

ftag="_vt2_"

arr_text = np.empty(len(S1), dtype=object)
for i in range(len(S1)):
    arr_text[i] = "#!/bin/bash \n"
    arr_text[i] += "source /mnt/users/switte/.bashrc \n"
    arr_text[i] += "cd .. \n"
    arr_text[i] += "alpha_max=1.0 \n"
    arr_text[i] += "alpha_min=0.03 \n"
    arr_text[i] += "alpha_pts=20 \n"
    arr_text[i] += "S1=\""+S1[i]+"\" \n"
    arr_text[i] += "S2=\""+S2[i]+"\" \n"
    arr_text[i] += "S3=\""+S3[i]+"\" \n"
    arr_text[i] += "S4=\""+S4[i]+"\" \n"
    arr_text[i] += "ftag=\""+ftag+"\" \n"
    arr_text[i] += "srun --exclusive julia Compute_all_rates.jl --alpha_max $alpha_max --alpha_min $alpha_min --alpha_pts $alpha_pts --S1 $S1 --S2 $S2 --S3 $S3 --S4 $S4 --ftag $ftag \n"
   

for i in range(len(S1)):
    text_file = open("Script_Rate_{:.0f}.sh".format(i), "w")
    text_file.write(arr_text[i])
    text_file.close()


    os.system("chmod u+x Script_Rate_{:.0f}.sh".format(i))
#    os.system("cd "+script_dir)
    os.system("addqueue -s -n 1 -c \"1 week\" -m 10 ./Script_Rate_{:.0f}.sh".format(i))
