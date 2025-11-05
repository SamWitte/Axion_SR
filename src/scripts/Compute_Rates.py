import numpy as np
import os
import math

nstart = 1
N_groups = 100
minIX=0 # this is an index which allows you to append to load_rate_input and run only N > minIX  
S1 = []
S2 = []
S3 = []
S4 = []

run_leaver=False
check_err=False
file1 = open("load_rate_input_Nmax_5.txt", 'r')
Lines = file1.readlines()

count = 0
for line in Lines:
    line_parts = line.split()
    S1.append(line_parts[0])
    S2.append(line_parts[1])
    S3.append(line_parts[2])
    S4.append(line_parts[3])
     

ftag="_GF_Cheby_"

arr_text = np.empty(len(S1), dtype=object)
for i in range(len(S1)):

    arr_text[i] = "S1=\""+S1[i]+"\" \n"
    arr_text[i] += "S2=\""+S2[i]+"\" \n"
    arr_text[i] += "S3=\""+S3[i]+"\" \n"
    arr_text[i] += "S4=\""+S4[i]+"\" \n"
    arr_text[i] += "ftag=\""+ftag+"\" \n"
    if run_leaver:
        arr_text[i] += "RL=true \n"
    else:
        arr_text[i] += "RL=false \n"
    if check_err:
        arr_text[i] += "check_err=true \n"
    else:
        arr_text[i] += "check_err=false \n"
    arr_text[i] += "srun --exclusive julia Compute_all_rates.jl --alpha_min $alpha_min --alpha_pts $alpha_pts --S1 $S1 --S2 $S2 --S3 $S3 --S4 $S4 --ftag $ftag --run_leaver $RL --check_err $check_err \n"

cnt = 0
n_per_file = int(math.ceil(len(S1) / float(N_groups)))
def split_into_groups(lst, n):
    k, m = divmod(len(lst), n)
    return [lst[i*k + min(i, m):(i+1)*k + min(i+1, m)] for i in range(n)]

  
groups = split_into_groups(list(range(len(S1))), N_groups)

for j in range(N_groups):
    text_file = open("Script_Rate_{:.0f}.sh".format(j), "w")
    text_file.write("#!/bin/bash \n")
    text_file.write("source /mnt/users/switte/.bashrc \n")
    text_file.write("cd .. \n")
    text_file.write("alpha_min=0.03 \n")
    text_file.write("alpha_pts=30 \n")

    for i in range(len(groups[j])):
        if cnt < len(arr_text):
            text_file.write(arr_text[cnt])
            text_file.write("wait \n")
            cnt += 1
    text_file.close()

    os.system("chmod u+x Script_Rate_{:.0f}.sh".format(j))
#    os.system("cd "+script_dir)                                                                                                                                                                                                               
    os.system("addqueue -s -n 1 -c \"1 week\" -m 40 ./Script_Rate_{:.0f}.sh".format(j))

