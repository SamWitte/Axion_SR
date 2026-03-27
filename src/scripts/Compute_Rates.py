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

# Filter out states whose rate files already exist
states_to_compute = []
for i in range(len(S1)):
    filename = "../rate_sve/" + S1[i] + "_" + S2[i] + "_" + S3[i] + "_" + S4[i] + ftag + ".dat"
    if not os.path.isfile(filename):
        states_to_compute.append(i)

print(f"Total states: {len(S1)}, Already computed: {len(S1) - len(states_to_compute)}, To compute: {len(states_to_compute)}")

arr_text = np.empty(len(states_to_compute), dtype=object)
for idx, i in enumerate(states_to_compute):

    arr_text[idx] = "S1=\""+S1[i]+"\" \n"
    arr_text[idx] += "S2=\""+S2[i]+"\" \n"
    arr_text[idx] += "S3=\""+S3[i]+"\" \n"
    arr_text[idx] += "S4=\""+S4[i]+"\" \n"
    arr_text[idx] += "ftag=\""+ftag+"\" \n"
    if run_leaver:
        arr_text[idx] += "RL=true \n"
    else:
        arr_text[idx] += "RL=false \n"
    if check_err:
        arr_text[idx] += "check_err=true \n"
    else:
        arr_text[idx] += "check_err=false \n"
    arr_text[idx] += "srun --exclusive julia Compute_all_rates.jl --alpha_min $alpha_min --alpha_pts $alpha_pts --S1 $S1 --S2 $S2 --S3 $S3 --S4 $S4 --ftag $ftag --run_leaver $RL --check_err $check_err \n"

cnt = 0
n_per_file = int(math.ceil(len(states_to_compute) / float(N_groups)))
def split_into_groups(lst, n):
    k, m = divmod(len(lst), n)
    return [lst[i*k + min(i, m):(i+1)*k + min(i+1, m)] for i in range(n)]


groups = split_into_groups(list(range(len(states_to_compute))), N_groups)

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

