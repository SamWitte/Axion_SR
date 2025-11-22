import glob
import numpy as np

file_out = "CygX1_N4_Rel.dat"
dir_out = "output_mcmc"

fileList = glob.glob(dir_out +"/LowMassRegion_TauMax_5.0e7_alpha_maxmin_10.0_0.001_Me_v*__Cyg_X1_Nsamps_500__add_n4__FullRel__mcmc.dat")
for i in range(len(fileList)):
    if i == 0:
        outputall = np.loadtxt(fileList[i])
    else:
        hold = np.loadtxt(fileList[i])
        outputall = np.vstack((outputall, hold))
        
np.savetxt(dir_out + "/" + file_out, outputall)
