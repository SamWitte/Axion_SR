import numpy as np
import os

massL = np.logspace(np.log10(3e-14), np.log10(7e-13), 30)
massL2 = np.logspace(np.log10(7.1e-13), np.log10(2e-12), 10)
chainsN = 10
ft="_"

for i in range(len(massL)):
    script_name = f"run_Ligo_oM_{i}.sh"
    submit_name = f"submit_Ligo_oM_{i}.sub"
    
    # Create wrapper shell script
    with open(script_name, "w") as script_file:
        script_file.write("#!/bin/bash\n")
        script_file.write("source ~/.bash_profile\n")
        script_file.write("cd /afs/cern.ch/user/g/gafranci/sjw/Axion_SR_TDE-main/src/ \n")
        script_file.write("dataname=\"GW231123_NRSur\"\n")
        script_file.write(f"ax_mass={massL[i]:.3e}\n")
        script_file.write("tau_max_override=1.0e6\n")
        script_file.write("delt_M=0.05\n")
        script_file.write("Nmax=3\n")
        script_file.write("numsamples_perwalker=2000\n")
        script_file.write("burnin=50\n")
        script_file.write("Ftag=\"_\"\n")
        script_file.write(f"julia --threads {chainsN} Ligo_1dL.jl "
                          "--dataname $dataname --ax_mass $ax_mass --tau_max_override $tau_max_override "
                          "--delt_M $delt_M --Nmax $Nmax --numsamples_perwalker $numsamples_perwalker "
                          "--burnin $burnin --Ftag $Ftag --one_BH true \n")

    os.chmod(script_name, 0o744)

    # Create HTCondor submit file
    with open(submit_name, "w") as submit_file:
        submit_file.write(f"executable = {script_name}\n")
        submit_file.write(f"arguments = $(ClusterId) $(ProcId) \n")
        submit_file.write(f"output = Ligo_{i}.out\n")
        submit_file.write(f"error = Ligo_{i}.err\n")
        submit_file.write(f"log = Ligo_{i}.log\n")
        submit_file.write("RequestCpus = {}\n".format(chainsN))
        # submit_file.write("+JobFlavour = testmatch \n")  # 1 week in seconds
        submit_file.write("+MaxRuntime = 259200 \n")
        submit_file.write("queue\n")

    os.system('condor_submit ' + submit_name)
