#!/bin/bash
#SBATCH --partition=allcpu --constraint="75F3|7542|7402|7H12" --nodes=1 --time=00:10:00 --requeue
#SBATCH --job-name=spt_chk --output=/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext/spin_tab/logs/check_%j.out
# Gate between the level energies and everything that needs them: fails (and holds the dependent jobs) if rows are missing.
python3 /data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext/spin_tab/check_levels.py
