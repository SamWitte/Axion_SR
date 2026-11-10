#!/bin/bash
source /mnt/users/switte/.bashrc
cd /mnt/users/switte/Axion_SR_gw_v2/src
julia -t 1 scripts/diag_gw_kerr_bigfloat.jl
