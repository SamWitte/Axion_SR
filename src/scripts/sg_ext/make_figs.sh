#!/bin/bash
# Spin + level-occupation figures of the extended-SG runs (plot script of the first campaign, sg_runs/plot_sg_states.py):
# per mass and alpha, SG_only vs noSG_noSI and SG_SI_fa_<f> vs SG_only (solid vs dashed).
#   ./make_figs.sh <set: p1|core|full> <gw: nonrel|rel> [eps_min=1e-14]
X=/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext
SET=$1; GW=$2; EMIN=${3:-1e-14}
PY=/software/python/3.12/bin/python3
cd /data/dust/user/samwitte/Axion_SR_gw_scan/sg_runs
for M in 10 10000; do
  O=$X/runs/$SET/gw_$GW/M_$M; F=$X/figs/$SET/gw_$GW/M_$M; mkdir -p $F
  plot() { [ -n "$(ls $O/$2/alpha_$1/Time_*.dat 2>/dev/null)" ] && [ -n "$(ls $O/$3/alpha_$1/Time_*.dat 2>/dev/null)" ] || return
           $PY plot_sg_states.py --alpha $1 --bh $M --tmax 1e9 --eps-min $EMIN --max-modes 24 --solid $2 --solid-dir $O/$2/alpha_$1 \
             --dashed $3 --dashed-dir $O/$3/alpha_$1 --out $F/states_M${M}_alpha$1_$2_vs_$3.png > $F/.log_$1_$2_vs_$3.txt 2>&1; }
  for a in 0.05 0.075 0.1 0.15 0.2 0.25; do
    plot $a SG_only noSG_noSI &
    for f in 1e16 5e16 1e17; do plot $a SG_SI_fa_$f SG_only & done
    wait
  done
  echo "M=$M: $(ls $F/*.png 2>/dev/null | wc -l) figures in $F"
done
