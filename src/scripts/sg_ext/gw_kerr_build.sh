#!/bin/bash
# Kerr (Teukolsky) GW correction tables for channels listed in channels.txt ("kind a b"), one process per channel.
source ~/.bashrc
B=/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext/gw_kerr_build
cd /data/dust/user/samwitte/Axion_SR_sg/src
while read -r kind a b; do
  case "$kind" in ''|'#'*) continue;; esac
  amin=0.03; [ "$kind" = "ann" ] && amin=0.05
  [ -e $B/logs/${kind}_${a}_${b}.done ] && continue
  ( ~/julia-1.8.0/bin/julia -t 1 scripts/build_gw_kerr_tables.jl --kind $kind --a $a --b $b --alpha_min $amin > $B/logs/${kind}_${a}_${b}.log 2>&1 && touch $B/logs/${kind}_${a}_${b}.done ) &
done < ${1:-$B/channels.txt}
wait
