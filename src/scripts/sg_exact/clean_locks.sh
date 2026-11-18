#!/bin/bash
# Release jobs claimed by workers whose Slurm job is gone (no result file, owner job not in the queue).
# usage: clean_locks.sh [outdir]      (default prod)
OUT=${1:-prod}
cd /data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext/sg_exact/$OUT/locks || exit 1
live=$(squeue -u $USER -h -o %i | tr '\n' ' ')
n=0
for d in */; do
  d=${d%/}
  [ -f ../res/$d.txt ] && continue
  [ -f ../fail/$d.txt ] && continue
  o=$(cut -d' ' -f1 $d/owner 2>/dev/null)
  [ "$o" = "held_for_Lextra" ] && continue
  case " $live " in *" $o "*) continue;; esac
  mv $d ../stale_$(date +%s)_$d 2>/dev/null && n=$((n+1))
done
echo "released $n stale claims"
