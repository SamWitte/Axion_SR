#!/bin/bash
# Release exact-exchange jobs whose worker died: lock older than 50 minutes, no result (a job takes 3-15 minutes).
# usage: reap_locks.sh <outdir> ...
cd /data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext/sg_exact || exit 1
for OUT in "$@"; do
  n=0
  for d in $(find $OUT/locks -mindepth 1 -maxdepth 1 -type d -mmin +50 2>/dev/null); do
    id=$(basename $d)
    [ -f $OUT/res/$id.txt ] && continue
    [ "$(cut -d' ' -f1 $d/owner 2>/dev/null)" = "held_for_Lextra" ] && continue
    mv $d $OUT/stale_$(date +%s%N)_$id && n=$((n+1))
  done
  echo "$OUT: released $n stuck claims"
done
