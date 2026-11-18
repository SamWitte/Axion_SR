#!/bin/bash
# Watch the exact-exchange production: release stuck claims every 10 minutes and return when all points are done,
# when fewer than MINW claims are active (workers died), or after MAXMIN minutes.   usage: watch_prod.sh [MAXMIN] [MINW]
MAXMIN=${1:-90}; MINW=${2:-8}
cd /data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext/sg_exact || exit 1
t0=$(date +%s)
while true; do
  ./reap_locks.sh prod2 prod_dip prod > /dev/null
  st=$(/software/python/3.12/bin/python3 status.py | grep '^total')
  miss=$(echo "$st" | awk '{print $NF}')
  act=0; for d in prod prod_dip prod2; do for l in $(find $d/locks -mindepth 1 -maxdepth 1 -type d -mmin -50); do [ -f $d/res/$(basename $l).txt ] || act=$((act+1)); done; done
  el=$(( ( $(date +%s) - t0 ) / 60 ))
  echo "$(date +%H:%M) $st | active claims $act | elapsed $el min"
  [ "$miss" = "0" ] && { echo "ALL DONE"; break; }
  [ $act -lt $MINW ] && { echo "FEW WORKERS"; break; }
  [ $el -ge $MAXMIN ] && { echo "TIME"; break; }
  sleep 600
done
squeue -u samwitte -h -o "%i %j %T %M %R" | grep sgx
for n in $(squeue -u samwitte -h -n sgx_exact -t R -o %N); do echo "$n $(ssh -o BatchMode=yes -o ConnectTimeout=5 $n 'ps -u samwitte -o args | grep -c "[j]ulia-1.8.0/bin/julia.*scan_gen"; free -g | awk "NR==2{print \"free GB\",\$4}"' | tr '\n' ' ')"; done
echo "n88: point files $(ls ../rates_out 2>/dev/null | wc -l)"
mmlsquota -j user.it.samwitte --block-size auto dust 2>/dev/null | tail -1
