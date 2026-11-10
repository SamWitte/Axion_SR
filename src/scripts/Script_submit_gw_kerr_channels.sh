#!/bin/bash
# Submit one Slurm job (addqueue, 1 core) per channel of a channel list file for the relativistic GW tables
# (scripts/build_gw_kerr_tables.jl --kind/--a/--b): tables land in rate_sve/gw_kerr/, logs in
# scripts/gw_kerr_logs/, per-channel wrappers in scripts/gw_kerr_jobs/.
#   bash scripts/Script_submit_gw_kerr_channels.sh LIST [extra build_gw_kerr_tables.jl args]
# LIST has one channel per line: "kind a b" (kind = tr | ann, levels as n-l-m); # comments allowed.
# Annihilations are computed for alpha >= 0.05, transitions for alpha >= 0.03 (as the default set).
# Rerunning resumes unfinished tables.
source /mnt/users/switte/.bashrc
cd /mnt/users/switte/Axion_SR_gw_v2/src
LIST=$1; shift; EXTRA="$*"
[ -f "$LIST" ] || { echo "usage: $0 LIST [extra args]"; exit 1; }
mkdir -p scripts/gw_kerr_jobs scripts/gw_kerr_logs rate_sve/gw_kerr
IDS=scripts/gw_kerr_jobs/ids_$(date +%Y%m%d_%H%M).txt
n=0
while read -r kind a b; do
  case "$kind" in ''|'#'*) continue;; esac
  amin=0.03; [ "$kind" = "ann" ] && amin=0.05
  scr=scripts/gw_kerr_jobs/Script_gwk_${kind}_${a}_${b}.sh
  cat > $scr <<EOF
#!/bin/bash
source /mnt/users/switte/.bashrc
cd /mnt/users/switte/Axion_SR_gw_v2/src
julia -t 1 scripts/build_gw_kerr_tables.jl --kind $kind --a $a --b $b --alpha_min $amin $EXTRA
EOF
  chmod +x $scr
  out=$(addqueue -s -n 1x1 -c "3 days" -m 6 -o scripts/gw_kerr_logs/gwk_${kind}_${a}_${b}-%j.out ./$scr 2>&1)
  id=$(echo "$out" | grep -oE "Script_gwk_${kind}_${a}_${b}.sh-[0-9]{7,}" | grep -oE '[0-9]{7,}$')
  [ -n "$id" ] || echo "submit failed for $kind $a $b: $out" >&2
  echo "$kind $a $b $id" | tee -a $IDS
  n=$((n + 1))
done < "$LIST"
echo "$n channels submitted; ids in $IDS"
