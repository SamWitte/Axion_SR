#!/bin/bash
# Submit one Slurm job (addqueue, 1 core) per channel of the relativistic GW table set
# (scripts/build_gw_kerr_tables.jl --default): tables land in rate_sve/gw_kerr/, logs in
# scripts/gw_kerr_logs/, per-channel wrappers in scripts/gw_kerr_jobs/.
#   bash scripts/Script_submit_gw_kerr_tables.sh [LMAX_LEVELS=5] [OVERTONES=1] [extra build_gw_kerr_tables.jl args]
# Annihilations are computed for alpha >= 0.05 (Float64 is accurate there; the ratio tends to a
# constant below), transitions for alpha >= 0.03. Rerunning resumes unfinished tables.
source /mnt/users/switte/.bashrc
cd /mnt/users/switte/Axion_SR_gw_v2/src
LMAX=${1:-5}; OVER=${2:-1}; shift 2 2>/dev/null; EXTRA="$*"
mkdir -p scripts/gw_kerr_jobs scripts/gw_kerr_logs rate_sve/gw_kerr
LIST=scripts/gw_kerr_jobs/channels_l${LMAX}_o${OVER}.txt
julia scripts/build_gw_kerr_tables.jl --default --lmax_levels $LMAX --overtones $OVER --list > $LIST || { echo "channel listing failed"; exit 1; }
echo "$(wc -l < $LIST) channels"
IDS=scripts/gw_kerr_jobs/ids_$(date +%Y%m%d_%H%M).txt
while read -r i kind a b; do
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
  echo "$kind $a $b $id" | tee -a $IDS
done < $LIST
echo "ids in $IDS"
