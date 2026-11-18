#!/bin/bash
# The SI-only runs of submit_si_only.sh (version sgtabSI: M = 10 Msun, alpha 0.1 and 0.15, f_a = 1e16, 5e16, 1e17, 1e18,
# non-relativistic GW, --self_grav false, NL_SI_FROM_SG=1) repeated with Gamma_SI from the Kerr-state overlay tables.
#   si_only_kerr.sh <jobid of a running allocation> [tables dir = tables3s_kerr_nr] [set name = si_only_kerr]
# Output: ../runs/<set>/gw_nonrel/M_10/<variant>/alpha_<a>/
JID=$1; TAB=${2:-tables3s_kerr_nr}; SET=${3:-si_only_kerr}
X=/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext
KL=$(paste -sd, $X/levels_sgx.txt)
L=$X/logs/runs/$SET; mkdir -p $L; T=$X/sg_exact/run_tasks_$SET.sh; : > $T
for a in 0.1 0.15; do
  for f in 1e16 5e16 1e17 1e18; do
    out=$X/runs/$SET/gw_nonrel/M_10/SI_only_fa_$f/alpha_$a
    [ -n "$(ls $out/Spin_*.dat 2>/dev/null)" ] && continue
    mkdir -p $out
    echo "~/julia-1.8.0/bin/julia run_Nlevels.jl --MassBH 10 --SpinBH 0.95 --f_a $f --alpha $a --Nmax 45 --tau_max 1e9 --gw_model nonrel --self_grav false --rate_input sgx_channels.txt --keep_levels \"$KL\" --gw_table_nmax sgx --regime_split true --solver radau --analytic_jac --timescale_saves false --maxiters 1e8 --n_times 20000 --outdir $out >> $L/SI_only_fa_${f}_alpha_$a.log 2>&1" >> $T
  done
done
n=$(wc -l < $T); [ $n -gt 0 ] || { echo "nothing to do"; exit 0; }
nohup srun --jobid=$JID --overlap -n1 -c$n bash -c "source ~/.bashrc; cd /data/dust/user/samwitte/Axion_SR_sg/src/Nlevels_Runs
  export JULIA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NL_PROGRESS_SEC=1800 NL_SI_FROM_SG=1 NL_SG_OVERLAY=$X/sg_exact/$TAB
  pids=(); while IFS= read -r line; do [ -z \"\$line\" ] && continue; bash -c \"\$line\" & pids+=(\$!); sleep 1; done < $T
  for p in \"\${pids[@]}\"; do wait \$p; done; echo finished \$(date)" > $L/overlap_$JID.out 2>&1 &
echo "started $n SI-only runs with $TAB inside job $JID"
