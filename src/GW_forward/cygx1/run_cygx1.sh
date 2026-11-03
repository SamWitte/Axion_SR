#!/bin/bash
# Cyg X-1 forward model + zoom, one (alpha, f_a) per call:
#
#   ./run_cygx1.sh <alpha> <f_a_GeV> [Nmax] [rounds] [N_round]
#
# e.g. the three test cases (one job each, single-threaded; submit them separately):
#   ./run_cygx1.sh 0.1 1e18 3 1 60
#   ./run_cygx1.sh 0.2 1e18 3 1 60
#   ./run_cygx1.sh 0.5 1e18 4 3 100
#
# Measurement (Miller-Jones et al. 2021): M = 21.2 +- 2.2 M_sun today, spin today
# > 0.85 (hard cut), age 4.8-7.6 Myr (log-uniform), d = 2.22 +- 0.18 kpc.
# alpha = G M m_a is defined at 21.2 M_sun. Birth prior: M0 in [12, 35] M_sun,
# a0 in [0, 0.998], f_Edd log-uniform in [1e-3, 0.05]. Orientation isotropic.
# f_a in GeV (it enters as f_a / M_pl, M_pl = 1.22e19 GeV).
#
# Output: cygx1/{samples,lines,weights,rounds}_alpha<alpha>_fa<f_a>.dat and a log.
# Each realization's solver is capped at NL_MAX_WALL_SEC (default 600 s).
# Plot: python ../plot_realizations.py . alpha<alpha>_fa<f_a>   (from cygx1/)
set -e
alpha=${1:?alpha}; fa=${2:?f_a [GeV]}; Nmax=${3:-3}; rounds=${4:-1}; Nround=${5:-60}
cd "$(dirname "$0")/.."
tag="alpha${alpha}_fa${fa}"
julia run_zoom.jl --outdir cygx1 --tag "$tag" --alpha "$alpha" --f_a "$fa" \
    --Nmax "$Nmax" --rounds "$rounds" --N_round "$Nround" \
    --M_obs 21.2 --M_obs_sigma 2.2 --a_obs_min 0.85 \
    --M_ref 21.2 --M0_min 12 --M0_max 35 --a0_min 0 --a0_max 0.998 \
    --log10_age_min 6.6812 --log10_age_max 6.8808 --fedd_min 1e-3 --fedd_max 0.05 \
    --d_mean 2.22 --d_sigma 0.18 \
    > "cygx1/${tag}.log" 2>&1
