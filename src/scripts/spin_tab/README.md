# Spin dependence of the scattering rates (November 2026)

Scripts that produced `src/rate_sve/spin_table_v1_20261123.zip` and `rates_m4fix_20261123.zip`. They were run on DESY
Maxwell and contain absolute paths of that setup (`.../Axion_SR_gw_scan/sg_ext/spin_tab` is where they live when they
run, next to the directories of `src/scripts/sg_ext`); change the paths at the top of each script to use them elsewhere.

## What the table is

The 2 -> 2 rate tables are computed at black-hole spin 0.95. For a channel whose fourth quantum is absorbed by the hole
the solver multiplies the tabulated rate by r+(a)/r+(0.95), which is the alpha -> 0 limit only. The table holds the
residual factor per channel and alpha,

    rate(alpha, a) = table(alpha, 0.95) x [r+(a)/r+(0.95) for BH-labelled channels] x R(alpha, a),

for the self-interaction and the self-gravity amplitude separately, at the spins 0, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99,
0.998 and the 14 alpha of the rate tables. For the populated low-n channels R at spin 0.5 is 0.91 at alpha = 0.05 and
0.65-0.77 at alpha = 0.1-0.25 (the r+ scaling alone is up to ~50 % too high, e.g. 322 x 322 -> 844 x BH at alpha = 0.25,
a = 0.46: R = 0.67). Channels to infinity keep R = 1. For horizon channels with m4 != 0 the tabulated quantity excludes kH = w4 - m4 Omega_H, as in the rate tables.

File format (`NL_SPIN_TABLE`, read by `src/spin_table.jl`):

    # spins 0 0.3 0.5 0.7 0.8 0.9 0.95 0.99 0.998
    S1 S2 S3 S4 alpha flag  lnR_SI[9]  lnR_SG[9]          (natural logarithms)

In the solver the factor is a correction of the rate used so far,
`rate + A (R_SI - 1) + X (sqrt(R_SI R_SG) - 1) + G (R_SG - 1)` with A, X, G the SI, cross and SG terms at the current
mass; monotone cubic (PCHIP) in spin on ln R, linear in ln alpha between rows; a channel without rows has R = 1.

## Method

`spin_fast.jl` is the production calculation (`gf_radial` of `src/solve_sr_rates.jl`: relativistic cloud modes for
n <= 8 from the cache, hydrogenic above, spheroidal angular functions, Klein-Gordon normalisation, SI and SG sources,
the same source grid and rate formula) with the two homogeneous radial solutions of the fourth quantum integrated in
double precision. At equal radial step it reproduces `gf_radial` to 1e-4 and is 100-1000 times faster, so the step in
r* is capped at 0.25, and at 0.1/|w4 - m4 Omega_H| for m4 != 0 (production: 0.12/alpha, i.e. up to 2.4). Against 3,048
points of the full calculation the spin ratios agree to 0.05 % (median), 0.2 % (90 %).

Two defects of the inputs are corrected in the table:

- `wf_bad.txt` (found by `wf_check.jl`): the relativistic radial function from `solve_radial` is wrong at a = 0.998, and
  only there, for 622, 633, 722, 733, 744, 833, 844, 855, 866 (100 level-alpha entries; overlap with the a = 0.95 function
  0.68-0.94, against >= 0.999998 for every other entry; the level energies are right). The 0.998 node of the 35,186 rows
  that contain one is computed with the a = 0.99 function (`wfsub=wf_bad.txt only_sub=1`). Anything computed directly at
  a >= 0.998 with `solve_radial` is affected; the rate tables (a = 0.95) are not.
- `lev_bad.txt` (found by `level_check.py`): 38 level energies at spins 0 and 0.3 (19-18-18, 21-20-20, 22-20-20,
  23-20-20, 24-20-20) converged to a much higher overtone (n_eff 80-130). Replaced by the a = 0.5 value plus the hyperfine
  term; 1,396 rows recomputed at these spins (`only=lev_bad.txt`).

## Order of the calculation

    make_tasks.py                      tasks/: levels, channel lists (9,520 channels whose fourth quantum can be bound)
    submit_all.sh                      spin_levels.jl (energies and n <= 8 radial functions at the 9 spins; gate
                                       check_levels.py), then spin_fast.jl on all channels (node_job.sh MODE=levels|fast)
    wf_check.jl, level_check.py        the two lists above; recomputation with node_job.sh MODE=fix and only=lev_bad.txt
    build_table.py                     fast_out/ -> spin_table_v1.txt and its report

## Horizon channels with m4 != 0

57 channels of `sgx_channels.txt` (another 210 need l4 = 24 or 30 and are zero). Their stored a = 0.95 points were wrong
wherever the production step did not resolve |w4 - m4 Omega_H|: 604 of 798 grid points by more than 10 %, 223 by more
than a factor 10. `m4_tables.py` builds corrected tables from a `spin_fast.jl kg=0` run, `m4_install.py` converts them
to the Klein-Gordon normalisation and installs them (`rates_m4fix_20261123.zip`); `scripts/sg_ext/merge_allW.py` reads
the corrected flat tables. In the evolution runs these channels carry < 2e-8 of the scattering flux.
`NL_BH_M4_SPIN=1` makes them write the angular momentum of the absorbed quantum to the spin row (off by default).

## Solver switches added with this work (`src/solve_system_unified.jl`)

| Environment variable | Default | Effect |
|---|---|---|
| `NL_SPIN_TABLE=<file>` | unset | spin factors of the table; unset reproduces earlier runs to the bit |
| `NL_BH_M4_SPIN=1` | off | spin row of horizon channels with m4 != 0 |
| `NL_DT_EFOLD=c` | 1 | a step stays below ~c e-folds of any growing level; 0: no such limit (before 2026-11-23) |
| `NL_RHS_TRACE=t0,t1` | unset | prints the mass and spin rows by term for every RHS call in the window |

The step limit removes a failure of the earlier runs: with steps of ~5e6 yr a level whose growth is the small difference
of a scattering gain and a GW loss (844 or 13-4-4 against 322, equal to 1-2 %) could gain 26 instead of 8 e-folds in one
step, after which the saved rows show the hole losing mass at fixed angular momentum and sometimes the state is
destroyed (4 of 408 runs of the November sets, all at f_a = 1e17 GeV, alpha 0.2-0.25). Whether such a level grows at all
is decided at the 1-2 % level, which is below the accuracy of the rates.
