# Extended self-gravity system (November 2026)

Scripts that produced the self-gravity (SG) rate tables in `src/rate_sve/rates_sg_ext.zip` and `rates_si_ext.zip`, and the
evolution runs that use them. They were run on DESY Maxwell and contain absolute paths of that setup
(`/data/dust/user/samwitte/Axion_SR_sg` for the code, `.../Axion_SR_gw_scan/sg_ext` for the campaign directory, which is
where these scripts live when they run); change the paths at the top of each script to use them elsewhere.

## System

215 levels (`src/rate_sve/sgx_levels.txt`): Table I of "Toward the minimal set of states in axion superradiance"
(10 Msun column), the 21 levels of the first SG campaign, and (n,4,4) for n = 5..20. Channels
(`src/rate_sve/sgx_channels.txt`): the paper's list `load_rate_input_Nmax_45.txt` restricted to these levels, plus every
rule-allowed channel (rules of `rate_sve/gen_custom_input.py`) that contains an (n,4,4) level, among the levels with
n <= 20. `states_lib.py` defines the level sets.

Since 2026-11-17 the list also holds every other rule-allowed channel among the levels with n <= 20 that contains an
(n,8,8) level (`plan_n88.py`): 17,130 more channels, 74,009 in total, on the same 215 levels. The first 56,879 lines
of `sgx_channels.txt` are the earlier list.

## Rate tables

Each channel is computed for self-interactions and self-gravity together by `src/Compute_sg_points.jl`
(a = 0.95, BH recoil, Gauss-Legendre angles, BHlmax 2) at
alpha = 0.03, 0.04, 0.05, 0.062, 0.075, 0.09, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.28, with

- `--rel_nmax 8`: relativistic (Leaver/Heun) bound states for levels with n <= 8, hydrogenic ones above. Hydrogenic low-n
  states are not accurate enough for the SG amplitude (322 322 -> 211 + infinity is off by 2x at alpha = 0.2); a fully
  relativistic point at n = 10 costs about 10 minutes and agrees with the mixed mode to 0.5%.
- `--alpha_h 0.12`: when the fourth quantum is bound, the step of its radial integration is capped at 0.12/alpha. The
  default step is not converged at low alpha (433 433 -> 10-6-6 off by 65% at alpha = 0.05, 322 322 -> 544 by 12%).
- `--cache DIR`: the eigenvalue searches (4.6 s each, three per point) and relativistic radial solves are cached on disk.

Tables (unzip in `src/rate_sve/`, after `Rates_master_Nmax18_LvrHc.zip` and `rates_sg.zip`):

- `rates_sg_ext.zip`: `<S1>_<S2>_<S3>_<BH|Inf>_LvrHc_SG.dat` with columns alpha, Gamma_SI, Gamma_SG, Gamma_x, kappa at
  f_a = M_pl. It replaces the 1,073 tables of `rates_sg.zip` (first campaign, linear alpha grid with only two or three
  nodes below 0.25 for channels whose lowest m is 3 or more, default step): their rows above alpha = 0.28 are kept, the
  rows below are the new ones. A regime-split channel has the same table under both labels.
- `rates_si_ext.zip`: `*_LvrHc_.dat` (alpha, Gamma_SI) for the channels that have no table in the master archive.

- `rates_ext_update_20261116.zip`: 499 SG and 138 SI tables that are new or changed since the two archives above (the
  channels finished after 2026-11-14). Unzip it after the two archives above.
- `rates_ext_n88_20261117.zip`: the SG and SI tables of the 17,130 (n,8,8) channels (34,260 files), same alpha grid
  and settings. Unzip it last.

New tables stop at alpha = 0.28. Channels whose rates are zero on this grid have no SG table (179, for example
(n,16,16) pairs to (n',8,8)). 128 channels between the (n,10,10) and (n,20,20) towers at n > 20 are incomplete: their
bound fourth quanta need more than 20 GB per process at low alpha, and two attempts ran out of memory. `sgx_channels.txt`
lists all 74,009 channels; one without an SG table uses its SI table only.

Ten of the (n,8,8) channels, all of the form 15-x 20-y -> 12-z, are labelled as emitting to infinity although the
energies of the computed points give a bound fourth quantum (nine at every alpha, 15-8-8 20-6-6 -> 12-8-8 for
alpha >= 0.2). In the hydrogenic spectrum the fourth quantum of these channels sits exactly at threshold
(1/15^2 + 1/20^2 = 1/12^2), so small corrections decide the regime and their rates are not reliable either way.
In the evolution runs of this campaign (M = 10 and 1e4 Msun, alpha_0 = 0.05-0.25) the (n,8,8) channels change final
spins by at most 1e-4.

## Using the tables

`run_Nlevels.jl --self_grav true --rate_input sgx_channels.txt --keep_levels <the levels> --gw_table_nmax sgx --Nmax 45`
(see `submit_sgx_runs.py` for the full command). `--keep_levels` matters: a level that is evolved but has no scattering
channel grows by superradiance unchecked (311 and 522 did so in the first 21-level campaign at f_a = 1e16).
`NL_SI_FROM_SG=1` makes `load_rates.jl` take Gamma_SI from the SG table (same computation, alpha grid and step as
Gamma_SG and kappa) inside that table's alpha range; without it the separate `*_LvrHc_.dat` table is used.
Since 19 November 2026 `--self_grav true` also includes the s-channel graviton term (`1/f_a^2 -> 1/f_a^2 + 12 pi G`,
computed in the loader from Gamma_SI and kappa); the runs of this campaign were made without it, which
`NL_S_CHANNEL=0` reproduces. It changes their final spins by at most 2e-4.
The GW table for these levels is not tracked: `julia build_gw_nr_custom.jl ../../rate_sve/sgx_levels.txt sgx`.

## Scripts

| script | purpose |
|---|---|
| `plan_rates.py` | level and channel lists, task chunks for the rate jobs (`--only-missing` for leftovers) |
| `rates_job.sh` | Slurm job: one node of single-thread workers that claim task chunks; resumable |
| `merge_rates.py` | installs the computed points as tables in `rate_sve`, writes the run list |
| `completeness.py`, `status.sh` | progress and gaps of the rate computation |
| `submit_sgx_runs.py`, `launch_overlap.sh` | evolution runs (2 masses x 6 alphas x 5 variants x 2 GW models) |
| `chain.sh`, `post.sh`, `finish_rates.sh` | merge -> runs -> summary and figures; finishing leftover channels |
| `summarize_runs.py`, `make_figs.sh`, `plot_sg_states.py` | summary table and figures of the runs |
| `build_gw_nr_custom.jl` | non-relativistic GW table for an explicit level list |
| `rank_gw_sgx.jl`, `gw_kerr_build.sh`, `gw_kerr_channels.txt` | GW channels lacking a Kerr table; building them |
