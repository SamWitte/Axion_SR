# Relativistic self-gravity tables for the extended system (November 2026)

Scripts that turn the relativistic exchange code (`src/ongoing_tests/sg_relativistic/exchange_gen.jl`, worker
`scan_gen.jl`) into rate tables, and that use those tables in evolution runs. They were run on DESY Maxwell and
contain absolute paths of that setup (`/data/dust/user/samwitte/Axion_SR_sg` for the code,
`.../Axion_SR_gw_scan/sg_ext/sg_exact` for the directory in which these scripts live when they run); change the paths
at the top of each script to use them elsewhere. No run output is kept in the repository.

## What is computed

The production self-gravity tables (`src/scripts/sg_ext`) use the Newtonian potential with the black-hole recoil, and
hydrogenic radial functions for levels above n = 8. The relativistic code differs in both: the exchange is the
linearised Einstein equation on Kerr, and every bound state is a Kerr mode. The two changes are kept apart:

| overlay set (in `rates_sg_relativistic_20261117.zip`) | exchange | bound states | builder mode |
|---|---|---|---|
| `relativistic_exchange` | relativistic | as in the production tables | `gr` |
| `kerr_states_newtonian_exchange` | Newtonian + recoil | Kerr for every n | `kerr_nr` |
| `relativistic_exchange_kerr_states` | relativistic | Kerr for every n | `kerr_gr` |

Points: every channel among the levels with peak occupation above 1e-6 in the campaign's evolution runs, at the run's
alpha (2,488 relativistic points); and every channel among the levels above 1e-20, at the run's alpha and the next lower
grid node, with the Newtonian amplitudes only (24,887 points, about a minute each). All at spin 0.95.

## Workflow

1. `select_levels.py`, `make_jobs.py <threshold> <prefix>`: job lists (`S1 S2 S3 alpha`) from the peak occupations of
   the evolution runs. `make_jobs_newton.py`, `make_jobs_newton_runs.py`: the lists of the Newtonian-only points.
2. `exact_job.sh` / `exact_job2.sh` (batch job, several workers per node), `extra_workers.sh` (more workers inside a
   running allocation), `newton_workers.sh`: run `scan_gen.jl` on the lists. A relativistic point needs 5-20 minutes
   and up to 200 GB (many multipoles for high-n levels); `reap_locks.sh` releases the claims of workers that died,
   `watch_prod.sh` and `status.py` report what is missing.
3. `energies.py` (level energies from the worker logs), `abs_norm.py` (absolute normalisation against the production
   tables), `build_tables2.py --mode gr|kerr_nr|kerr_gr --out <dir> [result dirs]`: the overlay tables. Only the
   computed nodes are exact; rows at other alphas are the production row times 1 + c alpha^2, with c from the computed
   nodes, when all computed factors of the channel lie in (0.5, 2), and are left as tabulated otherwise. Each set
   has a `SUMMARY.tsv` listing every computed point and its ratio to the production table. `affected.txt`
   (`affected_rcut.py`) lists the points that had to be recomputed with the source cut taken from the wavefunctions.
4. `submit_exact_runs.py --name <set> --tables <overlay dir> [--rate_input <list>] [--variants ...]` (or
   `launch_runs_overlap.sh` inside a running allocation): the evolution runs of the campaign with
   `NL_SG_OVERLAY=<overlay dir> NL_SI_FROM_SG=1`. `runsum.py`, `compare_sets.py`, `tower_extent.py`: final spin, final
   mass and populated levels of a set against a reference set; `check_run_files.py`: consistency of run directories.

To rebuild the tables without recomputing the points, unzip `src/rate_sve/sg_relativistic_points_20261117.zip` next to
these scripts (it holds `prod/res`, `prod_dip/res`, `prod2/res`, `prodN/res`, `energies.pkl` and the job lists) and run
step 3: `build_tables2.py --mode gr --out tablesF_gr prod prod_dip prod2` and
`--mode kerr_gr --out tablesF_kerr_gr prod prod_dip prod2 prodN` (`kerr_nr` likewise).

## Checks

- Absolute rates of the code (Newtonian + recoil amplitudes) against the production tables for channels with all
  n <= 8, where both use Kerr states: 1,738 points, 0.1-0.2 % of them differ by more than 3 %.
- Frame / black-hole recoil (`frame_check/`): with `SGX_HOLD_FIXED=1` (hole held at rest by the struts) the
  relativistic amplitude equals the Newtonian one without recoil to 0.1 %; the normal solution tends to Newton + recoil
  as alpha^2. Over the 743 points with a dipole exchange, where the recoil changes the Newtonian rate by more than 2x
  (574 points) the relativistic rate is within 30 % of Newton + recoil in 95.6 % and of the fixed-hole value in none.
- Size of the relativistic correction per point: median +0.2 % in Gamma_SG; 2.5 % of the points change by more than
  30 %, 1.0 % by more than 2x. Main channels: 211 x 211 -> 322 x BH x1.08 / 1.23 at alpha = 0.15 / 0.25,
  322 x 322 -> 211 x inf x1.12 / 1.39, 322 x 544 -> 211 x inf x0.66 / 0.30.
- Kerr instead of hydrogenic states above n = 8: 3-6 % of the points with such a level change Gamma_SI by more than
  30 % (6-9 % for Gamma_SG), next to zeros of an amplitude.
- Not reproduced: the absolute rates of 12 points absorbed by the hole with m4 != 0 (the ratio SG / SI agrees; those
  keep the tabulated Gamma_SI). Not verified: two channels at alpha = 0.15 with a relativistic rate 1e-4 of the
  Newtonian one or less (422 x 14-6-6 -> 988 x BH, 14-6-6 x 988 -> 766 x inf).
