# Exact (fully relativistic) one-graviton exchange for 2 -> 2 level transitions

Branch `relativistic-sg`.  Replaces the Newtonian self-gravity amplitude of `gf_radial(...; self_grav=true)` /
`ongoing_tests/sg_vs_si` by the exact linear-in-cloud-mass result on Kerr, for N1 x N1 -> N3 x {BH, Inf}.

## What is computed

1. **Transition stress tensor** from the Kerr quasi-bound modes (Leaver series in horizon-regular ingoing form,
   `modes.jl`; they satisfy the Klein-Gordon equation on Kerr to ~1e-12):
   `T31_mn = d_m phi1 d_n phi3* + d_n phi1 d_m phi3* - g_mn (d phi1 . d phi3* + mu^2 phi1 phi3*)`,
   frequency `w1 - w3`, azimuthal number `m1 - m3`.
2. **Linearised Einstein equations on Kerr**, `dG_mn[h] = 8 pi T31_mn`, solved as a 2D (r, theta) boundary-value
   problem (`solver.jl`): Lorenz-gauge-fixed operator, ingoing Kerr coordinates (regularity at the future horizon picks
   the ingoing solution), spin-weighted spherical harmonics for the 10 components (Galerkin in theta), 8th-order finite
   differences on a log-like radial grid, exterior complex scaling for the outgoing boundary condition.  The pointwise
   operators (`jets.jl`) are built from the metric by automatic differentiation and checked against exact identities
   (`test_ops.jl`: dG vanishes on mass and pure-gauge perturbations; the scalar coupling satisfies the Lie-derivative
   identity), the solver against a manufactured solution (`test_solver.jl`).
3. **Coupling to mode 4**: `(box - mu^2) phi4 = -delta box[h] phi1 = h^mn nabla_m nabla_n phi1 + Z^n d_n phi1`, projected
   on mode 4 and pushed through the *same* Kerr Green's function as the existing SI/SG rates (`gf_amplitudes` of
   `sg_vs_si.jl`).  Normalised so that `A_GR/A_N -> 1` in the Newtonian limit (`exchange.jl`).

The output is the complex ratio `rho = A_GR/A_N` per channel and alpha, with `A_N` the Newtonian (fixed-hole)
self-gravity amplitude of the existing pipeline, computed here with the same modes.  Hence
`Gamma_SG -> |rho|^2 Gamma_SG` and `Gamma_x -> Re(rho) Gamma_x` (A_SI and A_N are in phase).

## Black-hole recoil (dipole channels)

For transitions whose density has a dipole (`|m1 - m3| = 1`: 211x211->322xBH and 322x322->211xInf) the hole is not a
fixed centre: it is pulled back and forth by the oscillating dipole of the cloud.  In the hole's frame this is the
fictitious ("indirect") field `mu A.x`, with `A = G int rho y/|y|^3`, of the same order as the direct interaction
(equivalently the mass-polarisation term `p_i.p_j/M` of the three-body problem, which for gravity is not suppressed by
mu/M because both couplings are gravitational).  The Newtonian rates used so far omit it.  In multipole form it replaces
the L = 1 kernel `x.y / r_>^3` by `x.y (1/r_>^3 - 1/|y|^3)` (`newtonian_sources`, "NR").  Effect: `A_NR/A_N ~ -0.59`
(211 BH) and `~ -0.6` (322 -> 211 Inf) at small alpha, i.e. the rate drops to ~0.35-0.4 of the fixed-hole value and the
SG amplitude changes sign relative to SI.

The exact calculation contains the recoil automatically, but numerically it is delicate: with decaying boundary
conditions the hole's displacement is a near-zero mode of the linearised operator (eigenvalue ~ w^2, w = w1 - w3), and
the corresponding gauge field near the horizon is ~alpha^-4 larger than the cloud's own field, so discretisation errors
are amplified by ~1/w^2 (up to 1e10).  `pinned.jl` therefore solves in the hole's frame: the horizon is pinned (no
l = 1 displacement, and no tilt for Kerr) with auxiliary strut forces in a bordered system, a uniform-acceleration
(and rotating) frame field is added analytically, and its amplitude is fixed by requiring the struts to cancel.  The
result is an exact strut-free solution, gauge-equivalent to the centre-of-mass one, obtained without any 1/w^2
amplification (and independent of the choice of frame field, strut and window: checked).

## Validation

* Quadrupole channels (no recoil): `A_GR/A_N = 1 + O(alpha^2)`, e.g. 544x544->322xInf: 1.00095 (alpha 0.1),
  1.0074 (0.3); 322x322->544xBH: 1.00003 (0.03), 1.0004 (0.1), 1.0017 (0.2).
* Hole held fixed (strut not cancelled): `A/A_N` = 0.99995 (alpha 0.025), 0.99966 (0.05), 0.994 (0.2).
* Dipole channels: `A_GR/A_NR -> 1` with `O(alpha^2)` corrections: 322x322->211xInf (Schwarzschild) 1.004 (0.035),
  1.0065 (0.05), 1.026 (0.1), 1.12 (0.2); 211x211->322xBH 1.001 (0.03), ~1.0 (0.05).
* Gauge condition `|Z|/|grad h| ~ 1e-12` for the quadrupole channels; imaginary parts of rho <~ 1e-3.
* Resolution, angular truncation, frame-field shape, strut, window, FD order: changes <~ 0.2%.

## Files

| file | content |
|---|---|
| `swsh.jl` | spin-weighted harmonics + theta derivatives (BigFloat Goldberg formula) |
| `kerr.jl` | ingoing Kerr metric, Christoffels by AD, tortoise coordinates, KS normal |
| `jets.jl` | pointwise linearised Einstein / Lorenz operator and scalar coupling on 60 "jets" |
| `solver.jl` | angular basis, radial grid with complex scaling, sparse assembly |
| `modes.jl` | Kerr quasi-bound modes (Leaver, continued fraction) and spheroidal harmonics |
| `exchange.jl` | T31, Lorenz solve, mode-4 source, Newtonian (+recoil) reference sources |
| `pinned.jl` | hole-frame (pinned + strut + frame field) solver for dipole channels |
| `hybrid.jl` | earlier attempt (windowed displacement gauge); kept for reference, ill-conditioned |
| `run_exchange.jl` | driver: `julia run_exchange.jl 322 211 Inf 0.1 [spin] [Lextra] [drho]` |
| `scan.jl`, `jobs_*.txt` | batch runs -> `results_*.dat` |
| `plot_results.py` | figure + table from `results_prod.dat` |
| `make_tables.py` | `tables/<chan>_LvrHc_SGrel.dat` (exact rates in the format of `rate_sve/*_LvrHc_SG.dat`) |
| `plot_draft_figs.py` | draft figures: `fig_sg_exact`, `fig_fstar_exact`, `fig_211_exact` |
| `plot_evol_panels.py` | evolution panels (alpha x f_a grid, or alphas wrapped at one f_a) with several line sets |
| `test_*.jl`, `diag_*.jl`, `explore_*.jl` | checks and diagnostics used during development |

## Results (a = 0.95, alpha grids of the rate_sve tables; `results_prod.dat`, `fig_exact_vs_newton.pdf`)

`rho = A_GR/A_N` (Gamma_SG -> |rho|^2 Gamma_SG):

| channel | alpha = 0.03 | 0.1 | 0.2 | 0.3 | max alpha |
|---|---|---|---|---|---|
| 211x211->322xBH | -0.641 | -0.649 | -0.675 | -0.723 | -0.777 (0.373) |
| 322x322->211xInf | -0.589 | -0.595 | -0.612 | -0.641 | -0.674 (0.373) |
| 322x322->544xBH | 1.00003 | 1.0004 | 1.0021 | 1.0041 | 1.042 (0.746) |
| 544x544->322xInf | 1.00003 | 1.0004 | 1.0022 | 1.0044 | 1.040 (0.746) |

* Dipole channels: Gamma_SG is 0.41-0.60 (211 BH) and 0.35-0.45 (322 -> 211 Inf) of the Newtonian fixed-hole value,
  and the SG amplitude is opposite in sign to SI, so Gamma_x -> Re(rho) Gamma_x < 0: the two interfere destructively
  and cancel near the (new) crossover.  Since f* ∝ Gamma_SG^(-1/4), the crossover moves up by ~25%.
* Relativistic corrections beyond Newton + recoil: O(alpha^2), reaching +65% (211 BH) and +135% (322 Inf) in the rate
  at the largest alpha; the "Newton + recoil" approximation degrades fastest there.
* Quadrupole channels: +0.1% at alpha = 0.1, +4% (amplitude) / +8% (rate) at alpha = 0.75.

`make_tables.py` writes `tables/<chan>_LvrHc_SGrel.dat` (same columns as the Newtonian `*_LvrHc_SG.dat`) and
`tables/<chan>_rho.dat`.

## Recoil in the main pipeline (`gf_radial`)

`solve_sr_rates.jl`: `sg_coulomb_radial(...; recoil)`, `sg_transition(...; recoil)` and `gf_radial(...; sg_recoil=false)`
add the BH-frame (indirect) term to the L = 1 multipole of the Newtonian SG potential,
`Q_1(r) -> r^-2 int_0^r y^3 f - r int_0^r f`.  `Compute_all_rates.jl --sg_recoil true` writes recoil-corrected
`*_LvrHc_SG.dat` tables (Newton + recoil).  Only channels with a dipole transition density (|m_{1,2} - m3| = 1) change;
in the 5-level network these are 211x211->322xBH, 211x322->433xBH, 211x433->544xBH, 322x322->211xInf,
322x433->211xInf, 322x544->211xInf, 322x766->211xInf, 433x544->322xInf, 433x766->322xInf, 544x766->433xInf.
`run_Nlevels.jl --dlnM_sc` sets the mass-grid spacing of the tabulated scattering rates (needed to resolve the SI-SG
cancellation at f_a = f_a^*).

## General channels (November 2026): `exchange_gen.jl`, `scan_gen.jl`

`exchange_gen.jl` extends the calculation to N1 x N2 -> N3 x {BH, Inf} with different incoming states and any n: the
source of mode 4 is the symmetrised sum of the two exchanges, the bound states come from the radial ODE (the Leaver
series cancels catastrophically for states with radial nodes), and every multipole l4 summed by the production code is
projected. Per l4 it returns the contact (SI), Newtonian (fixed hole), Newtonian + recoil and relativistic amplitudes
with a common normalisation. `scan_gen.jl <jobfile> <outdir> [newton_only]` is the worker (job lines `S1 S2 S3 alpha`,
one result file per point, claims by lock directories so that many workers share a list).

Dipole exchanges (`|m_b - m3| = 1`) use the hole-frame solver of `pinned.jl` and need the larger angular basis
`Lextra_dip = 10`. The procedure (pin the horizon, hold it with a strut, add the accelerated-frame field, release the
hole by cancelling the strut) is written out in Section 3.5 of the note `sg_source_note` (notes are not kept in the repository). Test switch
`SGX_HOLD_FIXED=1`: the struts are left in, i.e. the hole is held at rest; the relativistic amplitude then reproduces
the Newtonian one WITHOUT recoil (0.1 % at alpha = 0.05-0.1), while the normal solution tends to Newton + recoil as
alpha^2 (`src/scripts/sg_exact/frame_check`).

Two additions of 19 November 2026 (note, Sections 3.6 and 3.7):
- The pair (s-channel) term. The exchange solved here is the t/u-channel one (transition stress tensor of an incoming
  and the outgoing level). The third pairing, the stress tensor of the two incoming levels (frequency 2 mu) acting on
  the outgoing one, is local and equivalent to `1/f_a^2 -> 1/f_a^2 + 12 pi G`. It is not in any table; `load_rates.jl`
  adds it from Gamma_SI and kappa whenever `self_grav = true` (`NL_S_CHANNEL=0` leaves it out).
- The order-alpha^2 coefficients of the exchange are confirmed by an analytic 1PN calculation in the same variables
  (`src/scripts/sg_1pn`): +0.040 against +0.039 and +0.042 against +0.041 in two channels without recoil, -1.09
  against -1.06 and -1.50 against -1.48 in two with recoil (coefficients of alpha^2 of (A - A_NR)/A_N). Test switches
  for it in `exchange_gen.jl`: `SGX_NEWTON_SHIFT`, `SGX_NEWTON_WX`, `SGX_NEWTON_WY`, `SGX_NEWTON_WRX`,
  `SGX_NEWTON_RECRET`, `SGX_NEWTON_EFAC` (Newtonian term in the harmonic radius, with the 1PN weights and frequency
  factors); `SGX_SPIN` in `scan_gen.jl`. Over eleven channels the complete 1PN and the relativistic rates agree in
  eight; in three (422 x 544, 422 x 844, 322 x 844 -> 211: the level acted on by the dipole exchange lies far outside
  the transition density) they disagree at order alpha^2 in the recoil part, and which of the two is wrong is open
  (note, Section 3.7; `src/scripts/sg_1pn/README.md`).

Products (archives in `src/rate_sve`, scripts in `src/scripts/sg_exact`):
- `sg_relativistic_points_20261117.zip`: the computed points (2,488 relativistic, 24,887 with Newtonian amplitudes and
  Kerr bound states), their job lists and level energies.
- `rates_sg_relativistic_20261117.zip`: three sets of overlay tables built from them (relativistic exchange;
  relativistic exchange + Kerr bound states; Kerr bound states with the Newtonian exchange). They have the names of
  the production tables: do not unzip them into `rate_sve`, read them with `NL_SG_OVERLAY=<dir>`.
