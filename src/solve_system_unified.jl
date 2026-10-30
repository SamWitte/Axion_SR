using DelimitedFiles
using Interpolations
include(joinpath(@__DIR__, "Core", "accretion.jl"))

# ============================================================================
# RHS / JACOBIAN KERNELS (standard mode)
# ============================================================================
# The hot loops of solve_system's RHS live here, behind function barriers with
# typed arguments. Inside the solve_system closure several captured variables
# are assigned more than once (idx_lvl, m_list, bn_list in both branches of
# `if spinone`, track_alpha reassigned), so Julia boxes them and every access in
# a loop is a dynamic dispatch: the scattering loop alone cost ~420 ns per
# channel (17 ms per RHS call at Nmax 8, 42k channels). Same operations in the
# same order, so the trajectory is unchanged to the bit.
#
# xfl: occupation at or below which a level counts as at the seed floor; a
# channel with >= 2 floor participants is skipped (-Inf: never skip).

function _scatter_terms!(du, X, rate_cache, sc_rates, track_alpha::Bool, idx_lvl::Int, massI::Int,
                         rP_ratio::Float64, alph::Float64, omH::Float64, xfl::Float64)
    for (kk, (idxV, sgn, is_bh_final, base_rate, m_drag)) in enumerate(rate_cache)
        nfl = 0
        for j in 1:length(sgn)
            k = idxV[j]
            (k <= idx_lvl && k > 0 && X[k] <= xfl) && (nfl += 1)
        end
        nfl >= 2 && continue
        if track_alpha
            base_rate = sc_rates[kk]
        end
        u_term_tot = 1.0
        for j in 1:length(sgn)
            if (idxV[j] <= idx_lvl) && (idxV[j] > 0)
                u_term_tot *= X[idxV[j]]
            end
        end
        kH = (is_bh_final && m_drag != 0) ? (alph - m_drag * omH) : 1.0
        rate_val = base_rate * (is_bh_final ? rP_ratio : 1.0) * kH
        for j in 1:length(sgn)
            idx_j = idxV[j]
            if idx_j == 0; continue; end
            if idx_j == -1; idx_j = massI; end
            du[idx_j] += sgn[j] * rate_val * u_term_tot
        end
    end
    return
end

# gw_rhs! (Numerics/gw_rates.jl) with the floor skip
function _gw_terms!(du, X, c::GWCache, α, a, xfl::Float64)
    @inbounds for i in eachindex(c.n)
        c.ω[i] = gw_omega(c.n[i], c.l[i], c.m[i], α, a)
    end
    @inbounds for k in eachindex(c.ann_i)
        i, j = c.ann_i[k], c.ann_j[k]
        (X[i] <= xfl && X[j] <= xfl) && continue
        r = c.ann_C[k] * α^c.ann_p[k] * X[i] * X[j]
        du[i] -= r
        du[j] -= r        # i == j: -2r, two quanta per annihilation
    end
    @inbounds for t in eachindex(c.tr_i)
        i, j = c.tr_i[t], c.tr_j[t]
        (X[i] <= xfl && X[j] <= xfl) && continue
        δ = c.ω[i] - c.ω[j]
        δ == 0 && continue
        r = gw_transition_rate_cached(c, t, α, abs(δ)) * X[i] * X[j]
        hi, lo = δ > 0 ? (i, j) : (j, i)
        du[hi] -= r
        du[lo] += r
    end
    return du
end

# sanitize_state! for the standard mode
function _sanitize_levels!(X, spinI::Int, n::Int, e_init::Float64, bn::Vector{Float64})
    if X[spinI] > maxSpin
        X[spinI] = maxSpin
    elseif X[spinI] < 0.0
        X[spinI] = 0.0
    end
    for i in 1:n
        if X[i] < e_init
            X[i] = e_init
        end
        if X[i] > bn[i]
            X[i] = bn[i]
        end
    end
    return
end

function _sr_terms!(du, X, SR::Vector{Float64}, ml::Vector{Int}, mu::Real, n::Int, spinI::Int, massI::Int)
    du[spinI] = 0.0
    du[massI] = 0.0
    for i in 1:n
        du[i] = SR[i] .* X[i] ./ mu
        du[spinI] += -ml[i] * SR[i] .* X[i] ./ mu
        du[massI] += -SR[i] .* X[i] ./ mu
    end
    return
end

# bosenova cap / seed floor switches and the unit conversion of the level rows
function _level_units!(du, u, ln_bn::Vector{Float64}, ln_e0::Float64, thr::Float64, mu::Real, n::Int)
    for i in 1:n
        if ((abs(u[i] - ln_bn[i]) < thr) || (u[i] > ln_bn[i])) && (du[i] > 0)
            du[i] = 0.0
        elseif (abs(u[i] - ln_e0) < thr) && (du[i] < 0)
            du[i] = 0.0
        else
            du[i] *= mu ./ hbar .* 3.15e7
        end
    end
    return
end

function _zero_off!(du, off, n::Int)
    for i in 1:n
        off[i] && (du[i] = 0.0)
    end
    return
end

# Jacobian of the scattering terms in ln N: d(prod N)/d ln N_k = prod N, once per
# occurrence of k. D accumulates the raw terms (as _scatter_terms! does into du).
function _scatter_jac!(J, D, X, live, rate_cache, sc_rates, track_alpha::Bool, idx_lvl::Int, massI::Int,
                       rP_ratio::Float64, alph::Float64, omH::Float64, xfl::Float64)
    for (kk, (idxV, sgn, is_bh_final, base_rate, m_drag)) in enumerate(rate_cache)
        nfl = 0
        for j in 1:length(sgn)
            k = idxV[j]
            (k <= idx_lvl && k > 0 && X[k] <= xfl) && (nfl += 1)
        end
        nfl >= 2 && continue
        if track_alpha
            base_rate = sc_rates[kk]
        end
        P = 1.0
        for j in 1:length(sgn)
            if (idxV[j] <= idx_lvl) && (idxV[j] > 0)
                P *= X[idxV[j]]
            end
        end
        kH = (is_bh_final && m_drag != 0) ? (alph - m_drag * omH) : 1.0
        v = base_rate * (is_bh_final ? rP_ratio : 1.0) * kH * P
        for jt in 1:length(sgn)
            it = idxV[jt]
            it == 0 && continue
            it == -1 && (it = massI)
            c = sgn[jt] * v
            D[it] += c
            for jj in 1:length(sgn)
                k = idxV[jj]
                if (k <= idx_lvl) && (k > 0) && live[k]
                    J[it, k] += c
                end
            end
        end
    end
    return
end

function _gw_jac!(J, D, X, live, c::GWCache, α, a, xfl::Float64)
    for i in eachindex(c.n)
        c.ω[i] = gw_omega(c.n[i], c.l[i], c.m[i], α, a)
    end
    for k in eachindex(c.ann_i)
        i, j = c.ann_i[k], c.ann_j[k]
        (X[i] <= xfl && X[j] <= xfl) && continue
        r = c.ann_C[k] * α^c.ann_p[k] * X[i] * X[j]
        D[i] -= r; D[j] -= r
        for o in (i, j)
            if live[o]
                J[i, o] -= r; J[j, o] -= r
            end
        end
    end
    for t in eachindex(c.tr_i)
        i, j = c.tr_i[t], c.tr_j[t]
        (X[i] <= xfl && X[j] <= xfl) && continue
        δ = c.ω[i] - c.ω[j]
        δ == 0 && continue
        r = gw_transition_rate_cached(c, t, α, abs(δ)) * X[i] * X[j]
        hi, lo = δ > 0 ? (i, j) : (j, i)
        D[hi] -= r; D[lo] += r
        for o in (i, j)
            if live[o]
                J[hi, o] -= r; J[lo, o] += r
            end
        end
    end
    return
end


"""
    solve_system(mu, fa_or_nothing, aBH, M_BH, t_max; spinone=false, ...)

Unified ODE solver for axion-BH superradiance evolution.

Handles both standard multi-level mode and spinone single-level mode via the spinone parameter.

# Arguments - Standard Mode (spinone=false)
- `mu::Float64`: Axion mass (eV)
- `fa::Float64`: Axion decay constant (1/GeV)
- `aBH::Float64`: Black hole spin (0 ≤ a ≤ 1)
- `M_BH::Float64`: Black hole mass (solar masses)
- `t_max::Float64`: Maximum integration time (years)
- `n_times::Int`: Number of output time steps (default 10000)
- `impose_low_cut::Float64`: Minimum coupling parameter threshold (default 0.01)
- `stop_on_a::Float64`: Termination spin threshold (default 0)
- `abstol::Float64`: Absolute tolerance for ODE solver (default 1e-30)
- `non_rel::Bool`: Use non-relativistic approximation (default true)
- `self_grav::Bool`: add the Newtonian self-gravity contribution (and its interference with self-interactions)
  to the scattering rates, from rate_sve/*_LvrHc_SG.dat (Compute_all_rates.jl --self_grav true); requires non_rel=false
- `si_sign::Float64`: sign of the quartic self-interaction for the SI-SG interference (+1 attractive/axion, -1 repulsive)
- `rate_input`: custom scattering-channel list (file in rate_sve/ or absolute path) instead of
  load_rate_input_Nmax_X.txt, e.g. from rate_sve/gen_custom_input.py; its states must be Nmax modes
- `gw_si_levels_only::Bool`: GW annihilations/transitions only between levels that appear in the
  scattering channel list (scattering_levels; the hard-coded rate keys for non_rel=true), so a level
  without self-interaction channels is not populated by GW either (default true). The Nmax level list
  and GW table contain more levels (e.g. the (19,18,18) truncation mode at Nmax = 18); false keeps
  them all, as before.
- `high_p::Bool`: Use high-precision tolerances (default true)
- `Nmax::Int`: Maximum principal quantum number 3-8 (default 3)
- `cheby::Bool`: Use Chebyshev interpolation (default true)
- `f_edd::Float64`: Eddington ratio of constant thin-disc accretion (default 0 = off);
  Mdot = f_edd * Mdot_Edd(M_BH) is fixed in time (Core/accretion.jl)
- `acc_eta::Float64`: radiative efficiency defining Mdot_Edd = L_Edd/(eta c^2) (default 0.1)
- `acc_dlnM_dt::Float64`: alternative to f_edd: accretion with Mdot = acc_dlnM_dt * M(t)
  [1/yr], i.e. a fixed Eddington ratio lam with acc_dlnM_dt = lam (1-eps)/(eps t_Edd)
  (default 0 = off; use either this or f_edd)

- `track_alpha`: let the SR rates (tabulated on a grid in M) and the scattering
  rates (interpolated in log alpha between grid masses) follow M(t). Default
  (`nothing`) = on (off in spinone mode); with it off they stay at the initial
  alpha, as in the original code. `track_alpha_dlnM` = (SR, scattering) grid
  spacing in ln M (default (0.0025, 0.05), converged in tests; 0.01 for SR shifted
  late spins by ~5e-3)); `track_alpha_refresh` = change in ln M that triggers a
  refresh of the scattering rates (default 0 = whenever M changes, i.e. every RHS
  call; one exp per channel). A finite value holds the rates fixed in between, and
  each refresh is a jump in the RHS: the levels relax to the new rates through a
  fast transient that collapses the step size, a sawtooth in ln N of ~1e-4-1e-3
  (BH 10, alpha 1, f_a 1e14, Nmax 8, refresh 1e-5), and the jump times depend on
  the solver's trial states, so the small modes are not reproducible. The GW rates always follow M(t);
  the bosenova caps and Emax2 always use the initial alpha. With track_alpha the
  cloud's back-reaction on M and a also uses the occupation normalisation
  u = N/(G M_BH^2) correctly (without it the old approximation M ~ M_BH is kept:
  the SR / scattering terms in dM/dt and da/dt are then too large by (M/M_BH)^2).
- `emax2_mode`: cap on the 211 cloud (211 superradiance is switched off above it).
  :fixed = Emax2, the saturated 211 cloud of a BH that starts at (M_BH, aBH) with
  no cloud (original behaviour). :evolving = current 211 cloud + the 211 cloud the
  BH can still produce from its current (M, a), i.e. Emax2 evaluated on the
  current BH, recomputed every step. It equals the fixed cap at a fresh start,
  stays constant (to first order) while only superradiance acts, grows when
  accretion spins the BH up, and depends only on the current state, so a run
  resumed from any intermediate state (u0_override) gets the same cap. Default:
  :evolving with accretion, :fixed otherwise. The evolving cap acts smoothly: the
  211 SR rate is scaled by clamp((cap - u_211) / (emax2_taper * cap), 0, 1).
  Under steady accretion the cloud sits at the cap (the spin is pinned just above
  the 211 threshold), and a hard on/off switch there makes the ODE crawl.
- `emax2_taper`: width of that taper, as a fraction of the cap (default 0.01).
- `save_factor`, `save_dlna`, `save_rel_floor`: save-on-change output. Besides the
  log-spaced saveat grid, an accepted step is saved whenever, since the last
  saved point, a relevant level's occupation changed by more than a factor
  save_factor (default 5) or the spin changed by |d ln a| > save_dlna (default
  1e-3). Relevant = within a factor save_rel_floor (default 1e-10) of the largest
  occupation, so levels growing up from the seed floor do not trigger saves. The
  solver already takes small steps where the state changes quickly; this only
  records them, and when a single step changes the state by more than the
  thresholds, extra points inside it are taken from the step's interpolant (as
  for saveat), at most save_max_sub per step. No extra RHS evaluations; the
  solution is unchanged. Set save_factor = save_dlna = Inf to switch it off.
- `skip_floor_channels`: skip scattering and GW channels with two or more
  participating levels at the seed floor (u <= e_init, within the solver's floor
  band). Such a channel adds at most ~e_init (relative) to any d ln N/dt, below
  double precision next to the other terms; identical trajectories in tests
  (Nmax 8 and 15), and most channels at large Nmax (default true).
- `solver`: :trbdf2 (default) or :radau (RadauIIA5, 5th order; far fewer steps at
  these tolerances, best with `analytic_jac`). Agreed with :trbdf2 to ~1e-8 in the
  final spin at Nmax 15 and 18; the final BH mass differed by 5e-8 (Nmax 15,
  alpha 0.2, f_a 1e14) and 3.2e-7 (Nmax 18, alpha 1.3, f_a 1e16, with either
  Jacobian), not yet understood.
  With :radau the saved points between step ends are bounded, level by level,
  by the values at the two ends of their step (see cback_saveat): the stage
  values of stiff levels are not accurate, and the raw interpolant overflowed.
- `analytic_jac`: exact Jacobian for the level columns (one pass over the
  channels; in ln N every term is a rate times a product of N's), finite
  differences for the spin and mass columns (2 RHS calls), instead of finite
  differences in every column (one RHS call per state). Matches finite
  differences to ~1e-5 (their own error) along an Nmax 8 run to 5e6 yr. A
  Jacobian error would only slow Newton; the solution is controlled by the error
  estimate on the RHS. Default false.

# Arguments - Spinone Mode (spinone=true)
- `mu::Float64`: Axion mass (eV)
- `fa_or_nothing::Any`: Ignored in spinone mode (for signature compatibility)
- `aBH::Float64`: Black hole spin
- `M_BH::Float64`: Black hole mass
- `t_max::Float64`: Maximum integration time
- `n_times::Int`: Number of output time steps

# Returns
- `(spinBH::Float64, MassB::Float64)`: Final spin and mass after evolution

# Details
Both modes integrate the coupled ODE system tracking:
- Energy populations E_nlm for axion cloud states
- Black hole spin a
- Black hole mass M

Standard mode: Multi-quantum level with scattering terms, bosenova boundaries, Emax2 cutoff
Spinone mode: Single quantum level with precomputed rates, simpler spin dynamics
"""
function solve_system(mu, fa_or_nothing, aBH, M_BH, t_max;
    n_times=10000, debug=false, impose_low_cut=0.01, return_all_info=false,
    eq_threshold=1e-100, stop_on_a=0, abstol=1e-30, non_rel=true, high_p=true,
    N_pts_interp=200, N_pts_interpL=200, Nmax=3, cheby=true, spinone=false, lm_only=false,
    u0_override=nothing, t_start=0.0,
    gw_model=:nonrel, gw_min_rate_per_yr=1e-10, gw_literature=true, gw_si_levels_only=true,
    f_edd=0.0, acc_eta=0.1, acc_dlnM_dt=0.0, track_alpha=nothing, emax2_mode=nothing,
    emax2_taper=0.01, save_factor=5.0, save_dlna=1e-3, save_rel_floor=1e-10, save_max_sub=100,
    track_alpha_dlnM=(0.0025, 0.05), track_alpha_refresh=0.0, self_grav=false, si_sign=1.0, rate_input=nothing,
    skip_floor_channels=true, solver=:trbdf2, analytic_jac=false)

    # ============================================================================
    # PARAMETER SETUP & VALIDATION
    # ============================================================================
    t_setup0 = time()
    alph = GNew .* M_BH .* mu
    (f_edd > 0 && acc_dlnM_dt > 0) && error("give either f_edd or acc_dlnM_dt, not both")
    accreting = f_edd > 0 || acc_dlnM_dt > 0
    if accreting && spinone
        error("accretion is not implemented for spinone mode")
    end
    Mdot_acc = f_edd > 0 ? f_edd * eddington_rate(M_BH; eta=acc_eta) : 0.0  # M_sun/yr, constant part
    acc_mdot(M) = Mdot_acc + acc_dlnM_dt * M                                 # M_sun/yr
    track_alpha = (track_alpha === nothing) ? !spinone : track_alpha
    emax2_mode = (emax2_mode === nothing) ? (accreting ? :evolving : :fixed) : emax2_mode
    emax2_mode in (:fixed, :evolving) || error("emax2_mode must be :fixed or :evolving")
    if track_alpha && spinone
        error("track_alpha is not implemented for spinone mode")
    end

    # Compute tolerances based on physical regime
    default_reltol, reltol_Thres = initialize_solver_tolerances(non_rel, high_p)

    # Override for testing (lines 93-95 in original)
    # default_reltol = 1e-7


    # ============================================================================
    # QUANTUM LEVEL SETUP (spinone vs standard mode)
    # ============================================================================
    if spinone
        # SPINONE MODE: Single quantum level
        idx_lvl, m_list, bn_list, modes = setup_quantum_levels_spinone()
        fa = nothing  # Not used in spinone mode
        Emax2 = nothing
        Mvars_keys = [:mu, :aBH, :M_BH]
    else
        # STANDARD MODE: Multiple quantum levels
        fa = fa_or_nothing  # Unpack the fa parameter
        idx_lvl, m_list, bn_list, modes = setup_quantum_levels_standard(Nmax, fa, M_pl, alph, aBH)

        # Compute Emax2 cutoff for 211 level
        Emax2 = 1.0
        OmegaH = aBH ./ (2 .* (GNew .* M_BH) .* (1 .+ sqrt.(1 .- aBH.^2)))
        if (OmegaH .> ergL(2, 1, 1, mu, M_BH, aBH))
            Emax2 = emax_211(M_BH, mu, aBH)
        end
        Mvars_keys = [:mu, :fa, :Emax2, :aBH, :M_BH, :impose_low_cut]
    end

    # ============================================================================
    # STATE VECTOR INITIALIZATION
    # ============================================================================
    e_init = 1.0 ./ (GNew .* M_BH.^2 .* M_to_eV)  # unitless
    spinI = idx_lvl + 1
    massI = spinI + 1

    y0, reltol = setup_state_vectors(idx_lvl, aBH, M_BH, e_init, default_reltol)
    # Resume: replace log-space initial state; rates still use original aBH/M_BH above/below.
    if u0_override !== nothing
        if length(u0_override) != length(y0)
            error("u0_override length $(length(u0_override)) != expected $(length(y0)) (idx_lvl=$idx_lvl)")
        end
        y0 = copy(u0_override)
        println("Resuming with overridden u0 at t_start=", t_start,
                "  spin=", exp(y0[spinI]), "  mass=", exp(y0[massI]))
        flush(stdout)
    end

    # ============================================================================
    # RATE SETUP
    # ============================================================================
    if spinone
        # Spinone: Precomputed rates
        wR, wI = precomputed_spin1(alph, aBH, M_BH)
        SR_rates = [2 .* wI]
        if SR_rates[1] < 1e-100
            SR_rates[1] = 1e-100
        end
        Mvars = [mu, aBH, M_BH]
        rates = Dict()  # Not used in spinone mode
        interp_funcs = Function[]
    else
        # Standard: Compute interpolated rates using smooth symlog interpolation
        SR_rates, interp_funcs, interp_dict = compute_sr_rates_smooth(modes, M_BH, aBH, alph, cheby=cheby)
        rates = load_rate_coeffs(mu, M_BH, aBH, fa, Nmax, SR_rates; non_rel=non_rel, lm_only=lm_only, self_grav=self_grav, si_sign=si_sign, rate_input=rate_input)
        Mvars = [mu, fa, Emax2, aBH, M_BH, impose_low_cut]
        rP_initial = 1.0 + sqrt(1.0 - aBH^2)
    end

    # ============================================================================
    # TRACK_ALPHA: rate grids in BH mass so that alpha can follow M(t)
    # ============================================================================
    # M_lo allows for SR spin-down mass loss (up to ~30% at large alpha); M_hi bounds the mass reached by
    # t_max. SR rates: log-spaced nodes (1%), linear in ln M between them.
    # Scattering rates: coarser nodes (5%), linear in ln(rate) vs ln M.
    M_hi = track_alpha ? 1.02 * (M_BH + Mdot_acc * t_max) * exp(acc_dlnM_dt * t_max) : M_BH
    lnM_lo, nM_sr, dlnM_sr, M_sr, sr_nodes, nM_sc, dlnM_sc, rates_sc = if track_alpha
        M_lo = 0.7 * M_BH
        n_sr = max(2, ceil(Int, log(M_hi / M_lo) / track_alpha_dlnM[1]) + 1)
        d_sr = log(M_hi / M_lo) / (n_sr - 1)
        Ms_sr = exp.(log(M_lo) .+ d_sr .* (0:(n_sr - 1)))
        n_sc = max(2, ceil(Int, log(M_hi / M_lo) / track_alpha_dlnM[2]) + 1)
        d_sc = log(M_hi / M_lo) / (n_sc - 1)
        Ms_sc = exp.(log(M_lo) .+ d_sc .* (0:(n_sc - 1)))
        if debug
            println("track_alpha: f_edd=$(f_edd), Mdot=$(Mdot_acc) M_sun/yr, dlnM/dt=$(acc_dlnM_dt)/yr; SR grid $(n_sr) masses, ",
                    "scattering grid $(n_sc) masses in [$(M_lo), $(M_hi)] M_sun")
        end
        (log(M_lo), n_sr, d_sr, collect(Ms_sr),
         [compute_sr_rates_smooth(modes, Mj, aBH, GNew * Mj * mu, cheby=cheby)[2] for Mj in Ms_sr],
         n_sc, d_sc,
         [load_rate_coeffs(mu, Mj, aBH, fa, Nmax, SR_rates; non_rel=non_rel, lm_only=lm_only, self_grav=self_grav, si_sign=si_sign, rate_input=rate_input) for Mj in Ms_sc])
    else
        (0.0, 2, 1.0, Float64[], Vector{Any}[], 2, 1.0, Dict[])
    end
    # node index and weight for linear interpolation in ln M
    function mass_bracket(M, nM, dlnM)
        x = (log(M) - lnM_lo) / dlnM
        x = isnan(x) ? 0.0 : clamp(x, 0.0, Float64(nM - 1))   # trial Newton states may be wild
        j = clamp(floor(Int, x) + 1, 1, nM - 1)
        return j, clamp(x - (j - 1), 0.0, 1.0)
    end

    # GW annihilations/transitions for every pair of levels that take part in the
    # scattering channels (gw_si_levels_only); gw_rhs! re-evaluates them each RHS
    # call at the current BH mass and spin (single assignment keeps the RHS
    # closure type-stable).
    gw_levels = if spinone || !gw_si_levels_only
        nothing
    elseif non_rel
        Set(lv for k in keys(rates) for lv in rate_key_levels(k))
    else
        scattering_levels(Nmax; lm_only=lm_only, rate_input=rate_input)
    end
    gw_cache = spinone ? gw_empty_cache(Tuple{Int,Int,Int}[]) :
               gw_build_cache(Nmax, modes, mu, M_BH, aBH; gw_model=gw_model,
                              min_rate_per_yr=gw_min_rate_per_yr, literature_overrides=gw_literature,
                              M_cut=M_hi, levels=gw_levels)
    if debug && !spinone
        lv_note = gw_levels === nothing ? "all $(length(modes)) levels" :
                  "$(count(md -> (md[1], md[2], md[3]) in gw_levels, modes)) of $(length(modes)) levels (those in scattering channels)"
        println("GW ($(gw_model)): $(length(gw_cache.ann_i)) annihilation + $(length(gw_cache.tr_i)) transition channels, $(lv_note)")
    end

    # ============================================================================
    # PRE-COMPUTE RATE INDEX CACHE (Fix 1: avoid O(N_modes × N_rates) per RHS call)
    # ============================================================================
    # key_to_indx is O(N_modes) with string allocations per call. For Nmax=15 with
    # ~22k rate keys × 560 modes this was ~18 min/Jacobian. Cache it once here.
    #
    # m_drag is the net azimuthal quantum number absorbed by the horizon for BH-type
    # keys, used below to weight the rate by the horizon superradiance factor. Only
    # computed for is_bh keys: non-BH keys (e.g. "*_*^Inf") don't have a
    # state label in the trailing slot, so get_m would fail to parse it there.
    get_m(s) = occursin("-", s) ? parse(Int, split(s, "-")[3]) : parse(Int, s[end:end])
    rate_cache = if !spinone
        map(collect(keys(rates))) do k
            if abs.(rates[k]) .> 1e20
                rates[k] = 0.0
            end
            idxV, sgn = key_to_indx(k, Nmax)
            is_bh = any(idxV .== -1)
            m_drag = 0
            if is_bh
                parts = split(k, "_", limit=2)
                caret = split(parts[2], "^")
                if length(caret) == 2
                    # "state1_state2^BH": both quanta absorbed by the horizon
                    m_drag = get_m(String(parts[1])) + get_m(String(caret[1]))
                elseif length(caret) == 3
                    # "state1_state2^state3^BH": net m into the horizon = m1 + m2 - m3
                    m_drag = get_m(String(parts[1])) + get_m(String(caret[1])) - get_m(String(caret[2]))
                end
            end
            (idxV, sgn, is_bh, rates[k], m_drag)
        end
    else
        Vector{Tuple{Vector{Int}, Vector{Int}, Bool, Float64, Int}}()
    end
    # Accretion: ln|rate| of every rate_cache entry at each scattering-grid mass
    # (same key order as rate_cache; sign taken from the initial-mass rate).
    rate_keys = spinone ? String[] : collect(keys(rates))
    lnR_sc = if track_alpha
        [log(max(abs(get(rates_sc[j], k, 0.0)) < 1e20 ? abs(get(rates_sc[j], k, 0.0)) : 0.0, 1e-300))
         for k in rate_keys, j in 1:nM_sc]
    else
        zeros(0, 0)
    end
    # Current scattering rates under track_alpha, recomputed only when ln M has
    # moved by more than 1e-5 since the last refresh (rates ~ alpha^p with p <~ 11,
    # so the error is <~ 1e-4), instead of one exp per rate in every RHS call.
    sc_rates_now = [r[4] for r in rate_cache]
    sc_lnM = Ref(-Inf)
    function refresh_sc_rates!(M)
        lnM = log(M)
        abs(lnM - sc_lnM[]) <= track_alpha_refresh && lnM == lnM && return
        j, w = mass_bracket(M, nM_sc, dlnM_sc)
        for kk in eachindex(sc_rates_now)
            b = rate_cache[kk][4]
            sc_rates_now[kk] = b == 0.0 ? 0.0 : sign(b) * exp((1 - w) * lnR_sc[kk, j] + w * lnR_sc[kk, j + 1])
        end
        sc_lnM[] = lnM
        return
    end

    # ============================================================================
    # ODE SETUP
    # ============================================================================
    if t_start >= t_max
        error("t_start ($t_start) >= t_max ($t_max); nothing to evolve")
    end
    tspan = (Float64(t_start), Float64(t_max))
    # Log-spaced save points so that early-time dynamics are resolved on the log x-axis.
    # Full-grid base uses (0, t_max) so resume hops splice onto the same grid.
    t_log_start = max(1.0, t_max * 1e-9)
    saveat_full = exp10.(range(log10(t_log_start), log10(t_max), length=n_times))
    saveat_full[end] = t_max      # exp10(log10(t_max)) can round above t_max and never be saved
    if t_start > 0.0
        saveat = saveat_full[saveat_full .> t_start]
        if isempty(saveat)
            saveat = [t_max]
        end
        println("Resume saveat: $(length(saveat)) / $(length(saveat_full)) points after t_start")
        flush(stdout)
    else
        saveat = saveat_full
    end

    # Spin callback: without accretion the spin cannot exceed its initial value
    # (reset to aBH beyond aBH + 0.01); with accretion it may spin up to maxSpin.
    spin_reset = accreting ? maxSpin : aBH
    spin_ceiling = accreting ? maxSpin : aBH + 0.01

    # Trackers for callbacks
    wait = 0
    turn_off = fill(false, idx_lvl)
    turn_off_M = false

    # Typed, singly-assigned copies for the RHS kernels (see the note above _scatter_terms!)
    n_lvl = Int(idx_lvl)
    m_vec = spinone ? Int[] : Int[m for m in m_list]
    bn_vec = spinone ? Float64[] : Float64[b for b in bn_list]
    ln_bn_vec = log.(bn_vec)
    ln_e0 = log(e_init)
    track_alpha_on = track_alpha === true
    x_floor = skip_floor_channels ? e_init * exp(SOLVER_TOLERANCES.bosenova_threshold) : -Inf

    # ============================================================================
    # SHARED HELPERS: used by both RHS and check_timescale to keep logic identical
    # ============================================================================
    function sanitize_state!(u_real)
        if u_real[spinI] > maxSpin
            u_real[spinI] = maxSpin
        elseif u_real[spinI] < 0.0
            u_real[spinI] = 0.0
        end

        for i in 1:idx_lvl
            if u_real[i] < e_init
                u_real[i] = e_init
            end
            if !spinone && u_real[i] > bn_list[i]
                u_real[i] = bn_list[i]
            end
        end
    end

    # 211 cloud the current BH (M, a) can still produce before saturating, in
    # u units (N / (G M_BH^2)); emax_211 is normalised to G M^2.
    function emax2_remaining(u_real)
        M, a = u_real[massI], clamp(u_real[spinI], 0.0, maxSpin)
        OmH = a / (2 * GNew * M * (1 + sqrt(1 - a^2)))
        OmH > ergL(2, 1, 1, mu, M, a) || return 0.0
        return max(emax_211(M, mu, a), 0.0) * (M / M_BH)^2
    end

    function compute_SR_rates_local(u_real)
        if spinone
            OmegaH = u_real[spinI] ./ (2 .* (GNew .* u_real[massI]) .* (1 .+ sqrt.(1 .- u_real[spinI].^2)))
            wR, wI = precomputed_spin1(alph, u_real[spinI], u_real[massI])
            if wR .> OmegaH
                return [0.0], true
            end
            SR_rates_local = [2 .* wI]
            if u_real[1] .>= u_real[spinI]
                SR_rates_local[1] = 0.0
            end
            return SR_rates_local, false
        else
            spin_val = clamp(u_real[spinI], 0.0, maxSpin)
            if track_alpha
                # rate [eV] at node mass M_j is (dimensionless) * 2/(G M_j)
                M_now = u_real[massI]
                j, w = mass_bracket(M_now, nM_sr, dlnM_sr)
                fl, fr = sr_nodes[j], sr_nodes[j + 1]
                SR_rates_local = [((1 - w) * fl[i](spin_val) * M_sr[j] + w * fr[i](spin_val) * M_sr[j + 1]) / M_now
                                  for i in 1:idx_lvl]
            else
                SR_rates_local = [func(spin_val) for func in interp_funcs]
            end
            if emax2_mode == :evolving
                if SR_rates_local[1] > 0
                    rem = emax2_remaining(u_real)
                    cap = u_real[1] + rem
                    SR_rates_local[1] *= clamp(rem / (emax2_taper * cap), 0.0, 1.0)
                end
            elseif (u_real[1] .> Emax2) && (SR_rates_local[1] > 0)
                SR_rates_local[1] = 0.0
            end
            return SR_rates_local, false
        end
    end

    # ============================================================================
    # RHS FUNCTION
    # ============================================================================
    function RHS_ax!(du, u, Mvars, t)
        u_real = exp.(u)
        if spinone
            sanitize_state!(u_real)
        else
            _sanitize_levels!(u_real, spinI, n_lvl, e_init, bn_vec)
        end

        SR_rates_local, should_zero = compute_SR_rates_local(u_real)

        if spinone && should_zero
            du .*= 0.0
            return
        end

        # Superradiance terms
        if spinone
            du[spinI] = 0.0
            du[massI] = 0.0
            for i in 1:idx_lvl
                du[i] = SR_rates_local[i] .* u_real[i] ./ mu
                du[spinI] += -m_list[i] * SR_rates_local[i] .* u_real[i] ./ mu
                du[massI] += -SR_rates_local[i] .* u_real[i] ./ mu
            end
        else
            _sr_terms!(du, u_real, SR_rates_local isa Vector{Float64} ? SR_rates_local : Float64.(SR_rates_local),
                       m_vec, mu, n_lvl, spinI, massI)
        end

        # Scattering and GW terms (standard mode only)
        if !spinone
            rP_now = 1.0 + sqrt(1.0 - u_real[spinI]^2)
            rP_ratio_now = rP_now / rP_initial
            a_now = u_real[spinI]
            omH_now = a_now / (rP_now^2 + a_now^2)
            alph_now = GNew * u_real[massI] * mu
            track_alpha && refresh_sc_rates!(u_real[massI])
            _scatter_terms!(du, u_real, rate_cache, sc_rates_now, track_alpha_on, n_lvl, massI,
                            rP_ratio_now, alph_now, omH_now, x_floor)
            _gw_terms!(du, u_real, gw_cache, alph_now, a_now, x_floor)
        end

        # Unit corrections
        if spinone
            for i in 1:idx_lvl
                if u_real[i] < e_init
                    du[i] = 0.0
                else
                    du[i] *= mu ./ hbar .* 3.15e7
                end
            end
        else
            _level_units!(du, u, ln_bn_vec, ln_e0, SOLVER_TOLERANCES.bosenova_threshold, mu, n_lvl)
        end

        du[spinI] *= mu ./ hbar .* 3.15e7
        if track_alpha
            # u = N / (G M_BH^2): dJ = -m dN -> da = -m du (M_BH/M)^2; dM = -mu dN = -mu G M_BH^2 du
            du[spinI] *= (M_BH / u_real[massI])^2
            du[massI] *= (mu .* M_BH) .* (mu .* GNew .* M_BH) ./ hbar .* 3.15e7
        else
            du[massI] *= (mu .* u_real[massI]) .* (mu .* GNew .* u_real[massI]) ./ hbar .* 3.15e7
        end

        if accreting   # already per year
            dM_acc, da_acc = accretion_rhs(u_real[massI], u_real[spinI], acc_mdot(u_real[massI]))
            du[massI] += dM_acc
            du[spinI] += da_acc
        end

        du ./= u_real

        _zero_off!(du, turn_off, n_lvl)
        return
    end

    # Analytic Jacobian (analytic_jac): level columns exact in ln N, spin and mass
    # columns by forward differences (SR interpolants, rate grids in M). Mirrors
    # RHS_ax!: D holds the raw terms before the unit factors.
    jac_f0 = zeros(massI); jac_f1 = zeros(massI); jac_up = zeros(massI)
    function jac_ax!(J, u, Mvars, t)
        n = length(u)
        fill!(J, 0.0)
        X0 = exp.(u)
        X = copy(X0)
        _sanitize_levels!(X, spinI, n_lvl, e_init, bn_vec)
        live = [X0[i] == X[i] for i in 1:n]      # clamped components do not move
        SR, _ = compute_SR_rates_local(X)
        D = zeros(n)
        for i in 1:n_lvl
            g = SR[i] * X[i] / mu
            D[i] += g; D[spinI] += -m_vec[i] * g; D[massI] += -g
            if live[i]
                J[i, i] += g; J[spinI, i] += -m_vec[i] * g; J[massI, i] += -g
            end
        end
        a_now = X[spinI]
        rP_now = 1.0 + sqrt(1.0 - a_now^2)
        omH_now = a_now / (rP_now^2 + a_now^2)
        alph_now = GNew * X[massI] * mu
        track_alpha && refresh_sc_rates!(X[massI])
        _scatter_jac!(J, D, X, live, rate_cache, sc_rates_now, track_alpha_on, n_lvl, massI,
                      rP_now / rP_initial, alph_now, omH_now, x_floor)
        _gw_jac!(J, D, X, live, gw_cache, alph_now, a_now, x_floor)
        K = mu ./ hbar .* 3.15e7
        thr = SOLVER_TOLERANCES.bosenova_threshold
        for i in 1:n_lvl
            if ((abs(u[i] - ln_bn_vec[i]) < thr) || (u[i] > ln_bn_vec[i])) && (D[i] > 0)
                D[i] = 0.0; J[i, :] .= 0.0
            elseif (abs(u[i] - ln_e0) < thr) && (D[i] < 0)
                D[i] = 0.0; J[i, :] .= 0.0
            else
                D[i] *= K; @views J[i, 1:n_lvl] .*= K
            end
        end
        fs = track_alpha ? K * (M_BH / X[massI])^2 : K
        fm = track_alpha ? (mu .* M_BH) .* (mu .* GNew .* M_BH) ./ hbar .* 3.15e7 :
                           (mu .* X[massI]) .* (mu .* GNew .* X[massI]) ./ hbar .* 3.15e7
        D[spinI] *= fs; @views J[spinI, 1:n_lvl] .*= fs
        D[massI] *= fm; @views J[massI, 1:n_lvl] .*= fm
        for i in 1:n                              # f = D / X
            @views J[i, 1:n_lvl] ./= X[i]
            (i <= n_lvl && live[i]) && (J[i, i] -= D[i] / X[i])
        end
        for i in 1:n_lvl
            turn_off[i] && (@views J[i, :] .= 0.0)
        end
        RHS_ax!(jac_f0, u, Mvars, t)
        jac_up .= u
        for j in (spinI, massI)
            h = 1e-7 * max(1.0, abs(u[j]))
            jac_up[j] = u[j] + h
            RHS_ax!(jac_f1, jac_up, Mvars, t)
            jac_up[j] = u[j]
            @views J[:, j] .= (jac_f1 .- jac_f0) ./ h
        end
        return
    end

    # ============================================================================
    # CALLBACKS (spinone vs standard mode)
    # ============================================================================

    # Shared: check_spin and affect_spin! (with spinone-specific differences)
    function check_spin(u, t, integrator)
        wait += 1
        u_real = exp.(u)

        if spinone
            # Spinone: No stop_on_a check
            if u_real[spinI] .> (aBH .+ 0.01)
                return true
            elseif u_real[spinI] .<= 0.0
                return true
            else
                return false
            end
        else
            # Standard: Full checks including stop_on_a
            if u_real[spinI] <= stop_on_a
                return true
            end
            if u_real[spinI] .> spin_ceiling
                return true
            elseif u_real[spinI] .<= 0.0
                return true
            else
                return false
            end
        end
    end

    function affect_spin!(integrator)
        u_real = exp.(integrator.u)
        if !spinone && u_real[spinI] <= stop_on_a
            terminate!(integrator)
        end
        if u_real[spinI] .> spin_reset
            integrator.u[spinI] = log(spin_reset)
        elseif u_real[spinI] .< 0.0
            integrator.u[spinI] = -10.0
        end
        set_proposed_dt!(integrator, integrator.dt .* 0.3)
    end

    # Standard mode only: check_timescale and affect_timescale!
    function check_timescale(u, t, integrator)
        u_real = exp.(u)
        sanitize_state!(u_real)
        SR_rates_local, _ = compute_SR_rates_local(u_real)

        for i in 1:idx_lvl
            if (u[i] < log(1e-75)) && (SR_rates_local[i] < 0)
                turn_off[i] = true
            elseif turn_off[i] && (SR_rates_local[i] > 0)
                turn_off[i] = false
            end
        end

        if u_real[massI] > (1.4 * M_BH)
            turn_off_M = true
        end

        integrator.opts.reltol = reltol
        
        du = get_du(integrator)

        tlist = Float64[]
        for i in 1:idx_lvl
            condBN = (abs(u[i] - log(bn_list[i])) < SOLVER_TOLERANCES.bosenova_threshold)
            if (u[i] > log(e_init)) && condBN && (du[i] != 0.0)
                push!(tlist, abs(1.0 ./ du[i]))
            end
        end

        if du[spinI] != 0.0
            push!(tlist, abs(def_spin_tol ./ du[spinI]))
        end

        if isempty(tlist); return false; end
        tmin = minimum(tlist)

        if (integrator.dt ./ tmin .>= 0.1); return true
        elseif (integrator.dt ./ tmin .<= 0.001); return true
        elseif (integrator.dt .<= 1e-12); return true
        else; return false
        end
    end

    function affect_timescale!(integrator)
        du = get_du(integrator)
        tlist = Float64[]
        indx_list = Int[]
        for i in 1:idx_lvl
            condBN = (abs(integrator.u[i] - log(bn_list[i])) < SOLVER_TOLERANCES.bosenova_threshold)
            if (integrator.u[i] > log(e_init)) && condBN && (du[i] != 0.0)
                push!(tlist, 1.0 ./ du[i])
                push!(indx_list, i)
            end
        end

        if du[spinI] != 0.0
            push!(tlist, def_spin_tol ./ du[spinI])
        end

        if isempty(tlist); return; end
        tmin = minimum(abs.(tlist))

        if (integrator.dt ./ integrator.t < 1e-6) && (wait % 1000 == 0) && (wait > 10000)
            for i in 1:idx_lvl
                if reltol[i] < reltol_Thres
                    reltol[i] *= 1.2
                    integrator.opts.reltol = reltol
                else
                    if integrator.opts.abstol < 1e-10
                        integrator.opts.abstol *= 2.0
                    end
                end
            end
        end

        if (integrator.dt ./ tmin .>= 1)
            set_proposed_dt!(integrator, tmin .* 0.1)
        elseif (integrator.dt ./ tmin .<= 1e-3) && (wait % 1000 == 0)
            set_proposed_dt!(integrator, integrator.dt .* 1.03)
        elseif ((integrator.dt ./ tmin .<= 1e-3) || (integrator.dt ./ integrator.t .<= 1e-4)) &&
                (wait % 50 == 0) && (wait > 5000)
            for i in 1:idx_lvl
                if reltol[i] < reltol_Thres
                    reltol[i] *= 1.2
                    integrator.opts.reltol = reltol
                else
                    if integrator.opts.abstol < 1e-10
                        integrator.opts.abstol *= 2.0
                    end
                end
            end
        elseif (integrator.dt .<= 1e-13)
            debug && println("affect_timescale!: dt = ", integrator.dt, " <= 1e-13 at t = ", integrator.t, "; stopping")
            terminate!(integrator)
        end
    end

    # Shared: wall-clock limit. Must be shorter than the Slurm wall so the
    # solver can terminate!(integrator), return, and write .dat output.
    # Override with NL_MAX_WALL_SEC (seconds). Default 4.5 days (astro3_long is 5d).
    max_real_time = let
        env = get(ENV, "NL_MAX_WALL_SEC", "")
        if isempty(env)
            4.5 * 24.0 * 60.0 * 60.0
        else
            parse(Float64, env)
        end
    end
    start_time = Dates.now()
    println("ODE wall-clock limit: $(max_real_time / 3600.0) h (then save and exit)")
    flush(stdout)

    function time_limit_callback(u, t, integrator)
        elapsed_time = Dates.now() - start_time
        if Dates.value(elapsed_time) > max_real_time * 1e3
            println("Terminating integration due to wall-clock limit at t=", t,
                    " after ", Dates.value(elapsed_time) / 3.6e6, " h")
            flush(stdout)
            return true
        else
            return false
        end
    end

    function affect_time!(integrator)
        terminate!(integrator)
    end

    # Shared: enforce occupation numbers >= initial value
    log_e_init = log(e_init)
    function check_lower_bound(u, t, integrator)
        for i in 1:idx_lvl
            if u[i] < log_e_init
                return true
            end
        end
        return false
    end

    function affect_lower_bound!(integrator)
        for i in 1:idx_lvl
            if integrator.u[i] < log_e_init
                integrator.u[i] = log_e_init
            end
        end
    end

    # :radau only: a step that ends with a level more than an e-fold above its bosenova cap is
    # rejected and retried with a smaller dt. Above the cap the RHS is clamped flat (X = bn, no
    # growth), so the error estimate cannot see such an overshoot: at BH 10, f_a = 1e12,
    # alpha = 1.3 one accepted step took (n,8,8) from ~1e-14 to ~1e190 and the run then stalled.
    # :trbdf2 keeps the original behaviour.
    radau_out_of_domain(u, p, t) = any(i -> u[i] > ln_bn_vec[i] + 1.0, 1:idx_lvl)

    # ============================================================================
    # BUILD CALLBACK SET
    # ============================================================================
    def_spin_tol = 1e-3
    dt_guess = abs.((maximum(SR_rates) ./ hbar .* 3.15e7)^(-1) ./ 5.0)
    cback_lower = DiscreteCallback(check_lower_bound, affect_lower_bound!, save_positions=(false, false))

    # :radau only: the saveat points, and the save-on-change points inside a step, are taken
    # here (cback_saveat, affect_change!) from the step's interpolant with every level bounded
    # by its values at the two ends of the step, instead of by solve(saveat). The RadauIIA5
    # interpolant runs through the stage values, which are not accurate for a stiff level held
    # in quasi-equilibrium (only the step ends are): at BH 10, f_a = 1e14, alpha = 1.6 the
    # accepted steps kept (15,10,10) at ln N = -32.1 while the interpolant inside the same steps
    # swung between -42 and -16, and in other steps far above the bosenova cap, up to overflow
    # (~1e5 saved points). Only the saved points change; the integration does not.
    bounded_saves = solver === :radau && !spinone
    function save_interp(integrator, tk)
        tk == integrator.t && return copy(integrator.u)
        v = integrator(tk)
        if bounded_saves
            for i in 1:idx_lvl
                a, b = integrator.uprev[i], integrator.u[i]
                v[i] = clamp(v[i], min(a, b), max(a, b))
            end
        end
        return v
    end
    next_save = Ref(1)
    check_saveat(u, t, integrator) = next_save[] <= length(saveat) && saveat[next_save[]] <= t
    function affect_saveat!(integrator)
        sol = integrator.sol
        while next_save[] <= length(saveat) && saveat[next_save[]] <= integrator.t
            tk = saveat[next_save[]]
            if (isempty(sol.t) || tk > sol.t[end]) && integrator.saveiter == length(sol.t)
                push!(sol.t, tk)
                push!(sol.u, save_interp(integrator, tk))
                integrator.saveiter += 1
            end
            next_save[] += 1
        end
        u_modified!(integrator, false)
    end
    cback_saveat = DiscreteCallback(check_saveat, affect_saveat!, save_positions=(false, false))

    # Save-on-change (see docstring), in log space against the last saved point
    # (saveat or callback). change_units = change in units of the thresholds; only
    # the part of a level's change above save_rel_floor * (largest occupation)
    # counts. First in the callback set: the other saving callbacks fire on exactly
    # the fast steps, and would otherwise save the end point first and hide the
    # step from this one. It never modifies u, so the integration is unaffected;
    # points it adds inside a step can precede pending saveat points and are put in
    # time order on output.
    ln_save_factor = log(save_factor)
    ln_save_floor = log(save_rel_floor)
    function change_units(ua, ub)
        thr = max(maximum(@view ua[1:idx_lvl]), maximum(@view ub[1:idx_lvl])) + ln_save_floor
        dl = 0.0
        for i in 1:idx_lvl
            dl = max(dl, abs(max(ua[i], thr) - max(ub[i], thr)))
        end
        return max(dl / ln_save_factor, abs(ua[spinI] - ub[spinI]) / save_dlna)
    end
    function check_change(u, t, integrator)
        sol = integrator.sol
        isempty(sol.t) && return false
        sol.t[end] >= t && return false                    # already saved at this time
        return change_units(u, sol.u[end]) > 1
    end
    function affect_change!(integrator)
        sol = integrator.sol
        t1 = integrator.t
        t0 = max(sol.t[end], integrator.tprev)             # the interpolant covers [tprev, t]
        u0 = t0 == sol.t[end] ? sol.u[end] : save_interp(integrator, t0)
        n = min(ceil(Int, change_units(integrator.u, u0)), save_max_sub)
        if n > 1 && integrator.saveiter == length(sol.t)
            for k in 1:(n - 1)
                tk = t0 + (t1 - t0) * k / n
                push!(sol.t, tk)
                push!(sol.u, save_interp(integrator, tk))
                integrator.saveiter += 1
            end
        end
        u_modified!(integrator, false)                     # the end point is saved by save_positions
    end
    cback_change = DiscreteCallback(check_change, affect_change!, save_positions=(false, true))

    if spinone
        # Spinone: minimal callbacks
        cbackspin = DiscreteCallback(check_spin, affect_spin!, save_positions=(false, true))
        cbset = CallbackSet(cback_change, cbackspin, cback_lower)
    else
        # Standard: full callback set
        callbackTIME = DiscreteCallback(time_limit_callback, affect_time!, save_positions=(false, false))
        # affect_timescale! never changes u. With :radau tell the integrator so
        # (a DiscreteCallback otherwise flags u as modified, which marks the
        # Jacobian stale); :trbdf2 keeps the original behaviour.
        cbackdt = solver === :radau ?
            DiscreteCallback(check_timescale, itg -> (affect_timescale!(itg); u_modified!(itg, false)),
                             save_positions=(false, false)) :
            DiscreteCallback(check_timescale, affect_timescale!, save_positions=(false, true))
        cbackspin = DiscreteCallback(check_spin, affect_spin!, save_positions=(false, true))
        cbset = bounded_saves ? CallbackSet(cback_saveat, cback_change, cbackspin, cbackdt, callbackTIME, cback_lower) :
                                CallbackSet(cback_change, cbackspin, cbackdt, callbackTIME, cback_lower)

    end

    # Minimum step. The OrdinaryDiffEq default, dtmin = eps(t_end), is ~2e-6 yr for t_end = 1e10 yr: longer than
    # the fastest transients (e.g. the 322 burst at alpha = 0.2, f_a = 1e16 GeV, e-folding ~1 min), which then
    # abort with DtLessThanMin. Use no fixed floor (dtmin = floatmin) and instead stop once a step can no longer
    # advance t, i.e. dt is below the floating-point resolution at the current t.
    cback_stall = DiscreteCallback((u, t, integrator) -> integrator.t + integrator.dt == integrator.t,
                                   integrator -> (debug && println("dt = ", integrator.dt, " no longer advances t = ",
                                                                   integrator.t, "; stopping"); terminate!(integrator)),
                                   save_positions=(false, false))
    cbset = CallbackSet(cbset, cback_stall)

    # ============================================================================
    # SOLVE ODE
    # ============================================================================
    if spinone
        # Spinone uses fixed reltol
        prob = ODEProblem(RHS_ax!, y0, tspan, Mvars, reltol=1e-7, abstol=1e-10)
    else
        # Standard uses adaptive reltol array
        # println(reltol)
        # prob = ODEProblem(RHS_ax!, y0, tspan, Mvars, reltol=reltol, abstol=1e-10)
        prob = ODEProblem(analytic_jac ? ODEFunction(RHS_ax!; jac=jac_ax!) : RHS_ax!, y0, tspan, Mvars,
                          reltol=reltol, abstol=1e-10)
    end
    solver in (:trbdf2, :radau) || error("solver must be :trbdf2 or :radau")
    (analytic_jac && spinone) && error("analytic_jac is not implemented for spinone mode")
    alg = solver === :radau ? RadauIIA5(autodiff=false) : TRBDF2(autodiff=false)
    t_solve0 = time()
    domain_kw = (solver === :radau && !spinone) ? (isoutofdomain=radau_out_of_domain,) : NamedTuple()
    # :radau saves the saveat grid itself (cback_saveat, bounded interpolant). As with solve(saveat), t_start
    # (not on the grid) is not saved, and the end point is (t_max, or the wall-clock stop)
    save_kw = bounded_saves ? (saveat=Float64[], save_everystep=false, save_start=false, save_end=true) : (saveat=saveat,)
    sol = solve(prob, alg, dt=dt_guess, dtmin=floatmin(Float64), callback=cbset,
                maxiters=5e6; save_kw..., domain_kw...)
    debug && println("solver retcode: ", sol.retcode, "  t_end = ", sol.t[end], " of ", t_max,
                     "  [setup ", round(t_solve0 - t_setup0, digits=1), " s, integration ",
                     round(time() - t_solve0, digits=1), " s, ", sol.stats.naccept, " steps, ",
                     sol.stats.nf, " RHS calls]")
    # ============================================================================
    # EXTRACT AND PROCESS OUTPUT
    # ============================================================================
    # Save-on-change can add points inside a step before that step's saveat points;
    # a stable sort restores time order (identity when nothing was added).
    perm = sortperm(sol.t, alg=MergeSort)
    sol_t = sol.t[perm]
    sol_u = sol.u[perm]
    state_out = []
    for j in 1:idx_lvl
        push!(state_out, [exp(sol_u[i][j]) for i in 1:length(sol_u)])
    end


    spinBH = [exp(sol_u[i][spinI]) for i in 1:length(sol_u)]
    MassB = [exp(sol_u[i][massI]) for i in 1:length(sol_u)]
    if return_all_info
        return sol_t, state_out, modes, spinBH, MassB
    end

    # Check for incomplete evolution
    if spinone
        if (sol_t[end] != t_max)
            return 0.0, MassB[end]
        end
    else
        if (sol_t[end] != t_max) && (spinBH[end] > stop_on_a)
            return 0.0, MassB[end]
        end
    end

    # Handle NaN and Inf
    if isnan(spinBH[end])
        spinBH = spinBH[.!isnan.(spinBH)]
    end
    if isinf(spinBH[end])
        spinBH = spinBH[.!isinf.(spinBH)]
    end

    if isnan(MassB[end])
        MassB = MassB[.!isnan.(MassB)]
    end
    if isinf(MassB[end])
        MassB = MassB[.!isinf.(MassB)]
    end

    return spinBH[end], MassB[end]

end
