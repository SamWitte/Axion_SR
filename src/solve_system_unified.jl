using DelimitedFiles
using Interpolations
include(joinpath(@__DIR__, "Core", "accretion.jl"))


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
- `high_p::Bool`: Use high-precision tolerances (default true)
- `Nmax::Int`: Maximum principal quantum number 3-8 (default 3)
- `cheby::Bool`: Use Chebyshev interpolation (default true)
- `f_edd::Float64`: Eddington ratio of constant thin-disc accretion (default 0 = off);
  Mdot = f_edd * Mdot_Edd(M_BH) is fixed in time (Core/accretion.jl)
- `acc_eta::Float64`: radiative efficiency defining Mdot_Edd = L_Edd/(eta c^2) (default 0.1)

- `track_alpha`: let the SR rates (tabulated on a grid in M) and the scattering
  rates (interpolated in log alpha between grid masses) follow M(t). Default
  (`nothing`) = on iff f_edd > 0, since accretion can change M by O(1); with it
  off they stay at the initial alpha, as before. The GW rates always follow M(t);
  the bosenova caps and Emax2 always use the initial alpha.

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
    gw_model=:nonrel, gw_min_rate_per_yr=1e-10, gw_literature=true,
    f_edd=0.0, acc_eta=0.1, track_alpha=nothing)

    # ============================================================================
    # PARAMETER SETUP & VALIDATION
    # ============================================================================
    alph = GNew .* M_BH .* mu
    accreting = f_edd > 0
    if accreting && spinone
        error("accretion (f_edd > 0) is not implemented for spinone mode")
    end
    Mdot_acc = accreting ? f_edd * eddington_rate(M_BH; eta=acc_eta) : 0.0  # M_sun/yr, constant
    track_alpha = (track_alpha === nothing) ? accreting : track_alpha
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
        # Under accretion, levels that are not superradiant at t=0 may become so
        # later, so their scattering rates are kept.
        SR_kill = accreting ? ones(length(SR_rates)) : SR_rates
        rates = load_rate_coeffs(mu, M_BH, aBH, fa, Nmax, SR_kill; non_rel=non_rel, lm_only=lm_only)
        Mvars = [mu, fa, Emax2, aBH, M_BH, impose_low_cut]
        rP_initial = 1.0 + sqrt(1.0 - aBH^2)
    end

    # ============================================================================
    # TRACK_ALPHA: rate grids in BH mass so that alpha can follow M(t)
    # ============================================================================
    # M_lo allows for SR spin-down mass loss; M_hi bounds the mass reached by
    # t_max. SR rates: log-spaced nodes (1%), linear in ln M between them.
    # Scattering rates: coarser nodes (5%), linear in ln(rate) vs ln M.
    M_hi = track_alpha ? 1.02 * (M_BH + Mdot_acc * t_max) : M_BH
    lnM_lo, nM_sr, dlnM_sr, M_sr, sr_nodes, nM_sc, dlnM_sc, rates_sc = if track_alpha
        M_lo = 0.8 * M_BH
        n_sr = max(2, ceil(Int, log(M_hi / M_lo) / 0.01) + 1)
        d_sr = log(M_hi / M_lo) / (n_sr - 1)
        Ms_sr = exp.(log(M_lo) .+ d_sr .* (0:(n_sr - 1)))
        n_sc = max(2, ceil(Int, log(M_hi / M_lo) / 0.05) + 1)
        d_sc = log(M_hi / M_lo) / (n_sc - 1)
        Ms_sc = exp.(log(M_lo) .+ d_sc .* (0:(n_sc - 1)))
        if debug
            println("track_alpha: f_edd=$(f_edd), Mdot=$(Mdot_acc) M_sun/yr; SR grid $(n_sr) masses, ",
                    "scattering grid $(n_sc) masses in [$(M_lo), $(M_hi)] M_sun")
        end
        (log(M_lo), n_sr, d_sr, collect(Ms_sr),
         [compute_sr_rates_smooth(modes, Mj, aBH, GNew * Mj * mu, cheby=cheby)[2] for Mj in Ms_sr],
         n_sc, d_sc,
         [load_rate_coeffs(mu, Mj, aBH, fa, Nmax, SR_kill; non_rel=non_rel, lm_only=lm_only) for Mj in Ms_sc])
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

    # GW annihilations/transitions for every pair of levels; gw_rhs! re-evaluates
    # them each RHS call at the current BH mass and spin (single assignment keeps
    # the RHS closure type-stable).
    gw_cache = spinone ? gw_empty_cache(Tuple{Int,Int,Int}[]) :
               gw_build_cache(Nmax, modes, mu, M_BH, aBH; gw_model=gw_model,
                              min_rate_per_yr=gw_min_rate_per_yr, literature_overrides=gw_literature,
                              M_cut=M_hi)
    if debug && !spinone
        println("GW ($(gw_model)): $(length(gw_cache.ann_i)) annihilation + $(length(gw_cache.tr_i)) transition channels")
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
            if (u_real[1] .> Emax2) && (SR_rates_local[1] > 0)
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
        sanitize_state!(u_real)
        

        SR_rates_local, should_zero = compute_SR_rates_local(u_real)
        
        if spinone && should_zero
            du .*= 0.0
            return
        end

        # Superradiance terms
        du[spinI] = 0.0
        du[massI] = 0.0

        for i in 1:idx_lvl
            du[i] = SR_rates_local[i] .* u_real[i] ./ mu
            du[spinI] += -m_list[i] * SR_rates_local[i] .* u_real[i] ./ mu
            du[massI] += -SR_rates_local[i] .* u_real[i] ./ mu
        end

        # Scattering terms (standard mode only)
        if !spinone
            rP_now = 1.0 + sqrt(1.0 - u_real[spinI]^2)
            rP_ratio_now = rP_now / rP_initial
            a_now = u_real[spinI]
            omH_now = a_now / (rP_now^2 + a_now^2)
            alph_now = GNew * u_real[massI] * mu
            if track_alpha
                jsc, wsc = mass_bracket(u_real[massI], nM_sc, dlnM_sc)
            end
            for (kk, (idxV, sgn, is_bh_final, base_rate, m_drag)) in enumerate(rate_cache)
                if track_alpha && base_rate != 0.0
                    base_rate = sign(base_rate) * exp((1 - wsc) * lnR_sc[kk, jsc] + wsc * lnR_sc[kk, jsc + 1])
                end
                u_term_tot = 1.0
                for j in 1:length(sgn)
                    if (idxV[j] <= idx_lvl) && (idxV[j] > 0)
                        u_term_tot *= u_real[idxV[j]]
                    end
                end
                kH = (is_bh_final && m_drag != 0) ? (alph_now - m_drag * omH_now) : 1.0
                rate_val = base_rate * (is_bh_final ? rP_ratio_now : 1.0) * kH
                for j in 1:length(sgn)
                    idx_j = idxV[j]
                    if idx_j == 0; continue; end
                    if idx_j == -1; idx_j = massI; end
                    du[idx_j] += sgn[j] * rate_val * u_term_tot
                end
            end
            gw_rhs!(du, u_real, gw_cache, alph_now, a_now)
        end

        # Unit corrections
        for i in 1:idx_lvl
            if spinone
                if u_real[i] < e_init
                    du[i] = 0.0
                else
                    du[i] *= mu ./ hbar .* 3.15e7
                end
            else
                if ((abs(u[i] - log(bn_list[i])) < SOLVER_TOLERANCES.bosenova_threshold) ||
                    (u[i] > log(bn_list[i]))) && (du[i] > 0)
                    du[i] = 0.0
                elseif (abs(u[i] - log(e_init)) < SOLVER_TOLERANCES.bosenova_threshold) && (du[i] < 0)
                    du[i] = 0.0
                else
                    du[i] *= mu ./ hbar .* 3.15e7
                end
            end
        end

        du[spinI] *= mu ./ hbar .* 3.15e7
        du[massI] *= (mu .* u_real[massI]) .* (mu .* GNew .* u_real[massI]) ./ hbar .* 3.15e7

        if accreting   # already per year
            dM_acc, da_acc = accretion_rhs(u_real[massI], u_real[spinI], Mdot_acc)
            du[massI] += dM_acc
            du[spinI] += da_acc
        end

        du ./= u_real

        for i in 1:idx_lvl
            if turn_off[i]
                du[i] = 0.0
            end
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

    # ============================================================================
    # BUILD CALLBACK SET
    # ============================================================================
    def_spin_tol = 1e-3
    dt_guess = abs.((maximum(SR_rates) ./ hbar .* 3.15e7)^(-1) ./ 5.0)
    cback_lower = DiscreteCallback(check_lower_bound, affect_lower_bound!, save_positions=(false, false))
    if spinone
        # Spinone: minimal callbacks
        cbackspin = DiscreteCallback(check_spin, affect_spin!, save_positions=(false, true))
        cbset = CallbackSet(cbackspin, cback_lower)
    else
        # Standard: full callback set
        callbackTIME = DiscreteCallback(time_limit_callback, affect_time!, save_positions=(false, false))
        cbackdt = DiscreteCallback(check_timescale, affect_timescale!, save_positions=(false, true))
        cbackspin = DiscreteCallback(check_spin, affect_spin!, save_positions=(false, true))
        cbset = CallbackSet(cbackspin, cbackdt, callbackTIME, cback_lower)

    end

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
        prob = ODEProblem(RHS_ax!, y0, tspan, Mvars, reltol=reltol, abstol=1e-10)
    end
    sol = solve(prob, TRBDF2(autodiff=false), dt=dt_guess, saveat=saveat, callback=cbset, maxiters=5e6)
    # ============================================================================
    # EXTRACT AND PROCESS OUTPUT
    # ============================================================================
    state_out = []
    for j in 1:idx_lvl
        push!(state_out, [exp(sol.u[i][j]) for i in 1:length(sol.u)])
    end


    spinBH = [exp(sol.u[i][spinI]) for i in 1:length(sol.u)]
    MassB = [exp(sol.u[i][massI]) for i in 1:length(sol.u)]
    if return_all_info
        return sol.t, state_out, modes, spinBH, MassB
    end

    # Check for incomplete evolution
    if spinone
        if (sol.t[end] != t_max)
            return 0.0, MassB[end]
        end
    else
        if (sol.t[end] != t_max) && (spinBH[end] > stop_on_a)
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
