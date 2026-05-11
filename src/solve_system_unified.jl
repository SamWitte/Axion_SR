using DelimitedFiles
using Interpolations


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
    N_pts_interp=200, N_pts_interpL=200, Nmax=3, cheby=true, spinone=false)

    # ============================================================================
    # PARAMETER SETUP & VALIDATION
    # ============================================================================
    alph = GNew .* M_BH .* mu

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
        rates = load_rate_coeffs(mu, M_BH, aBH, fa, Nmax, SR_rates; non_rel=non_rel)
        Mvars = [mu, fa, Emax2, aBH, M_BH, impose_low_cut]
        rP_initial = 1.0 + sqrt(1.0 - aBH^2)
    end

    # ============================================================================
    # PRE-COMPUTE RATE INDEX CACHE (Fix 1: avoid O(N_modes × N_rates) per RHS call)
    # ============================================================================
    # key_to_indx is O(N_modes) with string allocations per call. For Nmax=15 with
    # ~22k rate keys × 560 modes this was ~18 min/Jacobian. Cache it once here.
    rate_cache = if !spinone
        map(collect(keys(rates))) do k
            if abs.(rates[k]) .> 1e20
                rates[k] = 0.0
            end
            idxV, sgn = key_to_indx(k, Nmax)
            is_bh = any(idxV .== -1)
            (idxV, sgn, is_bh, rates[k])
        end
    else
        Vector{Tuple{Vector{Int}, Vector{Int}, Bool, Float64}}()
    end

    # ============================================================================
    # ODE SETUP
    # ============================================================================
    tspan = (0.0, t_max)
    # Log-spaced save points so that early-time dynamics are resolved on the log x-axis.
    # t_start avoids log(0); the floor at 1.0 is arbitrary but safe for all tau_max values.
    t_log_start = max(1.0, tspan[2] * 1e-9)
    saveat = exp10.(range(log10(t_log_start), log10(tspan[2]), length=n_times))

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
            SR_rates_local = [func(spin_val) for func in interp_funcs]
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
            for (idxV, sgn, is_bh_final, base_rate) in rate_cache
                u_term_tot = 1.0
                for j in 1:length(sgn)
                    if (idxV[j] <= idx_lvl) && (idxV[j] > 0)
                        u_term_tot *= u_real[idxV[j]]
                    end
                end
                rate_val = base_rate * (is_bh_final ? rP_ratio_now : 1.0)
                for j in 1:length(sgn)
                    idx_j = idxV[j]
                    if idx_j == 0; continue; end
                    if idx_j == -1; idx_j = massI; end
                    du[idx_j] += sgn[j] * rate_val * u_term_tot
                end
            end
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
            if u_real[spinI] .> (aBH .+ 0.01)
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
        if u_real[spinI] .> aBH
            integrator.u[spinI] = log(aBH)
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

    # Shared: time limit callback
    # max_real_time = 20.0 * 60  # Convert to seconds
    max_real_time = 10.0 .* 24.0 .* 60.0 .* 60.0
    start_time = Dates.now()

    function time_limit_callback(u, t, integrator)
        elapsed_time = Dates.now() - start_time
        if Dates.value(elapsed_time) > max_real_time * 1e3
            println("Terminating integration due to time limit")
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
    dt_guess = min(abs.((maximum(SR_rates) ./ hbar .* 3.15e7)^(-1) ./ 5.0), saveat)
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
