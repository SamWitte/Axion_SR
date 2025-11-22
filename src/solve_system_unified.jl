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
    default_reltol = 1e-7

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
        if !(3 <= Nmax <= 8)
            error("Nmax must be between 3 and 8 for standard mode, got $Nmax")
        end
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
        interp_funcs = []
    else
        # Standard: Compute interpolated rates
        SR_rates, interp_funcs, interp_dict = compute_sr_rates(modes, M_BH, aBH, alph, cheby=cheby)
        rates = load_rate_coeffs(mu, M_BH, aBH, fa, Nmax, SR_rates; non_rel=non_rel)
        Mvars = [mu, fa, Emax2, aBH, M_BH, impose_low_cut]
    end

    # ============================================================================
    # ODE SETUP
    # ============================================================================
    tspan = (0.0, t_max)
    saveat = (tspan[2] .- tspan[1]) ./ n_times

    # Trackers for callbacks
    wait = 0
    turn_off = fill(false, idx_lvl)
    turn_off_M = false

    # ============================================================================
    # RHS FUNCTION (shared logic with mode-specific branches)
    # ============================================================================
    function RHS_ax!(du, u, Mvars, t)
        u_real = exp.(u)

        # Unpack Mvars based on mode
        if spinone
            mu, aBH_i, M_BH_i = Mvars
        else
            mu, fa, Emax2, aBH_i, M_BH_i, impose_low_cut = Mvars
        end

        # ====================================================================
        # SPIN BOUNDARY CONDITION
        # ====================================================================
        if spinone
            # Spinone: Simple hard cutoff
            if u_real[spinI] > maxSpin
                u_real[spinI] = maxSpin
            end
        else
            # Standard: Compute rP factor (used elsewhere)
            if u_real[spinI] .> maxSpin
                u_real[spinI] = maxSpin
            elseif u_real[spinI] .< 0.0
                u_real[spinI] = 0.0
            end
        end

        # ====================================================================
        # BOSENOVA BOUNDARY (standard mode only)
        # ====================================================================
        if !spinone
            for i in 1:idx_lvl
                if u_real[i] < e_init
                    u_real[i] = e_init
                    u[i] = log(e_init)
                end

                if (abs(u[i] - log(bn_list[i])) < SOLVER_TOLERANCES.bosenova_threshold) ||
                   (u[i] > log(bn_list[i]))
                    u[i] = log(bn_list[i])
                    u_real[i] = bn_list[i]
                end
            end
        end

        # ====================================================================
        # RATE COMPUTATION (mode-specific)
        # ====================================================================
        OmegaH = u_real[spinI] ./ (2 .* (GNew .* u_real[massI]) .* (1 .+ sqrt.(1 .- u_real[spinI].^2)))

        if spinone
            wR, wI = precomputed_spin1(alph, u_real[spinI], u_real[massI])
            if wR .> OmegaH
                du .*= 0.0
                return
            end
            SR_rates_local = [2 .* wI]

            # Spinone saturation condition
            if u_real[1] .>= u_real[spinI]
                u[1] = log(u_real[spinI])
                u_real[1] = u_real[spinI]
                SR_rates_local = [0.0]
            end
        else
            # Standard: Interpolated rates
            SR_rates_local = [func(u_real[spinI]) for func in interp_funcs]

            # Emax2 cutoff for 211 level
            if (u_real[1] .> Emax2) && (SR_rates_local[1] > 0)
                SR_rates_local[1] *= 0.0
            end
        end

        # ====================================================================
        # SUPERRADIANCE (SR) TERMS
        # ====================================================================
        du[spinI] = 0.0
        du[massI] = 0.0
        for i in 1:idx_lvl
            du[i] = SR_rates_local[i] .* u_real[i] ./ mu
            du[spinI] += -m_list[i] * SR_rates_local[i] .* u_real[i] ./ mu
            du[massI] += -SR_rates_local[i] .* u_real[i] ./ mu
        end

        # ====================================================================
        # SCATTERING TERMS (standard mode only)
        # ====================================================================
        if !spinone
            rate_keys = collect(keys(rates))
            for i in 1:length(rate_keys)
                idxV, sgn = key_to_indx(rate_keys[i], Nmax)
                u_term_tot = 1.0

                for j in 1:length(sgn)
                    if (idxV[j] <= idx_lvl) && (idxV[j] > 0)
                        u_term_tot *= u_real[idxV[j]]
                    end
                end

                for j in 1:length(sgn)
                    if idxV[j] == 0
                        continue
                    end
                    if idxV[j] == -1
                        idxV[j] = massI
                    end

                    du[idxV[j]] += sgn[j] * rates[rate_keys[i]] * u_term_tot
                end
            end
        end

        # ====================================================================
        # UNIT CORRECTIONS & BOUNDARY CHECKS
        # ====================================================================
        for i in 1:idx_lvl
            if spinone
                # Spinone: simpler energy floor check
                if u_real[i] < e_init
                    u_real[i] = e_init
                    u[i] = log(e_init)
                end
                du[i] *= mu ./ hbar .* 3.15e7
            else
                # Standard: bosenova boundary checks
                if ((abs(u[i] - log(bn_list[i])) < SOLVER_TOLERANCES.bosenova_threshold) ||
                    (u[i] > log(bn_list[i]))) && (du[i] > 0)
                    du[i] *= 0.0
                elseif (abs(u[i] - log(e_init)) < SOLVER_TOLERANCES.bosenova_threshold) && (du[i] < 0)
                    du[i] *= 0.0
                else
                    du[i] *= mu ./ hbar .* 3.15e7
                end
            end
        end

        # Spin and mass unit corrections (shared)
        du[spinI] *= mu ./ hbar .* 3.15e7
        du[massI] *= (mu .* u_real[massI]) .* (mu .* GNew .* u_real[massI]) ./ hbar .* 3.15e7

        du ./= u_real

        # Apply turn_off flags
        for i in 1:idx_lvl
            if turn_off[i]
                du[i] *= 0.0
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
        du = get_du(integrator)
        u_real = exp.(u)
        rate_keys = collect(keys(rates))

        all_contribs = zeros(idx_lvl)
        test = zeros(idx_lvl)
        u_fake = u_real * 1.1

        for i in 1:idx_lvl
            if (abs(u[i] - log(bn_list[i])) < SOLVER_TOLERANCES.bosenova_threshold) ||
               (u[i] > log(bn_list[i]))
                u[i] = log(bn_list[i])
                u_real[i] = bn_list[i]
                du[i] *= 0
            end
            all_contribs[i] = SR_rates[i] .* u_real[i] ./ mu
            test[i] = SR_rates[i] .* u_fake[i] ./ mu
        end

        for i in 1:idx_lvl
            if (u[i] .< log(1e-75)) && (SR_rates[i] < 0)
                turn_off[i] = true
            end
        end
        if u_real[massI] > (1.4 * M_BH)
            turn_off_M = true
        end

        for i in 1:length(rate_keys)
            idxV, sgn = key_to_indx(rate_keys[i], Nmax)

            u_term_tot = 1.0
            u_term_tot_fake = 1.0
            for j in 1:length(sgn)
                if (idxV[j] <= idx_lvl) && (idxV[j] > 0)
                    if j == 3
                        u_term_tot *= (u_real[idxV[j]] .+ e_init)
                        u_term_tot_fake *= (u_fake[idxV[j]] .+ e_init)
                    else
                        u_term_tot *= u_real[idxV[j]]
                        u_term_tot_fake *= u_fake[idxV[j]]
                    end
                end
            end

            for j in 1:length(sgn)
                if (idxV[j] <= 0)
                    continue
                end
                all_contribs[idxV[j]] += sgn[j] * rates[rate_keys[i]] * u_term_tot
                test[idxV[j]] += sgn[j] * rates[rate_keys[i]] * u_term_tot_fake
            end
        end

        integrator.opts.reltol = reltol

        tlist = []
        for i in 1:idx_lvl
            condBN = (abs(u[i] - log(bn_list[i])) < SOLVER_TOLERANCES.bosenova_threshold)

            if (u[i] > log(e_init)) && condBN && (du[i] != 0.0)
                append!(tlist, abs(1.0 ./ du[i]))
            end
        end

        def_spin_tol = 1e-3
        append!(tlist, def_spin_tol ./ du[spinI])
        tmin = minimum(abs.(tlist))

        if (integrator.dt ./ tmin .>= 0.1)
            return true
        elseif (integrator.dt ./ tmin .<= 0.001)
            return true
        elseif (integrator.dt .<= 1e-10)
            return true
        else
            return false
        end
    end

    function affect_timescale!(integrator)
        du = get_du(integrator)
        tlist = []
        indx_list = []
        for i in 1:idx_lvl
            condBN = (abs(integrator.u[i] - log(bn_list[i])) < SOLVER_TOLERANCES.bosenova_threshold)
            if (integrator.u[i] > log(e_init)) && condBN
                append!(tlist, (1.0 ./ du[i]))
                append!(indx_list, i)
            end
        end

        def_spin_tol = 1e-3
        append!(tlist, def_spin_tol ./ du[spinI])
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
    max_real_time = 20.0 * 60  # Convert to seconds
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

    # ============================================================================
    # BUILD CALLBACK SET
    # ============================================================================
    def_spin_tol = 1e-3
    dt_guess = (maximum(SR_rates) ./ hbar .* 3.15e7)^(-1) ./ 5.0

    if spinone
        # Spinone: minimal callbacks
        cbackspin = DiscreteCallback(check_spin, affect_spin!, save_positions=(false, true))
        cbset = CallbackSet(cbackspin)
    else
        # Standard: full callback set
        callbackTIME = DiscreteCallback(time_limit_callback, affect_time!, save_positions=(false, false))
        cbackdt = DiscreteCallback(check_timescale, affect_timescale!, save_positions=(false, true))
        cbackspin = DiscreteCallback(check_spin, affect_spin!, save_positions=(false, true))
        cbset = CallbackSet(cbackspin, cbackdt, callbackTIME)
    end

    # ============================================================================
    # SOLVE ODE
    # ============================================================================
    if spinone
        # Spinone uses fixed reltol
        prob = ODEProblem(RHS_ax!, y0, tspan, Mvars, reltol=1e-7, abstol=1e-10)
    else
        # Standard uses adaptive reltol array
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
