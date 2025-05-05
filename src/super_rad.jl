using Random
using OrdinaryDiffEq
using Statistics
using Distributions
using DelimitedFiles
using Dates
using Interpolations
include("Constants.jl")
include("solve_sr_rates.jl")
include("load_rates.jl")
using Printf

function super_rad_check(M_BH, aBH, massB, f_a; spin=0, tau_max=1e4, alpha_max_cut=100.0, debug=false, impose_low_cut=0.01, stop_on_a=0, eq_threshold=1e-100, abstol=1e-30, non_rel=true, high_p=true, N_pts_interp=15, N_pts_interpL=10, Nmax=3, cheby=false)
   
    alph = GNew .* M_BH .* massB #
    if debug
        print("Alpha \t", alph, "\n")
    end
    
    if alph .> alpha_max_cut
        if debug
            print("Need higher-level system... \n")
        end
        return aBH, M_BH
    elseif alph .< impose_low_cut
        return aBH, M_BH
    end
    

    print("Solving system... \n")

    final_spin, final_BH = solve_system(massB, f_a, aBH, M_BH, tau_max, debug=debug, impose_low_cut=impose_low_cut, stop_on_a=stop_on_a, eq_threshold=eq_threshold, abstol=abstol, non_rel=non_rel, high_p=high_p, N_pts_interp=N_pts_interp, N_pts_interpL=N_pts_interpL, Nmax=Nmax, cheby=cheby)

    
    print("Spin diff.. \t ", aBH, "\t", final_spin, "\t", alph, "\n")
    print("Mass diff.. \t ", M_BH, "\t", final_BH, "\t", alph, "\n")
    return final_spin, final_BH
    
end

function emax_211(MBH, mu, aBH)
    alph = GNew .* MBH .* mu
    
    emax_N = 1 .- 8 .* alph.^2 .+ 8 .* alph.^3 .* aBH .- sqrt.(abs.(1 .- 16 .* alph.^2 .+ 32 .* aBH .* alph.^3 - 16 .* aBH.^2 .* alph.^4)) # abs just in case...
    emax_D = 8 .* (-alph.^3 .+ aBH .* alph.^4)
    return (emax_N ./ emax_D)
end
    
function isapproxsigfigs(a, b, precision)
    return round(a, sigdigits=precision) == round(b, sigdigits=precision)
end

function solve_system(mu, fa, aBH, M_BH, t_max; n_times=10000, debug=false, impose_low_cut=0.01, return_all_info=false, eq_threshold=1e-100, stop_on_a=0, abstol=1e-30, non_rel=true, high_p=true, N_pts_interp=10, N_pts_interpL=5,  Nmax=3, cheby=false)

    alph = GNew .* M_BH .* mu
    if debug
        println("Entering...")
    end
    if non_rel
        default_reltol = 1e-5
    else
        default_reltol = 1e-3
    end
    if high_p
        reltol_Thres = 1e-3
    else
        reltol_Thres = 0.1
    end
    
    if (3 <= Nmax <= 8)
        idx_lvl = 0
        m_list = []
        bn_list = []
        modes = []
        
        e2_maxBN = 1024 * pi * (fa / M_pl).^2 ./ (9 * alph.^3)

        for nn in 1:Nmax, l in 1:(nn - 1),  m in 1:l
            idx_lvl += 1
            max_alph = aBH .* m ./ (2 .* (1 .+ sqrt.(1 .- aBH.^2))) .* 1.1
            push!(modes, (nn, l, m, max_alph))
            push!(m_list, m)
            push!(bn_list, e2_maxBN .* (nn ./ 2).^4)
        end
        
        ## need to add levels which truncate high m states
        for nn in 1:Nmax, l in 1:(nn - 1)
            m = l
            m_new = (2 * m)
            if m_new < Nmax
                continue # already included
            else
                idx_lvl += 1
                max_alph = aBH .* m_new ./ (2 .* (1 .+ sqrt.(1 .- aBH.^2))) .* 1.1
                push!(modes, (m_new + 1, m_new, m_new, max_alph))
                push!(m_list, m_new)
                push!(bn_list, e2_maxBN .* ((m_new + 1) ./ 2).^4)
               
            end
        end

    end

    ### only for testing...
    # default_reltol = 1e-7
    # reltol_Thres = 1e-7

    spinI = idx_lvl + 1
    massI = spinI + 1
    
    
    e_init = 1.0 ./ (GNew .* M_BH.^2 .* M_to_eV) # unitless
    
    
    y0 = []
    reltol = []
    for i in 1:idx_lvl
        append!(y0, e_init)
        append!(reltol, default_reltol)
    end
    append!(y0, aBH)
    append!(y0, M_BH)
    y0 = log.(y0)
    
    def_spin_tol = 1e-5
    append!(reltol, def_spin_tol)
    append!(reltol, def_spin_tol)
    
    wait = 0 # tracker
    
    
    Emax2 = 1.0
    OmegaH = aBH ./ (2 .* (GNew .* M_BH) .* (1 .+ sqrt.(1 .- aBH.^2)))
    if (OmegaH .> ergL(2, 1, 1, mu, M_BH, aBH))
        Emax2 = emax_211(M_BH, mu, aBH)
    end
    
    Mvars = [mu, fa, Emax2, aBH, M_BH, impose_low_cut]
    tspan = (0.0, t_max)
    saveat = (tspan[2] .- tspan[1]) ./ n_times
    
    
    xtol_slv = 1e-15
    iter_slv = 50
    
    SR_rates = zeros(idx_lvl)
    
    
    idx_fill = 1
    delt_a = 0.0001
    
    SR_rates, interp_funcs, interp_dict = compute_sr_rates(modes, M_BH, aBH, alph, cheby=cheby);
    
    if debug
        println(idx_lvl)
        println(modes)
        print("SR rates @ prod \t", SR_rates, "\n")
    end
    

    rates = load_rate_coeffs(mu, M_BH, aBH, fa, Nmax, SR_rates; non_rel=non_rel)
    
    if debug
        println("Rates loaded...")
        println(rates)
    end
    
   
    turn_off = []
    turn_off_M = false
    for i in 1:idx_lvl
        append!(turn_off, false)
    end
    

    #### DEFINING FUNCTION FOR EVOLUTION ######
    function RHS_ax!(du, u, Mvars, t)
    
        # [e211, e322, ... , aBH, MBH]
        # mu [mass eV], fa[1/GeV]
        u_real = exp.(u)
        
        mu, fa, Emax2, aBH_i, M_BH_i, impose_low_cut  = Mvars

        rP = nothing
        if u_real[spinI] .> maxSpin
            rP = 1.0 .+ sqrt.(1 - maxSpin .^2)
            u_real[spinI] = maxSpin
        elseif u_real[spinI] .< 0.0
            rP = 2.0
            u_real[spinI] = 0.0
        else
            rP = 1.0 .+ sqrt.(1 - u_real[spinI].^2)
        end
        
        for i in 1:idx_lvl
            if u_real[i] < e_init
                u_real[i] = e_init
                u[i] = log.(e_init)
            end
              
            if (abs.(u[i] .- log.(bn_list[i])) < 1e-2)||(u[i] > log.(bn_list[i]))
                u[i] = log.(bn_list[i])
                u_real[i] = bn_list[i]
                # print(u[i], "\t", log.(bn_list[i]), "\n")
            end
        end
        
        
        OmegaH = u_real[spinI] ./ (2 .* (GNew .* u_real[massI]) .* (1 .+ sqrt.(1 .- u_real[spinI].^2)))
        
        SR_rates = [func(u_real[spinI]) for func in interp_funcs]
        
        if (u_real[1] .> Emax2)&&(SR_rates[1] > 0)
            SR_rates[1] *= 0.0
        end
        
    
        # SR terms
        du[spinI] = 0.0
        du[massI] = 0.0
        for i in 1:idx_lvl
            du[i] = SR_rates[i] .* u_real[i] ./ mu
            du[spinI] += - m_list[i] * SR_rates[i] .*  u_real[i] ./ mu
            du[massI] += - SR_rates[i] .*  u_real[i] ./ mu
        end
        
        # Scattering terms
        rate_keys = collect(keys(rates))
        for i in 1:length(rate_keys)
            idxV, sgn = key_to_indx(rate_keys[i], Nmax)
            u_term_tot = 1.0
            
            for j in 1:length(sgn)
                if (idxV[j] <= idx_lvl)&&(idxV[j] > 0)
                    u_term_tot *= u_real[idxV[j]]
                end
            end
            
            for j in 1:length(sgn)
                # energy lost to infinity?
                if idxV[j] == 0
                    continue
                end
                if idxV[j] == -1
                    idxV[j] = massI
                end
                
                du[idxV[j]] += sgn[j] * rates[rate_keys[i]] * u_term_tot
            end
        end
        
        

        # check bosenova and correct units
        for i in 1:idx_lvl
      
            if ((abs.(u[i] .- log.(bn_list[i])) < 1e-2)||(u[i] > log.(bn_list[i])))&&(du[i] > 0)
                du[i] *= 0.0
            elseif (abs.(u[i] .- log.(e_init)) < 1e-2)&&(du[i] < 0)
                du[i] *= 0.0
            else
                du[i] *= mu ./ hbar .* 3.15e7
            end
            
        end
        
        
        
        du[spinI] *= mu ./ hbar .* 3.15e7
        du[massI] *= (mu .*  u_real[massI]) .* (mu .* GNew .*  u_real[massI]) ./ hbar .* 3.15e7

        du ./= u_real
        for i in 1:idx_lvl
            if turn_off[i]
                du[i] *= 0.0
            end
            
        end
        
        return
    end
    
 
    
    #### make sure the integrator timescale doesnt grow too large
    function check_timescale(u, t, integrator)
        du = get_du(integrator)
        
        u_real = exp.(u)
        rate_keys = collect(keys(rates))

        all_contribs = zeros(idx_lvl)
        test = zeros(idx_lvl)
        u_fake = u_real * 1.1 # was 2!
        
        for i in 1:idx_lvl
            if (abs.(u[i] .- log.(bn_list[i])) < 1e-2)||(u[i] > log.(bn_list[i]))
                u[i] = log.(bn_list[i])
                u_real[i] = bn_list[i]
                du[i] *= 0
            end
            all_contribs[i] = SR_rates[i] .* u_real[i] ./ mu
            test[i] = SR_rates[i] .* u_fake[i] ./ mu
        end
        
        for i in 1:idx_lvl
            if (u[i] .< log.(1e-75))&&(SR_rates[i] < 0)
                turn_off[i] = true
            end
        end
        if u_real[massI] > (1.4 * M_BH) # few safety nets in place...
            turn_off_M = true
        end
        
        
        for i in 1:length(rate_keys)
            idxV, sgn = key_to_indx(rate_keys[i], Nmax)
            
            u_term_tot = 1.0
            u_term_tot_fake = 1.0
            for j in 1:length(sgn)
                if (idxV[j] <= idx_lvl)&&(idxV[j] > 0)
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
        

        
        integrator.opts.reltol =  reltol
        
        tlist = []
        # print(integrator.sol.u, "\n")
        for i in 1:idx_lvl
            condBN =  (abs.(u[i] .- log.(bn_list[i])) < 1e-2)
            
            if (u[i] > log.(e_init)) && condBN &&  (du[i] != 0.0)
                append!(tlist, abs.(1.0 ./ du[i]))
            end
        end
        
        append!(tlist, def_spin_tol ./ du[spinI])
        tmin = minimum(abs.(tlist))
       
        # print(du, "\n")
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
           
            condBN =  (abs.(integrator.u[i] .- log.(bn_list[i])) < 1e-2)
            if (integrator.u[i] > log.(e_init)) && condBN
                append!(tlist, (1.0 ./ du[i]))
                append!(indx_list, i)
            end
        end
            
        append!(tlist, def_spin_tol ./ du[spinI])
        tmin = minimum(abs.(tlist))
        
        
        if (integrator.dt ./ integrator.t < 1e-6)&&(wait % 1000 == 0)&&(wait > 10000)
            for i in 1:idx_lvl
                if reltol[i] < reltol_Thres
                    reltol[i] *= 1.2
                    integrator.opts.reltol =  reltol
                    
                else
                    if integrator.opts.abstol < 1e-10
                        integrator.opts.abstol *= 2.0
                    end
                end
            end
        end
        
       

        if (integrator.dt ./ tmin .>= 1)
            set_proposed_dt!(integrator, tmin .* 0.1)
        elseif (integrator.dt ./ tmin .<= 1e-3)&&(wait % 1000 == 0)
            set_proposed_dt!(integrator, integrator.dt .* 1.03)
        elseif ((integrator.dt ./ tmin .<= 1e-3)||(integrator.dt ./ integrator.t .<= 1e-4))&&(wait % 50 == 0)&&(wait > 5000)
            # print(integrator.dt ./ tmin, "\t", wait, "\t", reltol, "\t", integrator.opts.abstol,  "\n")
            for i in 1:idx_lvl
                if reltol[i] < reltol_Thres
                    reltol[i] *= 1.2
                    integrator.opts.reltol =  reltol
                else
                    if integrator.opts.abstol < 1e-10
                        integrator.opts.abstol *= 2.0
                    end
                end
            end
            
        elseif (integrator.dt .<= 1e-13)
            print("time step too small!! \n")
            terminate!(integrator)
        end
       
    end
    
    # if running stat analysis, can truncate evolution when spin is small
    function check_spin(u, t, integrator)
        wait += 1
        u_real = exp.(u)

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
    function affect_spin!(integrator)
        u_real = exp.(integrator.u)
        if u_real[spinI] <= stop_on_a
            terminate!(integrator)
        end
        if u_real[spinI] .> aBH
            integrator.u[spinI] = log.(aBH)
        elseif u_real[spinI] .< 0.0
            integrator.u[spinI] = -10.0
        end
        set_proposed_dt!(integrator, integrator.dt .* 0.3)
    end
    

    dt_guess = (maximum(SR_rates) ./ hbar .* 3.15e7).^(-1) ./ 5.0

    if debug
        print("Time guess \t", dt_guess, "\n")
    end
    
    
    max_real_time = 20.0 # Set the maximum allowed time for integration (in minutes)
    max_real_time *= 60 # convert to seconds
    start_time = Dates.now() # Initialize the start time
    
    # Define callback -- if integrater very stuck, don't let it go on for infinite time! it might try!
    function time_limit_callback(u, t, integrator)
        elapsed_time = Dates.now() - start_time
        if Dates.value(elapsed_time) > max_real_time * 1e3  # Convert seconds to milliseconds
            println("Terminating integration due to time limit")
            return true
        else
            return false
        end
    end
    function affect_time!(integrator)
        terminate!(integrator)
    end
    

    # Set up the user_data with the start time
    callbackTIME = DiscreteCallback(time_limit_callback, affect_time!, save_positions=(false, false))
    cbackdt = DiscreteCallback(check_timescale, affect_timescale!, save_positions=(false, true))
    cbackspin = DiscreteCallback(check_spin, affect_spin!, save_positions=(false, true))
 
    
    
    cbset = CallbackSet(cbackspin, cbackdt, callbackTIME)
    prob = ODEProblem(RHS_ax!, y0, tspan, Mvars, reltol=reltol, abstol=1e-10)
    sol = solve(prob, TRBDF2(autodiff=false), dt=dt_guess, saveat=saveat, callback=cbset, maxiters=5e6)
    # sol = solve(prob, Rosenbrock23(autodiff=false), dt=dt_guess, saveat=saveat, callback=cbset, maxiters=5e6, dtmin=(dt_guess / 1e5), force_dtmin=true)
    
    state_out = []
    for j in 1:idx_lvl
        push!(state_out, [exp.(sol.u[i][j]) for i in 1:length(sol.u)])
    end
    
    spinBH = [exp.(sol.u[i][spinI]) for i in 1:length(sol.u)]
    MassB = [exp.(sol.u[i][massI]) for i in 1:length(sol.u)]

    
    if debug
        idx_211 = get_state_idx("211", Nmax)
        idx_322 = get_state_idx("322", Nmax)
        print("Initial \t", state_out[idx_211][1], "\t", state_out[idx_322][1],  "\t", spinBH[1], "\t", MassB[1], "\n")
        print("Final \t", state_out[idx_211][end], "\t", state_out[idx_322][end], "\t", spinBH[end], "\t", MassB[end], "\n")
       
    end

    if return_all_info
        return sol.t, state_out, idx_lvl, spinBH, MassB
    end
    
    if (sol.t[end] != t_max)&&(spinBH[end] > stop_on_a)
        print("Fail? Final time only \t", sol.t[end], "\n")
        return 0.0, MassB[end]
    end
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


function compute_sr_rates(qtm_cfigs, M_BH, aBH, alph; delt_a=0.0001, cheby=false)
    # Define the quantum numbers to iterate over
    # Each entry is [n, l, m, alph_threshold]
    
    # Initialize output array
    SR_rates = zeros(length(qtm_cfigs))
    interp_functions = []
    
    # Process each configuration
    for (idx, (n, l, m, alph_threshold)) in enumerate(qtm_cfigs)
        if m > 5  # file not computed!
            LinearInterpolation([0.0, 1.0], [1e-100, 1e-100], extrapolation_bc=Interpolations.Line())
            push!(interp_functions, interp_func)
            SR_rates[idx] = interp_func(aBH)
            continue
        end
        
        # Get minimum guess and pre-computed rates
        amin_guess, run_high, a_list_high, out_high, run_low, a_list_low, out_low = 
            pre_computed_sr_rates(n, l, m, alph, M_BH; n_high=200, n_low=200, delt_a=delt_a, cheby=cheby)
        
        # Handle upper interpolation (values >= threshold)
        if alph < alph_threshold
            if !run_high
                a_list_high = [amin_guess, amin_guess * 1.01]
                out_high = [1e-100, 1e-100]
            elseif length(a_list_high) == 1.0
                a_list_high = [alist[1], alist[1] * 1.01]
                out_high = [out_high[1], out_high[1]]
            end
        else
            amin_guess = 1.0
            a_list_high = [amin_guess, amin_guess * 1.01]
            out_high = [1e-100, 1e-100]
        end
        
        # Create upper interpolation
        itp_upper = LinearInterpolation(a_list_high, log10.(out_high), extrapolation_bc=Interpolations.Line())
        
        # Handle lower interpolation
        if !run_low
            a_list_low = [amin_guess * 0.99, amin_guess]
            out_low = [1e-100, 1e-100]
        elseif length(a_list_low) == 1.0
            a_list_low = [a_list_low[1] * 0.99, a_list_low[1]]
            out_low = [out_low[1], out_low[1]]
        end
        
        # Create lower interpolation
        itp_lower = LinearInterpolation(a_list_low, log10.(out_low), extrapolation_bc=Interpolations.Line())
        
        # Define boundaries
        abndry_upper = a_list_high[1]
        abndry_lower = a_list_low[end]
        
        # Special case for (3,2,2) which has a different buffer region behavior
        is_special_case = (n == 3 && l == 2 && m == 2)
        
        # Create an interpolation function with closure over the variables we need
        function create_interp_function(itp_u, itp_l, a_u, a_l, special)
            function interp_func(aspin)
                if aspin >= a_u
                    return 10.0 ^ itp_u(aspin)
                elseif aspin <= a_l
                    return -10.0 ^ itp_l(aspin)
                else
                    # Buffer region 
                    if special
                        return (-10.0 ^ itp_l(aspin) + 10.0 ^ itp_u(aspin)) / 2.0
                    else
                        return 0.0
                    end
                end
            end
            return interp_func
        end
        
        # Create interpolation function with proper closures
        interp_func = create_interp_function(itp_upper, itp_lower, abndry_upper, abndry_lower, is_special_case)
        
        # Store this function for later use
        push!(interp_functions, interp_func)
        
        # Calculate and store the rate for the current aBH
        SR_rates[idx] = interp_func(aBH)
    end
    
    # Create a Dict to access interpolation functions by quantum numbers
    interp_dict = Dict()
    for (i, (n, l, m, _)) in enumerate(qtm_cfigs)
        interp_dict[(n, l, m)] = interp_functions[i]
    end
    
    return SR_rates, interp_functions, interp_dict
end
