using Random
using OrdinaryDiffEq
using Statistics
using Distributions
using DelimitedFiles
using Interpolations
include("Constants.jl")
include("solve_sr_rates.jl")
include("load_rates.jl")


function super_rad_check(M_BH, aBH, massB, f_a; spin=0, tau_max=1e4, alpha_max_cut=10.0, debug=false, solve_322=true, solve_n4=false, impose_low_cut=0.01, input_data="Masha", stop_on_a=0, eq_threshold=1e-100, abstol=1e-30)
   
    alph = GNew .* M_BH .* massB #
    if debug
        print("Alpha \t", alph, "\n")
    end
    
    if alph .> alpha_max_cut
        if debug
            print("Need higher-level system... \n")
        end
        return aBH
    elseif alph .< impose_low_cut
        return aBH
    end
    
    
    if input_data == "Masha"
        OmegaH = aBH ./ (2 .* (GNew .* M_BH) .* (1 .+ sqrt.(1 .- aBH.^2)))
        if (ergL(2, 1, 1, massB, M_BH, aBH) .>= OmegaH)&&(f_a .< 2.6e17)
            return aBH
        end
    end
    
    final_spin, final_BH = solve_system(massB, f_a, aBH, M_BH, tau_max, debug=debug, solve_322=solve_322, impose_low_cut=impose_low_cut, input_data=input_data, solve_n4=solve_n4, stop_on_a=stop_on_a, eq_threshold=eq_threshold, abstol=abstol)
    # print("Spin diff.. \t ", aBH, "\t", final_spin, "\t", alph, "\n")
    return final_spin, final_BH
    
end

function emax_211(MBH, mu, aBH)
    alph = GNew .* MBH .* mu
    
    emax_N = 1 .- 8 .* alph.^2 .+ 8 .* alph.^3 .* aBH .- sqrt.(abs.(1 .- 16 .* alph.^2 .+ 32 .* aBH .* alph.^3 - 16 .* aBH.^2 .* alph.^4)) # abs just in case...
    emax_D = 8 .* (-alph.^3 .+ aBH .* alph.^4)
    return (emax_N ./ emax_D)
end
    
function solve_system(mu, fa, aBH, M_BH, t_max; n_times=10000, debug=true, solve_322=true, impose_low_cut=0.01, return_all_info=false, input_data="Masha", solve_n4=false, eq_threshold=1e-4, stop_on_a=0, abstol=1e-30, non_rel=true)
    e_init = 1.0 ./ (GNew .* M_BH.^2 .* M_to_eV) # unitless
    if !solve_n4
        y0 = [e_init, e_init, aBH, M_BH]
    else
        y0 = [e_init, e_init, e_init, e_init, e_init, aBH, M_BH]
    end
    wait = 0

    alph = GNew .* M_BH .* mu
    
    u1_eq = false
    u1_fix = nothing
    u2_eq = false
    u2_fix = nothing
    u3_fix = nothing
    u4_fix = nothing
    u5_fix = nothing
    u6_fix = nothing
    u2_kill = false
    
    alph = GNew .* M_BH .* mu
    if input_data != "Doddy"
        e2_maxBN = 1024 * pi * (fa / M_pl).^2 ./ (9 * alph.^3)
        # print("Max E211 \t", e2_maxBN, "\n")
    else
        e2_maxBN = 1024 * pi * (fa / M_pl).^2 ./ (9 * alph.^3)
    end
    e3_maxBN = e2_maxBN .* (3 ./ 2).^4
    e4_maxBN = e2_maxBN .* (4 ./ 2).^4
    
    
    
    Emax2 = 1.0
    OmegaH = aBH ./ (2 .* (GNew .* M_BH) .* (1 .+ sqrt.(1 .- aBH.^2)))
    if (OmegaH .> ergL(2, 1, 1, mu, M_BH, aBH))
        Emax2 = emax_211(M_BH, mu, aBH)
    end
    
    
    Mvars = [mu, fa, Emax2, solve_322, aBH, M_BH, impose_low_cut]
    tspan = (0.0, t_max)
    saveat = (tspan[2] .- tspan[1]) ./ n_times
    
    
    
    N_pts_interp = 5
    xtol_slv = 1e-15
    iter_slv = 50
    
    if !solve_n4
        spinI = 3
        massI = 4
    else
        spinI = 6
        massI = 7
    end
    idx_lvl = spinI - 1
    
    SR_rates = zeros(idx_lvl)
    m_list = [1 2 1 2 3] # [211, 322, 411, 422, 433]
    bn_list = [e2_maxBN e3_maxBN e4_maxBN e4_maxBN e4_maxBN]
    
    
    n = 2
    l = 1
    m = 1
    amin_guess = 8 * m * n.^2 .* alph .* (2 .* n.^2 .+ alph.^2) ./ (16 .* n.^4 .* alph.^2 .+ m.^2 .* (2 .* n.^2 .+ alph.^2).^2)
    alist, pts211 = compute_gridded(mu, M_BH, aBH, n, l, m; iter=iter_slv, xtol=xtol_slv, npts=N_pts_interp, amin=(amin_guess .* 0.95))
    itp_211 = LinearInterpolation(alist, log10.(pts211), extrapolation_bc=Line())
    SR211 = 10 .^itp_211(aBH)
    SR_rates[1] = SR211

    
    n = 3
    l = 2
    m = 2
    amin_guess = 8 * m * n.^2 .* alph .* (2 .* n.^2 .+ alph.^2) ./ (16 .* n.^4 .* alph.^2 .+ m.^2 .* (2 .* n.^2 .+ alph.^2).^2)
    alist, pts322 = compute_gridded(mu, M_BH, aBH, n, l, m; iter=iter_slv, xtol=xtol_slv, npts=N_pts_interp, amin=(amin_guess .* 0.95))
    itp_322 = LinearInterpolation(alist, log10.(pts322), extrapolation_bc=Line())
    SR322 = 10 .^itp_322(aBH)
    SR_rates[2] = SR322
    
    if solve_n4
        n = 4
        l = 1
        m = 1
        amin_guess = 8 * m * n.^2 .* alph .* (2 .* n.^2 .+ alph.^2) ./ (16 .* n.^4 .* alph.^2 .+ m.^2 .* (2 .* n.^2 .+ alph.^2).^2)
        alist, pts411 = compute_gridded(mu, M_BH, aBH, n, l, m; iter=iter_slv, xtol=xtol_slv, npts=N_pts_interp, amin=(amin_guess .* 0.95))
        itp_411 = LinearInterpolation(alist, log10.(pts411), extrapolation_bc=Line())
        SR411 = 10 .^itp_411(aBH)
        SR_rates[3] = SR411
        
        n = 4
        l = 2
        m = 2
        amin_guess = 8 * m * n.^2 .* alph .* (2 .* n.^2 .+ alph.^2) ./ (16 .* n.^4 .* alph.^2 .+ m.^2 .* (2 .* n.^2 .+ alph.^2).^2)
        alist, pts422 = compute_gridded(mu, M_BH, aBH, n, l, m; iter=iter_slv, xtol=xtol_slv, npts=N_pts_interp, amin=(amin_guess .* 0.95))
        itp_422 = LinearInterpolation(alist, log10.(pts422), extrapolation_bc=Line())
        SR422 = 10 .^itp_422(aBH)
        SR_rates[4] = SR422
        
        n = 4
        l = 3
        m = 3
        amin_guess = 8 * m * n.^2 .* alph .* (2 .* n.^2 .+ alph.^2) ./ (16 .* n.^4 .* alph.^2 .+ m.^2 .* (2 .* n.^2 .+ alph.^2).^2)
        alist, pts433 = compute_gridded(mu, M_BH, aBH, n, l, m; iter=iter_slv, xtol=xtol_slv, npts=N_pts_interp, amin=(amin_guess .* 0.95))
        itp_433 = LinearInterpolation(alist, log10.(pts433), extrapolation_bc=Line())
        SR433 = 10 .^itp_433(aBH)
        SR_rates[5] = SR433
    end
    
    
    rates = load_rate_coeffs(mu, M_BH, aBH, fa; non_rel=true, input_data=input_data, solve_n4=solve_n4)
    
    
    du_eq = zeros(length(y0))
    equilib = false
    post_eq = 0.0

    #### DEFINING FUNCTION FOR EVOLUTION ######
    function RHS_ax!(du, u, Mvars, t)
    
        # [e211, e322, aBH, MBH] or [e211, e322, e411, aBH, MBH]
        # mu [mass eV], fa[1/GeV]
        
        
        mu, fa, Emax2, solve_322, aBH_i, M_BH_i, impose_low_cut  = Mvars
            
        if equilib
            du[:] .= du_eq
        else
            
            # alph = GNew .* u[massI] .* mu #
            rP = nothing
            if u[spinI] .> maxSpin
                rP = 1.0 .+ sqrt.(1 - maxSpin .^2)
                u[spinI] = maxSpin
            elseif u[spinI] .< 0.0
                rP = 2.0
                u[spinI] = 0.0
            else
                rP = 1.0 .+ sqrt.(1 - u[spinI].^2)
            end
            
            for i in 1:idx_lvl
                if u[i] < e_init
                    u[i] = e_init
                end
            end
            
            OmegaH = u[spinI] ./ (2 .* (GNew .* u[massI]) .* (1 .+ sqrt.(1 .- u[spinI].^2)))
            
            
            if (OmegaH .< ergL(2, 1, 1, mu, u[massI], u[spinI]))
                SR211 = 0.0
            else
                SR211 = 10 .^itp_211(u[spinI])
            end
            SR_rates = [SR211]
            if (2 .* OmegaH .< ergL(3, 2, 2, mu, u[massI], u[spinI]))
                SR322 = 0.0
            else
                SR322 = 10 .^itp_322(u[spinI])
            end
            append!(SR_rates, SR322)
            if solve_n4
                if (u[3] >= e_init)&&(OmegaH .> ergL(4, 1, 1, mu, u[massI], u[spinI]))
                    SR411 = 10 .^itp_411(u[spinI])
                else
                    SR411 = 0.0
                end
                append!(SR_rates, SR411)
                if (u[4] >= e_init)&&(OmegaH * 2 .> ergL(4, 2, 2, mu, u[massI], u[spinI]))
                    SR422 = 10 .^itp_422(u[spinI])
                else
                    SR422 = 0.0
                end
                append!(SR_rates, SR422)
                if (u[5] >= 1e_init)&&(OmegaH * 3 .> ergL(4, 3, 3, mu, u[massI], u[spinI]))
                    SR433 = 10 .^itp_433(u[spinI])
                else
                    SR433 = 0.0
                end
                append!(SR_rates, SR433)
            end
            
            
            if u[1] .> Emax2
                SR_rates[1] *= 0.0
            end
            
        
            # SR terms
            du[spinI] = 0.0
            du[massI] = 0.0
            for i in 1:idx_lvl
                du[i] = SR_rates[i] .* u[i] ./ mu
                du[spinI] += - m_list[i] * SR_rates[i] .* u[i] ./ mu
            end
            
            # Scattering terms
            rate_keys = collect(keys(rates))
            for i in 1:length(rate_keys)
                idxV, sgn = key_to_indx(rate_keys[i]; solve_n4=solve_n4)
                u_term_tot = 1.0
                
                for j in 1:length(sgn)
                    if (idxV[j] <= idx_lvl)&&(idxV[j] > 0)
                        u_term_tot *= u[idxV[j]]
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
            
            if solve_322 == false
                du[2] = 0.0
            end
            
            # check bosenova and correct units
            for i in 1:idx_lvl
                if (u[i] > bn_list[i])&&(du[i] > 0)
                    du[i] = 0.0
                else
                    du[i] *= mu ./ hbar .* 3.15e7
                end
                if (u[i] <= e_init)&&(du[i] < 0)
                    du[i] = 0.0
                end
            end

            du[spinI] *= mu ./ hbar .* 3.15e7
            du[massI] *= (mu .* u[massI]) .* (mu .* GNew .* u[massI]) ./ hbar .* 3.15e7
            
        
        end
        return
    end
    
    
#    isapproxsigfigs(a, b, precision) = round(a, sigdigits=precision) == round(b, sigdigits=precision)
    ## CALLBACK 1
    function check_eq(u, t, integrator)
        
        du = get_du(integrator)
        tlist = []
        idxlist = []
        thresh = 3 # sig figs for check
        
        if !equilib
            for i in 1:idx_lvl
                if u[i] > e_init
                    append!(tlist, (u[i] ./ du[i]))
                    append!(idxlist, i)
                end
            end
            
            append!(tlist, (u[spinI] ./ du[spinI]))
            minTime = minimum(abs.(tlist))
            
            n_step = minTime ./ integrator.dt
            if n_step > 1e7
                post_eq = integrator.t + minTime / 10.0
                du_eq = copy(du)
                for i in 1:idx_lvl
                    du_eq[i] = 0.0
                end
                return true
            else
                return false
            end
        
        else
            if integrator.t > post_eq
                print("No more Eq!!! \n\n")
                equilib = false
            end
        
            return false
        end
    end
    
    function affect_eq!(integrator)
        print("Eqilibrium!!! \n\n")
        equilib = true
    end
    
    threshold_step = 1e-1
    #### CALLBACK 2
    function check_timescale(u, t, integrator)
        du = get_du(integrator)
        tlist = []
        
        for i in 1:idx_lvl
            if u[i] > 1e-100
                append!(tlist, (u[i] ./ du[i]))
            end
        end
            
        append!(tlist, u[spinI] ./ du[spinI])
        tmin = minimum(abs.(tlist))
        # print(integrator.dt, "\t", tlist, "\n")
        
        if ((integrator.dt < threshold_step * tmin)&&(wait % 10 == 0))||(integrator.dt > threshold_step * tmin)
            return true
        else
            return false
        end
            
    end
    function affect_timescale!(integrator)
        du = get_du(integrator)
        tlist = []
        
        for i in 1:idx_lvl
            if integrator.u[i] > 1e-100
                append!(tlist, (integrator.u[i] ./ du[i]))
            end
        end
            
        append!(tlist, integrator.u[spinI] ./ du[spinI])
        tmin = minimum(abs.(tlist))
        if (integrator.dt < threshold_step * tmin)
            set_proposed_dt!(integrator, integrator.dt .* 1.05)
        else
            set_proposed_dt!(integrator, integrator.dt .* 0.5)
        end
    end
    

    function check_spin(u, t, integrator)
        if debug && (wait%10000==0)
            if solve_n4
                print(t, "\t", u[1], "\t", u[2], "\t", u[3], "\t", u[4], "\t", u[5], "\t", u[6], "\t", u[7], "\n")
            else
                print(t, "\t", u[1], "\t", u[2], "\t", u[3], "\t", u[4], "\n")
            end
        end
        if u[spinI] <= stop_on_a
            return true
        end
        if u[spinI] .> (aBH .+ 0.01)
            return true
        elseif u[spinI] .<= 0.0
            return true
        else
            return false
        end
    end
    function affect_spin!(integrator)

        if integrator.u[spinI] <= stop_on_a
            terminate!(integrator)
        end
        if integrator.u[spinI] .> aBH
            integrator.u[spinI] = aBH
        elseif integrator.u[spinI] .< 0.0
            integrator.u[spinI] = 0.0
        end
        set_proposed_dt!(integrator, integrator.dt .* 0.3)
    end
    

    
    
    rP = 1.0 .+ sqrt.(1 - aBH.^2)
    t1 = (SR211 ./ hbar .* 3.15e7).^(-1)
    t2 = (SR322 ./ hbar .* 3.15e7).^(-1)
    if solve_n4
        t3 = (SR433 ./ hbar .* 3.15e7).^(-1)
        tlist = [t1 t2 t3]
    else
        tlist = [t1 t2]
    end
    dt_guess = minimum(tlist)
    if debug
        print("Time guess \t", dt_guess, "\n")
    end

    cbackdt = DiscreteCallback(check_timescale, affect_timescale!, save_positions=(false, true))
    cback_equil = DiscreteCallback(check_eq, affect_eq!, save_positions=(false, true))
    cbackspin = DiscreteCallback(check_spin, affect_spin!, save_positions=(false, true))
    
    if solve_n4
    
        # cbset = CallbackSet(cbackdt, cbackspin)
        cbset = CallbackSet(cback_equil, cbackspin, cbackdt)
        
        
        
        prob = ODEProblem(RHS_ax!, y0, tspan, Mvars, reltol=1e-3, abstol=abstol)
        sol = solve(prob, Rosenbrock23(), dt=dt_guess, saveat=saveat, callback=cbset, maxiters=5e6)
        
        # sol = solve(prob, AutoTsit5(Rodas4()), dt=dt_guess, saveat=saveat, callback=cbset)
        # sol = solve(prob, KenCarp5(), dt=dt_guess, saveat=saveat, callback=cbset)
        # sol = solve(prob, Rodas4(), dt=dt_guess, saveat=saveat, callback=cbset)
        # sol = solve(prob, Euler(), dt=dt_guess, saveat=saveat, callback=cbset)
    else
        
        # abstol = 1e-30
        # abstol = 1e-15
        
        cbset = CallbackSet(cback_equil, cbackdt, cbackspin)
        # prob = ODEProblem(RHS_ax!, y0, tspan, Mvars, reltol=1e-5, abstol=1e-15)
        # sol = solve(prob, Rodas4(), dt=dt_guess, saveat=saveat, callback=cbset)
        
        
        # sol = solve(prob, Euler(), dt=dt_guess, saveat=saveat, callback=cbset)
        if input_data != "Doddy"
            prob = ODEProblem(RHS_ax!, y0, tspan, Mvars, reltol=1e-3, abstol=abstol, maxiters=5e5)
            sol = solve(prob, Rosenbrock23(), dt=dt_guess, saveat=saveat, callback=cbset)
            # sol = solve(prob, Euler(), dt=dt_guess, saveat=saveat, callback=cbset)
        else
            # prob = ODEProblem(RHS_ax!, y0, tspan, Mvars, reltol=1e-5, abstol=abstol, maxiters=5e7)
            prob = ODEProblem(RHS_ax!, y0, tspan, Mvars, reltol=1e-7, abstol=abstol, maxiters=5e6)
            # sol = solve(prob, Rodas4(), dt=dt_guess, saveat=saveat, callback=cbset)
            sol = solve(prob, Euler(), dt=dt_guess, saveat=saveat, callback=cbset)
        end
    end

    # sol = solve(prob, Euler(), dt=dt_guess, saveat=saveat, callback=cbset)
    
    state211 = [sol.u[i][1] for i in 1:length(sol.u)]
    state322 = [sol.u[i][2] for i in 1:length(sol.u)]
    if !solve_n4
        spinBH = [sol.u[i][3] for i in 1:length(sol.u)]
        MassB = [sol.u[i][4] for i in 1:length(sol.u)]
    else
        state411 = [sol.u[i][3] for i in 1:length(sol.u)]
        state422 = [sol.u[i][4] for i in 1:length(sol.u)]
        state433 = [sol.u[i][5] for i in 1:length(sol.u)]
       
        spinBH = [sol.u[i][6] for i in 1:length(sol.u)]
        MassB = [sol.u[i][7] for i in 1:length(sol.u)]
    end
    
    # alph = GNew .* M_BH .* mu
    # bose_thresh_e = 1024 .* pi * (fa ./ M_pl).^2 ./ ( 9 .* alph.^3 )
    # print("Condition two \t", state211[end] ./ bose_thresh_e, "\n")
    # alph_out = GNew .* 14.28 .* mu
    # maxa = 4 .* alph_out ./ (1 .+ 4 .* alph_out.^2)
    # print("CHECK \t", maxa, "\n\n")
    if debug
        if !solve_n4
            # print("Initial \t", state211[1], "\t", state322[1], "\t", spinBH[1], "\t", MassB[1], "\n")
            # print("Final \t", state211[end], "\t", state322[end], "\t", spinBH[end], "\t", MassB[end], "\n")
            print("Initial \t", state211[1], "\t", state322[1], "\t", spinBH[1], "\t", MassB[1], "\n")
            print("Final \t", state211[end], "\t", state322[end], "\t", spinBH[end], "\t", MassB[end], "\n")
        else
            print("Initial \t", state211[1], "\t", state322[1] , "\t",  state411[1] , "\t", spinBH[1], "\t", MassB[1], "\n")
            print("Final \t", state211[end] , "\t", state322[end] , "\t",  state411[1] , "\t", spinBH[end], "\t", MassB[end], "\n")
        end
    end

    if return_all_info
        if !solve_n4
            return sol.t, state211, state322, spinBH, MassB
        else
            return sol.t, state211, state322, state411, state422, state433, spinBH, MassB
        end
    end

    if isnan(spinBH[end])
        spinBH = spinBH[.!isnan.(spinBH)]
    end
    return spinBH[end], MassB[end]
 
end
