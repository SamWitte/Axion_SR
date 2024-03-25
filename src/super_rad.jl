using Random
using OrdinaryDiffEq
using Statistics
using Distributions
using DelimitedFiles
# using Dierckx
include("Constants.jl")
include("solve_sr_matching.jl")



function super_rad_check(M_BH, aBH, massB, f_a; spin=0, tau_max=1e4, alpha_max_cut=0.2, debug=false, solve_322=true, solve_n4=false, impose_low_cut=0.01, input_data="Masha")
   
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
        if (ergL(2, 1, 1, massB, M_BH) .>= OmegaH)&&(f_a .< 2.6e17)
            return aBH
        end
    end
    
    final_spin = solve_system(massB, f_a, aBH, M_BH, tau_max, debug=debug, solve_322=solve_322, impose_low_cut=impose_low_cut, input_data=input_data, solve_n4=solve_n4)
    # print("Spin diff.. \t ", aBH, "\t", final_spin, "\t", alph, "\n")
    return final_spin
    
end

function emax_211(MBH, mu, aBH)
    alph = GNew .* MBH .* mu
    
    emax_N = 1 .- 8 .* alph.^2 .+ 8 .* alph.^3 .* aBH .- sqrt.(abs.(1 .- 16 .* alph.^2 .+ 32 .* aBH .* alph.^3 - 16 .* aBH.^2 .* alph.^4)) # abs just in case...
    emax_D = 8 .* (-alph.^3 .+ aBH .* alph.^4)
    return (emax_N ./ emax_D)
end
    
function solve_system(mu, fa, aBH, M_BH, t_max; n_times=100, debug=true, solve_322=true, impose_low_cut=0.01, return_all_info=false, input_data="Masha", solve_n4=false)
    e_init = 1.0 ./ (GNew .* M_BH.^2 .* M_to_eV)
    if !solve_n4
        y0 = [e_init, e_init, aBH, M_BH]
    else
        y0 = [e_init, e_init, e_init, aBH, M_BH]
    end
    wait = 0


    u1_eq = false
    u1_fix = nothing
    u2_eq = false
    u2_fix = nothing
    u3_fix = nothing
    u4_fix = nothing
    u2_kill = false

    
    Emax2 = 1.0
    OmegaH = aBH ./ (2 .* (GNew .* M_BH) .* (1 .+ sqrt.(1 .- aBH.^2)))
    if (OmegaH .> ergL(2, 1, 1, mu, M_BH))
        Emax2 = emax_211(M_BH, mu, aBH)
    end
    
    
    Mvars = [mu, fa, Emax2, solve_322, aBH, M_BH, impose_low_cut]
    tspan = (0.0, t_max)
    saveat = (tspan[2] .- tspan[1]) ./ n_times
    
    SR211 = find_im_part(mu, M_BH, aBH, 2, 1, 1, Ntot=400) ./ (GNew * M_BH)
    SR322 = find_im_part(mu, M_BH, aBH, 3, 2, 2, Ntot=800) ./ (GNew * M_BH)
    a_prev = aBH
    
    function RHS_ax!(du, u, Mvars, t)
    
        # [e211, e322, aBH, MBH] or [e211, e322, e411, aBH, MBH]
        # mu [mass eV], fa[1/GeV]
        mu, fa, Emax2, solve_322, aBH_i, M_BH_i, impose_low_cut  = Mvars
        
        if !solve_n4
            spinI = 3
            massI = 4
        else
            spinI = 4
            massI = 5
        end
        
        alph = GNew .* u[massI] .* mu #
        rP = nothing
        if u[spinI] .> 0.998
            rP = 1.0 .+ sqrt.(1 - 0.998 .^2)
            u[spinI] = 0.998
        elseif u[spinI] .< 0.0
            rP = 2.0
            u[spinI] = 0.0
        else
            rP = 1.0 .+ sqrt.(1 - u[spinI].^2)
        end
        
        if u[1] < 0
            u[1] = 0.0
        end
        if u[2] < 0
            u[2] = 0.0
        end
        if solve_n4
            if u[3] < 0
                u[3] = 0.0
            end
        end


        if abs.(u[3] - a_prev) > 0.001            
            SR211 = find_im_part(mu, u[massI], u[spinI], 2, 1, 1, Ntot=400) ./ (GNew * u[spinI])
            SR322 = find_im_part(mu, u[massI], u[spinI], 3, 2, 2, Ntot=800) ./ (GNew * u[spinI])
       
            if solve_n4
                SR411 = find_im_part(mu, u[massI], u[spinI], 4, 1, 1, Ntot=400) ./ (GNew * u[spinI])
            end
            a_prev = u[3]
        end

        
        if u[1] .> Emax2
            SR211 *= 0.0
            kSR_211 *= 0.0
        end
        
        if input_data == "Doddy"
            k322BH = 0.0
            k2I_333 = 0.0
            kGW_22 = 0.0
            kGW_3t2 = 0.0
            kI_222 = 0.0
            kSR_322 = 0.0
            kGW_33 = 0.0
        else
            k322BH = 4e-7  # k^322xBH_211x211
            k2I_333 = 1e-8 # k^211xInfinity_322x322
            kGW_22 = 1e-2 # k^GW_211x211
            kGW_3t2 = 5e-6 # k^GW_{322->211}
            kI_222 = 1.5e-8 # k^Infinity_(211)^3
            kSR_322 = 8e-5
            kGW_33 = 3e-8
        end
        
        kGW_23 = 0.0 # k^GW_211x322
        kI_223 = 0.0 # ??? k^Infinity_{(211)^2 x 322}
        kI_233 = 0.0 # ??? k^Infinity_{(211)x 322^2}
        kI_333 = 0.0 # ???
        
        OmegaH = u[spinI] ./ (2 .* (GNew .* u[massI]) .* (1 .+ sqrt.(1 .- u[spinI].^2)))
        if (2 .* OmegaH .< ergL(3, 2, 2, mu, u[massI]))
            kSR_322 *= 0.0
            SR322 *= 0.0
        end
        if isnothing(u1_fix)
        #     du[1] = kSR_211 .* alph.^8 .* (u[spinI] .- 2 .* alph .* rP) .* u[1]
            du[1] = SR211 .* u[1] ./ mu
            du[1] += - 2 .* k322BH .* alph.^11 .* (M_pl ./ fa).^4 .* rP .* u[1].^2 .* u[2]
            du[1] += k2I_333 .* alph.^8 .* (M_pl ./ fa).^4 .* u[2].^2 .* u[1]
            du[1] += -2 .* kGW_22 .* alph.^14 .* u[1].^2
            du[1] += -3 .* kI_222 .* alph.^21 .* (M_pl ./ fa).^4 .* u[1].^3
            du[1] += kGW_3t2 .* alph.^10 .* u[1] .* u[2]
    ######  not being used
    #     du[1] += -kGW_23 .* alph.^16 .* u[1] .* u[2]
    #     du[1] += -2 .* kI_223 .* alph.^23 .* (M_pl ./ fa).^4 .* u[1].^2 .* u[2]
    #     du[1] += kI_233 .* alph.^25 .* (M_pl ./ fa).^4 .* u[1] .* u[2].^2
        else
            du[1] = 0.0
        end
        

        
        if isnothing(u2_fix)
            # du[2] = kSR_322 .* alph.^12 .* (u[spinI] .- alph .* rP) .* u[2]
            du[2] = SR322 .* u[2] ./ mu
            du[2] += k322BH .* alph.^11 .* (M_pl ./ fa).^4 .* rP .* u[1].^2 .* u[2]
            du[2] += -2 .* k2I_333 .* alph.^8 .* (M_pl ./ fa).^4 .* u[2].^2 .* u[1]
            du[2] += -2 .* kGW_33 .* alph.^18 .* u[2].^2
            du[2] += -kGW_23 .* alph.^16 .* u[1] .* u[2] .- kGW_3t2 .* alph.^10 .* u[1] .* u[2]
            # du[2] += -3 .* kI_333 .* alph.^27 .* (M_pl ./ fa).^4 .* u[2].^3
            # du[2] += -kI_223 .* alph.^23 .* (M_pl ./ fa).^4 .* u[1].^2 .* u[2]
            # du[2] += -2 .* kI_233 .* alph.^25 .* (M_pl ./ fa).^4 .* u[1] .* u[2].^2
            if solve_322 == false
                du[2] *= 0.0
            end
        else
            du[2] = 0.0
        end
        
        if solve_n4
            if (OmegaH .< ergL(4, 1, 1, mu, u[massI]))
                SR411 *= 0.0
            end
            du[3] = SR411 .* u[3] ./ mu
            
            du[3] += -2.5e-8 .* alph.^11 .* (M_pl ./ fa).^4 .* rP .* u[1] .* u[2] .* u[3] # 211x411^{322 x BH}
            du[1] += -2.5e-8 .* alph.^11 .* (M_pl ./ fa).^4 .* rP .* u[1] .* u[2] .* u[3]
            du[2] += 2.5e-8 .* alph.^11 .* (M_pl ./ fa).^4 .* rP .* u[1] .* u[2] .* u[3]
            
            du[3] += -2 .* 9.8e-11 .* alph.^11 .* (M_pl ./ fa).^4 .* rP .* u[2] .* u[3].^2 # 411x411^{322 x BH}
            du[2] += 9.8e-11 .* alph.^11 .* (M_pl ./ fa).^4 .* rP .* u[2] .* u[3].^2
            
            # to infinity
            du[3] += -3.8e-9 .* alph.^8 .* (M_pl ./ fa).^4 .* u[1] .* u[2] .* u[3] # 322x411^{211 x Inf}
            du[2] += -3.8e-9 .* alph.^8 .* (M_pl ./ fa).^4 .* u[1] .* u[2] .* u[3]
            du[1] += 3.8e-9 .* alph.^8 .* (M_pl ./ fa).^4 .* u[1] .* u[2] .* u[3]
            # print(SR411 .* u[3] ./ mu ./ (-3.8e-9 .* alph.^8 .* (M_pl ./ fa).^4 .* u[1] .* u[2] .* u[3]), "\n\n")
          
        end
        
    
        if !solve_n4
            du[3] = - SR211 .* u[1] ./ mu .- 2 .* SR322 .* u[2] ./ mu
            du[4] = - SR211 .* u[1] ./ mu .- SR322 .* u[2] ./ mu
            du[4] += k322BH .* alph.^11 .* (M_pl ./ fa).^4 .* rP .* u[1].^2 .* u[2]
            
            
            du[1] *= mu ./ hbar .* 3.15e7
            du[2] *= mu ./ hbar .* 3.15e7
            du[3] *= mu ./ hbar .* 3.15e7
            du[4] *= mu.^2 .* (GNew .* u[4].^2) ./ hbar .* 3.15e7
        else
            du[spinI] = - SR211 .* u[1] ./ mu .- 2 .* SR322 .* u[2] ./ mu .- SR411 .* u[3] ./ mu
            du[massI] = - SR211 .* u[1] ./ mu .-  SR322 .* u[2] ./ mu .- SR411 .* u[3] ./ mu
            du[massI] += k322BH .* alph.^11 .* (M_pl ./ fa).^4 .* rP .* u[1].^2 .* u[2]
            du[massI] += 9.8e-11 .* alph.^11 .* (M_pl ./ fa).^4 .* rP .* u[2] .* u[3].^2

            
            du[1] *= mu ./ hbar .* 3.15e7
            du[2] *= mu ./ hbar .* 3.15e7
            du[3] *= mu ./ hbar .* 3.15e7
            
            du[spinI] *= mu ./ hbar .* 3.15e7
            du[massI] *= mu.^2 .* (GNew .* u[4].^2) ./ hbar .* 3.15e7
        end
        
        
        # print("Test \t", u[1] ./ du[1], "\t", u[2] ./ du[2], "\n")
        return
    end


    function check_timescale(u, t, integrator)
        
        
        if !solve_n4
            spinI = 3
            massI = 4
        else
            spinI = 4
            massI = 5
        end
        
        alph = GNew .* u[massI] .* mu
        rP = nothing
        if u[spinI] .> 0.998
            rP = 1.0 .+ sqrt.(1 - 0.998 .^2)
        elseif u[spinI] .< 0.0
            rP = 2.0
            u[spinI] = 0.0
        else
            rP = 1.0 .+ sqrt.(1 - u[spinI].^2)
        end
        du = get_du(integrator)
        
      
        t1 = abs.(u[1] ./ du[1])
        t2 = abs.(u[2] ./ du[2])
        if solve_n4
            t5 = abs.(u[3] ./ du[3])
        end
        
        if (u1_eq)||(u2_eq)      
            t1 = 1e100
            integrator.u[1] = u1_fix
            t2 = 1e100
            integrator.u[2] = u2_fix
            if u2_eq && solve_n4
                t5 = 1e100
                integrator.u[3] = u3_fix

            end
            
            du[spinI] = - SR211 .* u[1] ./ mu .-  2 .* SR322 .* u[2] ./ mu
            if solve_n4
                SR411 = sr_rates(4, 1, 1, mu, u[massI], u[spinI], impose_low_cut=impose_low_cut, solve_322=solve_322)
                du[spinI] += - SR411 .* u[3] ./ mu
            end
        end
        t3 = abs.(u[spinI] ./ du[spinI])
        t4 = abs.(u[massI] ./ du[massI])
        
       
        
        if !solve_n4
            tcheck = minimum([t1 t2 t3])
        else
            tcheck = minimum([t1 t2 t3 t5])
        end
        
        wait += 1
        
        
        if debug && (wait%10==0)
            print(t, "\t", u[1], "\t", u[2], "\t", u[3], "\t", u[4], "\n")
        end
        
        if (tcheck .>= 100.0 .* integrator.dt) && (wait % 10 == 0)
            return true
        elseif (tcheck .<= integrator.dt)
            return true
        elseif (integrator.dt .<= 1e-7)
            return true
        else
            return false
        end
    end
    function affect_timescale!(integrator)
        
        if !solve_n4
            spinI = 3
            massI = 4
        else
            spinI = 4
            massI = 5
        end
        du = get_du(integrator)

        t1 = abs.(integrator.u[1] ./ du[1])
        t2 = abs.(integrator.u[2] ./ du[2])
        if solve_n4
            t5 = abs.(integrator.u[3] ./ du[3])
        end

        if u1_eq || u2_eq
            
            
            t1 = 1e100
            integrator.u[1] = u1_fix
            t2 = 1e100
            integrator.u[2] = u2_fix
            if u2_eq
                t5 = 1e100
                integrator.u[3] = u3_fix
            end
           
            du[spinI] = - SR211 .* integrator.u[1] ./ mu .-  2 .* SR322 .* integrator.u[2] ./ mu
            if solve_n4
                du[spinI] += - SR411 .* integrator.u[3] ./ mu
            end
        end
        
        t3 = abs.(integrator.u[spinI] ./ du[spinI])
        if !solve_n4
            tcheck = minimum([t1 t2 t3])
        else
            tcheck = minimum([t1 t2 t3 t5])
        end

        
        if (tcheck .<= integrator.dt)
            set_proposed_dt!(integrator, tcheck .* 0.1)
        elseif (integrator.dt .<= 1e-7)
            terminate!(integrator)
        elseif (wait % 5 == 0)
            set_proposed_dt!(integrator, integrator.dt .* 1.05)
        end
    end
    
    function check_eq(u, t, integrator)
        
        if !solve_n4
            spinI = 3
            massI = 4
        else
            spinI = 4
            massI = 5
        end
        
        alph = GNew .* u[massI] .* mu
        du = get_du(integrator)
        eq_threshold = 1e-4
        
        
        if u2_kill
            terminate!(integrator)
        end
        if u2_eq
            return true
        end

        if u1_eq
            # check if 411 can change equilibrium!
            if !solve_n4
                return true
            else
                rP = 1.0 .+ sqrt.(1 - u[spinI].^2)
                hold_n1 = (-2.5e-8 .* alph.^11 .* (M_pl ./ fa).^4 .* rP .* u[1] .* u[2] .* u[3] + 3.8e-9 .* (alph.^8) .* (M_pl ./ fa).^4 .* u[1] .* u[2] .* u[3]) .* integrator.dt
                if abs.(hold_n1 ./ u[1]) .> eq_threshold
                    u1_eq = false
                    return false
                else
                    return true
                end
            end
        end
        
        
        t3 = abs.(integrator.u[spinI] ./ du[spinI])

        k322BH = 4e-7  # k^322xBH_211x211
        GR_322 = SR322 .* integrator.u[2] .+ k322BH .* alph.^11 .* (M_pl ./ fa).^4 .* rP .* integrator.u[1].^2 .* integrator.u[2] .* mu
        if solve_n4
            GR_411 = SR411 .* integrator.u[3]
            stable411 = abs.(du[3] ./ (GR_411 ./ hbar .* 3.15e7))
        end
        
        
        stable211 = abs.(du[1] ./ (SR211 .* integrator.u[1] ./ hbar .* 3.15e7))
        stable322 = abs.(du[2] ./ (GR_322 ./ hbar .* 3.15e7))

        
        if (stable211 .< eq_threshold)&&(stable322 .< eq_threshold)
            u1_eq = true
            u1_fix = integrator.u[1]
            u2_fix = integrator.u[2]
            u4_fix = integrator.u[massI]
            if solve_n4
                if (stable411 .< eq_threshold)
                    u2_eq = true
                    u3_fix = integrator.u[3]
                end
            end
        else
            u1_eq = false
        end
        
        
        if u1_eq||u2_eq
            return true
        else
            return false
        end
    end
    function affect_eq!(integrator)
        
        if !solve_n4
            spinI = 3
            massI = 4
        else
            spinI = 4
            massI = 5
        end
        du = get_du(integrator)
        
        if u1_eq && u2_eq
            du[1] = 0.0
            du[2] = 0.0
            du[massI] = 0.0
            if !solve_n4
                set_u!(integrator, [u1_fix, u2_fix, integrator.u[spinI], u4_fix])
            else
                set_u!(integrator, [u1_fix, u2_fix, u3_fix, integrator.u[spinI], u4_fix])
            end
        elseif u1_eq
            du[1] = 0.0
            du[2] = 0.0
            du[massI] = 0.0
            if !solve_n4
                set_u!(integrator, [u1_fix, u2_fix, integrator.u[spinI], u4_fix])
            else
                set_u!(integrator, [u1_fix, u2_fix, integrator.u[3], integrator.u[spinI], u4_fix])
            end
        end
        
    end
    

    function check_spin(u, t, integrator)
        if !solve_n4
            BHs_idx = 3
        else
            BHs_idx = 4
        end
        if u[BHs_idx] .> (aBH .+ 0.01)
            return true
        elseif u[BHs_idx] .<= 0.0
            return true
        else
            return false
        end
    end
    function affect_spin!(integrator)
        if !solve_n4
            BHs_idx = 3
        else
            BHs_idx = 4
        end
        if integrator.u[BHs_idx] .> aBH
            integrator.u[BHs_idx] = aBH
        elseif integrator.u[BHs_idx] .< 0.0
            integrator.u[BHs_idx] = 0.0
        end
        set_proposed_dt!(integrator, integrator.dt .* 0.3)
    end
    
    
        

    cbackdt = DiscreteCallback(check_timescale, affect_timescale!, save_positions=(false, true))
    cback_equil = DiscreteCallback(check_eq, affect_eq!, save_positions=(false, true))
    cbackspin = DiscreteCallback(check_spin, affect_spin!, save_positions=(false, true))
    
    cbset = CallbackSet(cback_equil, cbackdt, cbackspin)
    
    prob = ODEProblem(RHS_ax!, y0, tspan, Mvars, reltol=1e-6, abstol=1e-6)
    # prob = ODEProblem(RHS_ax!, y0, tspan, Mvars, reltol=1e-6, abstol=1e-10)
    
    rP = 1.0 .+ sqrt.(1 - aBH.^2)
    if (aBH .- 2 .* (GNew .* M_BH .* mu) .* rP) .> 0.0
        dt_guess = abs.(4e-2 .* (GNew .* M_BH .* mu ).^8 .* (aBH .- 2 .* (GNew .* M_BH .* mu) .* rP) .* mu / hbar .* 3.15e7).^(-1)
    else
        dt_guess = abs.(8e-5 .* (GNew .* M_BH .* mu ).^12 .* (aBH .- (GNew .* M_BH .* mu ) .* rP) .* mu / hbar .* 3.15e7).^(-1)
    end
    if debug
        print("Time guess \t", dt_guess, "\n")
    end
    
    sol = solve(prob, Euler(), dt=dt_guess, saveat=saveat, callback=cbset)
    # sol = solve(prob, Vern6(), saveat=saveat, callback=cbset)
    
    
    state211 = [sol.u[i][1] for i in 1:length(sol.u)]
    state322 = [sol.u[i][2] for i in 1:length(sol.u)]
    if !solve_n4
        spinBH = [sol.u[i][3] for i in 1:length(sol.u)]
        MassB = [sol.u[i][4] for i in 1:length(sol.u)]
    else
        state411 = [sol.u[i][3] for i in 1:length(sol.u)]
        spinBH = [sol.u[i][4] for i in 1:length(sol.u)]
        MassB = [sol.u[i][5] for i in 1:length(sol.u)]
    end
    
    # alph = GNew .* M_BH .* mu
    # bose_thresh_e = 1024 .* pi * (fa ./ M_pl).^2 ./ ( 9 .* alph.^3 )
    # print("Condition two \t", state211[end] ./ bose_thresh_e, "\n")
    
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
            return sol.t, state211, state322, state411, spinBH, MassB
        end
    end

    if isnan(spinBH[end])
        spinBH = spinBH[.!isnan.(spinBH)]
    end
    return spinBH[end]
 
end
