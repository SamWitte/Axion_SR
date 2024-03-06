using Random
using OrdinaryDiffEq
using Statistics
using Distributions
using DelimitedFiles
using Dierckx
include("Constants.jl")




function super_rad_check(M_BH, aBH, massB, f_a; spin=0, tau_max=1e4, alpha_max_cut=0.2, debug=false, solve_322=true, impose_low_cut=0.01, input_data="Masha")
   
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
    if input_data != "Doddy"
        elseif (alph .> 0.3)&&(f_a .< 2.6e17)
        return aBH
    end
    
    final_spin = solve_system(massB, f_a, aBH, M_BH, tau_max, debug=debug, solve_322=solve_322, impose_low_cut=impose_low_cut, input_data=input_data)
    # print("Spin diff.. \t ", aBH, "\t", final_spin, "\t", alph, "\n")
    return final_spin
    
end

function emax_211(MBH, mu, aBH)
    alph = GNew .* MBH .* mu
    
    emax_N = 1 .- 8 .* alph.^2 .+ 8 .* alph.^3 .* aBH .- sqrt.(abs.(1 .- 16 .* alph.^2 .+ 32 .* aBH .* alph.^3 - 16 .* aBH.^2 .* alph.^4)) # abs just in case...
    emax_D = 8 .* (-alph.^3 .+ aBH .* alph.^4)
    return (emax_N ./ emax_D)
end
    
function solve_system(mu, fa, aBH, M_BH, t_max; n_times=100, debug=true, solve_322=true, impose_low_cut=0.01, return_all_info=false, input_data="Masha")
    
    y0 = [1.0 ./ (GNew .* M_BH.^2 .* M_to_eV), 1.0 ./ (GNew .* M_BH.^2 .* M_to_eV), aBH, M_BH]
    wait = 0


    u1_eq = false
    u1_fix = nothing
    u2_eq = false
    u2_fix = nothing
    u4_fix = nothing
    u2_kill = false
    u2_rough = []
    u1_rough = []

    
    Emax2 = 1.0
    OmegaH = aBH ./ (2 .* (GNew .* M_BH) .* (1 .+ sqrt.(1 .- aBH.^2)))
    if (OmegaH .> ergL(2, 1, 1, mu, M_BH))
        Emax2 = emax_211(M_BH, mu, aBH)
    end
    
    
    Mvars = [mu, fa, Emax2, solve_322, aBH, M_BH, impose_low_cut]
    tspan = (0.0, t_max)
    saveat = (tspan[2] .- tspan[1]) ./ n_times
    
    function RHS_ax!(du, u, Mvars, t)
    
        # [e211, e322, aBH, MBH]
        # mu [mass eV], fa[1/GeV]
        
        mu, fa, Emax2, solve_322, aBH_i, M_BH_i, impose_low_cut  = Mvars
       
        alph = GNew .* u[4] .* mu #
        rP = nothing
        if u[3] .> 0.998
            rP = 1.0 .+ sqrt.(1 - 0.998 .^2)
            u[3] = 0.998
        elseif u[3] .< 0.0
            rP = 2.0
            u[3] = 0.0
        else
            rP = 1.0 .+ sqrt.(1 - u[3].^2)
        end
        
        if u[1] < 0
            u[1] = 0.0
        end
        if u[2] < 0
            u[2] = 0.0
        end

        
        kSR_211 = 4e-2 # kSR_211
        
        
        SR211 = sr_rates(2, 1, 1, mu, u[4], u[3], impose_low_cut=impose_low_cut, solve_322=solve_322)
        SR322 = sr_rates(3, 2, 2, mu, u[4], u[3], impose_low_cut=impose_low_cut, solve_322=solve_322)
        # print((SR211 ./ hbar .* 3.15e7).^(-1) , "\t", (SR322 ./ hbar .* 3.15e7).^(-1),"\n")
        
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
        
        OmegaH = u[3] ./ (2 .* (GNew .* u[4]) .* (1 .+ sqrt.(1 .- u[3].^2)))
        if (2 .* OmegaH .< ergL(3, 2, 2, mu, u[4]))
            kSR_322 *= 0.0
            SR322 *= 0.0
        end
        if isnothing(u1_fix)
        #     du[1] = kSR_211 .* alph.^8 .* (u[3] .- 2 .* alph .* rP) .* u[1]
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
            # du[2] = kSR_322 .* alph.^12 .* (u[3] .- alph .* rP) .* u[2]
            du[2] = SR322 .* u[2] ./ mu
            du[2] += k322BH .* alph.^11 .* (M_pl ./ fa).^4 .* rP .* u[1].^2 .* u[2]
            du[2] += -2 .* k2I_333 .* alph.^8 .* (M_pl ./ fa).^4 .* u[2].^2 .* u[1]
            du[2] += -2 .* kGW_33 .* alph.^18 .* u[2].^2
            du[2] += -kGW_23 .* alph.^16 .* u[1] .* u[2] .- kGW_3t2 .* alph.^10 .* u[1] .* u[2]
            du[2] += -3 .* kI_333 .* alph.^27 .* (M_pl ./ fa).^4 .* u[2].^3
            du[2] += -kI_223 .* alph.^23 .* (M_pl ./ fa).^4 .* u[1].^2 .* u[2]
            du[2] += -2 .* kI_233 .* alph.^25 .* (M_pl ./ fa).^4 .* u[1] .* u[2].^2
            if solve_322 == false
                du[2] *= 0.0
            end
        else
            du[2] = 0.0
        end
        
        
    
        # du[3] = - SR211 .* u[1] ./ mu .- 2 .* SR322 .* u[2] ./ mu
        du[3] = - SR211 .* u[1] ./ mu .- SR322 .* u[2] ./ mu
        du[4] = - SR211 .* u[1] ./ mu .- SR322 .* u[2] ./ mu
        du[4] += k322BH .* alph.^11 .* (M_pl ./ fa).^4 .* rP .* u[1].^2 .* u[2]
        
        # print(t, "\t", u[1] ./ Emax2, "\t", u[2] .* mu .* (GNew .* u[4]), "\t", u[3], "\t", u[4], "\n")
        
        du[1] *= mu ./ hbar .* 3.15e7
        du[2] *= mu ./ hbar .* 3.15e7
        du[3] *= mu ./ hbar .* 3.15e7
        du[4] *= mu.^2 .* (GNew .* u[4].^2) ./ hbar .* 3.15e7
        
        
        # print("Test \t", u[1] ./ du[1], "\t", u[2] ./ du[2], "\n")
        return
    end

    function c_e2max(u, t, integrator)
        return u[1] .- Emax2
    end
    function affect_e2max!(integrator)
        integrator.u[1] = Emax2
        if !u1_eq
            set_proposed_dt!(integrator, integrator.dt .* 0.2)
        end
    end
    function check_timescale(u, t, integrator)
        
        alph = GNew .* u[4] .* mu
        rP = nothing
        if u[3] .> 0.998
            rP = 1.0 .+ sqrt.(1 - 0.998 .^2)
        elseif u[3] .< 0.0
            rP = 2.0
            u[3] = 0.0
        else
            rP = 1.0 .+ sqrt.(1 - u[3].^2)
        end
        du = get_du(integrator)
        
      
        t1 = abs.(u[1] ./ du[1])
#        if (OmegaH .< ergL(2, 1, 1, mu, u[4]))||(u[1] < (GNew .* u[4].^2 .* M_to_eV).^(-1))
#            t1 = 1e100
#        end
        t2 = abs.(u[2] ./ du[2])
#        if (2 .* OmegaH .< ergL(3, 2, 2, mu, u[4]))
#            t2 = 1e100
#        end

        if u1_eq || u2_eq
            SR211 = sr_rates(2, 1, 1, mu, u[4], u[3], impose_low_cut=impose_low_cut, solve_322=solve_322)
            SR322 = sr_rates(3, 2, 2, mu, u[4], u[3], impose_low_cut=impose_low_cut, solve_322=solve_322)
            
            
            if u1_eq
                t1 = 1e100
                integrator.u[1] = u1_fix
            end
            if u2_eq
                t2 = 1e100
                integrator.u[2] = u2_fix
            end
            
            du[3] = - SR211 .* u[1] ./ mu .- 2 .* SR322 .* u[2] ./ mu
        end
        t3 = abs.(u[3] ./ du[3])
        t4 = abs.(u[4] ./ du[4])
        
       
#        tcheck = minimum([t1 t2 t3 t4])
        tcheck = minimum([t1 t2 t3])
        wait += 1
        
        
        if debug && (wait%10==0)
            # print("CHECK \t", integrator.dt, "\t", u[1] ./ du[1], "\t", u[2] ./ du[2], "\n")
            # print("CHECK \t", integrator.dt, "\t", t1, "\t", t2, "\t", t3, "\t", t4, "\n")
            # k322BH = 4e-7  # k^322xBH_211x211
            # GR_322 = sr_rates(3, 2, 2, mu, u[4], u[3], impose_low_cut=impose_low_cut, solve_322=solve_322) .* u[2] .+ k322BH .* alph.^11 .* (M_pl ./ fa).^4 .* rP .* u[1].^2 .* u[2]
            # print(t, "\t", u[1], "\t", u[2], "\t", u[3], "\t", u[4], "\t", u1_eq, "\t", u2_eq, "\t", du[2] ./ (GR_322 ./ hbar .* 3.15e7), "\n\n")
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
        
        du = get_du(integrator)

        t1 = abs.(integrator.u[1] ./ du[1])
        # OmegaH = integrator.u[3] ./ (2 .* (GNew .* integrator.u[4]) .* (1 .+ sqrt.(1 .- integrator.u[3].^2)))
#        if (OmegaH .< ergL(2, 1, 1, mu, integrator.u[4]))||(integrator.u[1] < (GNew .* integrator.u[4].^2 .* M_to_eV).^(-1))
#            t1 = 1e100
#        end
        t2 = abs.(integrator.u[2] ./ du[2])
#        if (2 .* OmegaH .< ergL(3, 2, 2, mu, integrator.u[4]))
#            t2 = 1e100
#        end
        if u1_eq || u2_eq
            if u1_eq
                t1 = 1e100
                integrator.u[1] = u1_fix
            end
            if u2_eq
                t2 = 1e100
                integrator.u[2] = u2_fix
            end
            SR211 = sr_rates(2, 1, 1, mu, integrator.u[4], integrator.u[3], impose_low_cut=impose_low_cut, solve_322=solve_322)
            SR322 = sr_rates(3, 2, 2, mu, integrator.u[4], integrator.u[3], impose_low_cut=impose_low_cut, solve_322=solve_322)
            
            du[3] = - SR211 .* integrator.u[1] ./ mu .- 2 .* SR322 .* integrator.u[2] ./ mu
        end
        
        t3 = abs.(integrator.u[3] ./ du[3])
#        if u1_eq&&u2_eq&&(t3 > 100.0 * t_max)
#            terminate!(integrator)
#        end
        # t4 = abs.(integrator.u[4] ./ du[4])
    
        
#        tcheck = minimum([t1 t2 t3 t4])
        tcheck = minimum([t1 t2 t3])
        
        
        if (tcheck .<= integrator.dt)
            set_proposed_dt!(integrator, tcheck .* 0.1)
        elseif (integrator.dt .<= 1e-7)
            terminate!(integrator)
#        else
        elseif (wait % 5 == 0)
            set_proposed_dt!(integrator, integrator.dt .* 1.05)
        end
    end
    
    function check_eq(u, t, integrator)
        if u2_kill
            terminate!(integrator)
        end
        if u1_eq&&u2_eq
            return true
        end
        du = get_du(integrator)
        t3 = abs.(integrator.u[3] ./ du[3])
        # watch out for stable equilibrium of 322 state
        # cVal1 = log.(u[1])
        # cVal2 = log.(u[2])
        alph = GNew .* u[4] .* mu
        SR211 = sr_rates(2, 1, 1, mu, integrator.u[4], integrator.u[3], impose_low_cut=impose_low_cut, solve_322=solve_322)
        SR322 = sr_rates(3, 2, 2, mu, integrator.u[4], integrator.u[3], impose_low_cut=impose_low_cut, solve_322=solve_322)
        k322BH = 4e-7  # k^322xBH_211x211
        GR_322 = SR322 .* integrator.u[2] .+ k322BH .* alph.^11 .* (M_pl ./ fa).^4 .* rP .* integrator.u[1].^2 .* integrator.u[2]
        if (abs.(du[1] ./ (SR211 .* integrator.u[1] ./ hbar .* 3.15e7)) .< 1e-4)&&(abs.(du[2] ./ (GR_322 ./ hbar .* 3.15e7)) .< 1e-4)
            u1_eq = true
            u2_eq = true
            u1_fix = integrator.u[1]
            u2_fix = integrator.u[2]
            u4_fix = integrator.u[4]
        end
        
        if u1_eq||u2_eq
            return true
        else
            return false
        end
    end
    function affect_eq!(integrator)
    
        du = get_du(integrator)
        if u1_eq && !u2_eq
            du[1] = 0.0
            du[4] = 0.0
            set_u!(integrator, [u1_fix, integrator.u[2], integrator.u[3], u4_fix])
        elseif u2_eq && !u1_eq
            du[2] = 0.0
            du[4] = 0.0
            set_u!(integrator, [integrator.u[1], u2_fix, integrator.u[3], u4_fix])
        elseif u1_eq && u2_eq
            du[1] = 0.0
            du[2] = 0.0
            du[4] = 0.0
            set_u!(integrator, [u1_fix, u2_fix, integrator.u[3], u4_fix])
            
        end
    end
    

    function check_spin(u, t, integrator)
        if u[3] .> (aBH .+ 0.01)
            return true
        elseif u[3] .<= 0.0
            return true
        else
            return false
        end
    end
    function affect_spin!(integrator)
        if integrator.u[3] .> aBH
            integrator.u[3] = aBH
        elseif integrator.u[3] .< 0.0
            integrator.u[3] = 0.0
        end
        set_proposed_dt!(integrator, integrator.dt .* 0.3)
    end
    
    
        
    cbackEmax = ContinuousCallback(c_e2max, affect_e2max!, interp_points=10, abstol=0.01)
    cbackdt = DiscreteCallback(check_timescale, affect_timescale!, save_positions=(false, true))
    cback_equil = DiscreteCallback(check_eq, affect_eq!, save_positions=(false, true))
    cbackspin = DiscreteCallback(check_spin, affect_spin!, save_positions=(false, true))
    
    # cbset = CallbackSet(cback_equil, cbackEmax, cbackdt, cbackspin)
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
    spinBH = [sol.u[i][3] for i in 1:length(sol.u)]
    MassB = [sol.u[i][4] for i in 1:length(sol.u)]
    if debug
        print("Initial \t", state211[1] ./ Emax2, "\t", state322[1] ./ Emax2, "\t", spinBH[1], "\t", MassB[1], "\n")
        print("Final \t", state211[end] ./ Emax2, "\t", state322[end] ./ Emax2, "\t", spinBH[end], "\t", MassB[end], "\n")
    end

    if return_all_info
        return sol.t, state211, state322, spinBH, MassB
    end

    if isnan(spinBH[end])
        spinBH = spinBH[.!isnan.(spinBH)]
    end
    return spinBH[end]
 
end

function ergL(n, l, m, massB, MBH)
    alph = GNew .* MBH .* massB
    # return massB .* (1.0 .- alph.^2 ./ (2 .* (n .+ l .+ 1).^2))
    return massB .* (1.0 .- alph.^2 ./ (2 .* n.^2))
end


function sr_rates(n, l, m, massB, MBH, aBH; impose_low_cut=0.01, solve_322=true)
    if (n==3)&&(l==2)&&(m==1)&&(solve_322==false)
        return 0.0
    end
    
    alph = GNew .* MBH .* massB
    
    if (alph ./ l < impose_low_cut)&&(MBH < 1e2)
        # We expect binaries to disrupt.
        return 0.0
    end
    
    
    rP = nothing
    if aBH .> 0.998
        rP = 1.0 .+ sqrt.(1 - 0.998 .^2)
        aBH = 0.998
    elseif aBH .< 0.0
        rP = 2.0
        aBH = 0.0
    else
        rP = 1.0 .+ sqrt.(1 - aBH.^2)
    end
    rP *= (GNew .* MBH)
    OmegaH = aBH ./ (2 .* (GNew .* MBH) .* (1 .+ sqrt.(1 .- aBH.^2)))
    Anl = 2 .^(4 .* l .+ 1) .* factorial(Int(l .+ n)) ./ (n.^(2 .* l .+ 4) .* factorial(n .- l .- 1))
    Anl *= (factorial(Int(l)) ./ (factorial(Int(2 .* l)) .* factorial(Int(2 .* l .+ 1)))).^2
    Chilm = 1.0
    # erg = ergL(n, l, m, massB, MBH)
    for k in 1:Int(l)
        Chilm *= (k.^2 .* (1.0 .- aBH.^2) .+ (aBH * m .- 2 .* (rP ./ (GNew .* MBH)) .* alph).^2)
        # Chilm *= (k.^2 .* (1.0 .- aBH.^2) .+ (aBH * m .- 2 .* (rP ./ (GNew .* MBH)) .* alph).^2)
    end
    Gamma_nlm = 2 * massB .* rP .* (m .* OmegaH .- ergL(n, l, m, massB, MBH)) .* alph.^(4 .* l + 4) .* Anl .* Chilm
    
#    if (n==2)&&(l==1)&&(m==1)
#        Gamma_nlm = 4e-2 .* alph.^8 .* (aBH .- 2 .* alph .* (1 .+ sqrt.(1.0 .- aBH.^2))) .* massB
#    end
    if Gamma_nlm > 0.0
        return Gamma_nlm
    else
        return 0.0
    end
end





#### TESTING ZONE
M_BH = 9.0850432104
aBH = 0.907515615
massB = 2.813925014e-12
f_a = 6.15407859e13
tau_max = 1e8
alpha_max_cut = 1.0
solve_322 = true
impose_low_cut=0.01
debug=true
# fs = super_rad_check(M_BH, aBH, massB, f_a, tau_max=tau_max, alpha_max_cut=alpha_max_cut, debug=debug, solve_322=solve_322, impose_low_cut=impose_low_cut)
# print("Final spin \t", fs, "\n")
########################

