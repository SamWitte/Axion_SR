using Random
using OrdinaryDiffEq
using Statistics
using Distributions
using DelimitedFiles
include("Constants.jl")




function super_rad_check(M_BH, aBH, massB, f_a; spin=0, tau_max=1e4, alpha_max_cut=0.2, debug=false)
   
    alph = GNew .* M_BH .* massB #
    if debug
        print("Alpha \t", alph, "\n")
    end
    if alph .> alpha_max_cut
        if debug
            print("Need higher-level system... \n")
        end
        return aBH
    end
    
    final_spin = solve_system(massB, f_a, aBH, M_BH, tau_max, debug=debug)
    # print("Spin diff.. \t ", aBH, "\t", final_spin, "\t", alph, "\n")
    return final_spin
    
end

function emax_211(MBH, mu, aBH)
    alph = GNew .* MBH .* mu

    emax_N = 1 .- 8 .* alph.^2 .+ 8 .* alph.^3 .* aBH .- sqrt.(1 .- 16 .* alph.^2 .+ 32 .* aBH .* alph.^3 - 16 .* aBH.^2 .* alph.^4)
    emax_D = 8 .* (-alph.^3 .+ aBH .* alph.^4)
    return (emax_N ./ emax_D)
end
    
function solve_system(mu, fa, aBH, M_BH, t_max; n_times=100, debug=true)
    
    
    y0 = [1.0 ./ (GNew .* M_BH.^2 .* M_to_eV), 1.0 ./ (GNew .* M_BH.^2 .* M_to_eV), aBH, M_BH]
    wait = 0
    Emax2 = emax_211(M_BH, mu, aBH)
    
    Mvars = [mu, fa, Emax2]
    tspan = (0.0, t_max)
    saveat = (tspan[2] .- tspan[1]) ./ n_times
    # print("T max \t", t_max, "\n")
    function c_e2max(u, t, integrator)
        return u[1] .- Emax2
    end
    function affect_e2max!(integrator)
        integrator.u[1] = Emax2
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
        t2 = abs.(u[2] ./ du[2])
        t3 = abs.(u[3] ./ du[3])
        t4 = abs.(u[4] ./ du[4])
        tcheck = minimum([t1 t2 t3 t4])
        wait += 1
        
    
        # print("CHECK \t", integrator.dt, "\t", u[1] ./ du[1], "\t", u[2] ./ du[2], "\t", tcheck, "\n")
        
        if (tcheck .>= 50.0 .* integrator.dt) && (wait % 100 == 0)
            # print("here \t", integrator.dt, "\t", wait, "\n")
            return true
        elseif (tcheck .<= integrator.dt)
            return true
        else
            return false
        end
    end
    function affect_timescale!(integrator)
        du = get_du(integrator)
        t1 = abs.(integrator.u[1] ./ du[1])
        t2 = abs.(integrator.u[2] ./ du[2])
        t3 = abs.(integrator.u[3] ./ du[3])
        t4 = abs.(integrator.u[4] ./ du[4])
        tcheck = minimum([t1 t2 t3 t4])
        if (tcheck .<= integrator.dt)
            set_proposed_dt!(integrator, integrator.dt .* 0.5)
        else
            set_proposed_dt!(integrator, integrator.dt .* 1.1)
        end
    end
    
    function check_spin(u, t, integrator)
        if u[3] .> aBH
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
        
    cbackEmax = ContinuousCallback(c_e2max, affect_e2max!, interp_points=10, abstol=1e-2)
    cbackdt = DiscreteCallback(check_timescale, affect_timescale!)
    cbackspin = DiscreteCallback(check_spin, affect_spin!)
    cbset = CallbackSet(cbackEmax, cbackdt, cbackspin)
    prob = ODEProblem(RHS_ax!, y0, tspan, Mvars, reltol=1e-6, abstol=1e-6)
    # sol = solve(prob, Vern6(), saveat=saveat, callback=cbset)
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
    # sol = solve(prob, Vern6(), dt=dt_guess, saveat=saveat, callback=cbset)
    
    state211 = [sol.u[i][1] for i in 1:length(sol.u)]
    state322 = [sol.u[i][2] for i in 1:length(sol.u)]
    spinBH = [sol.u[i][3] for i in 1:length(sol.u)]
    MassB = [sol.u[i][4] for i in 1:length(sol.u)]
    if debug
        print("Initial \t", state211[1], "\t", state322[1], "\t", spinBH[1], "\t", MassB[1], "\n")
        print("Final \t", state211[end], "\t", state322[end], "\t", spinBH[end], "\t", MassB[end], "\n")
        print("Giving back \t",spinBH[end], "\n" )
    end
    if isnan(spinBH[end])
        spinBH = spinBH[.!isnan.(spinBH)]
    end
    return spinBH[end]
end




function RHS_ax!(du, u, Mvars, t)

    # [e211, e322, aBH, MBH]
    # mu [mass eV], fa[1/GeV]
    
    mu, fa, Emax2  = Mvars
    alph = GNew .* u[4] .* mu #
    rP = nothing
    if u[3] .> 0.998
        rP = 1.0 .+ sqrt.(1 - 0.998 .^2)
    elseif u[3] .< 0.0
        rP = 2.0
        u[3] = 0.0
    else
        rP = 1.0 .+ sqrt.(1 - u[3].^2)
    end
    
    
    
    kSR_211 = 4e-2 # kSR_211
    if u[1] .> Emax2
        kSR_211 *= 0.0
    end
    k322BH = 4e-7  # k^322xBH_211x211
    k2I_333 = 1e-8 # k^211xInfinity_322x322
    kGW_22 = 1e-2 # k^GW_211x211
    kGW_3t2 = 5e-6 # k^GW_{322->211}
    kGW_23 = 0.0 # k^GW_211x322
    kI_222 = 1.5e-8 # k^Infinity_(211)^3
    kI_223 = 1.5e-8 # ??? k^Infinity_{(211)^2 x 322}
    kI_233 = 1.5e-8 # ??? k^Infinity_{(211)x 322^2}
    kSR_322 = 8e-5
    kGW_33 = 3e-8
    kI_333 = 1.5e-8  # ???
    
    
    du[1] = kSR_211 .* alph.^8 .* (u[3] .- 2 .* alph .* rP) .* u[1]
    du[1] += - 2 .* k322BH .* alph.^11 .* (M_pl ./ fa).^4 .* rP .* u[1].^2 .* u[2]
    du[1] += k2I_333 .* alph.^8 .* (M_pl ./ fa).^4 .* u[2].^2 .* u[1]
    du[1] += -2 .* kGW_22 .* alph.^14 .* u[1].^2
    du[1] += -kGW_23 .* alph.^16 .* u[1] .* u[2] .+ kGW_3t2 .* alph.^10 .* u[1] .* u[2]
    du[1] += -3 .* kI_222 .* alph.^21 .* (M_pl ./ fa).^4 .* u[1].^3
    du[1] += -3 .* kI_222 .* alph.^21 .* (M_pl ./ fa).^4 .* u[1].^3
    du[1] += -2 .* kI_223 .* alph.^23 .* (M_pl ./ fa).^4 .* u[1].^2 .* u[2]
    du[1] += kI_233 .* alph.^25 .* (M_pl ./ fa).^4 .* u[1] .* u[2].^2
    
    
    du[2] = kSR_322 .* alph.^12 .* (u[3] .- alph .* rP) .* u[2]
    du[2] += k322BH .* alph.^11 .* (M_pl ./ fa).^4 .* rP .* u[1].^2 .* u[2]
    du[2] += -2 .* k2I_333 .* alph.^8 .* (M_pl ./ fa).^4 .* u[2].^2 .* u[1]
    du[2] += -2 .* kGW_33 .* alph.^18 .* u[2].^2
    du[2] += -kGW_23 .* alph.^16 .* u[1] .* u[2] .- kGW_3t2 .* alph.^10 .* u[1] .* u[2]
    du[2] += -3 .* kI_333 .* alph.^27 .* (M_pl ./ fa).^4 .* u[2].^3
    du[2] += -kI_223 .* alph.^23 .* (M_pl ./ fa).^4 .* u[1].^2 .* u[2]
    du[2] += -2 .* kI_233 .* alph.^25 .* (M_pl ./ fa).^4 .* u[1] .* u[2].^2
    
    du[3] = - kSR_211 .* alph.^8 .* (u[3] .- 2 .* alph .* rP) .* u[1] .- 2 .* kSR_322 .* alph.^12 .* (u[3] .- alph .* rP) .* u[2]
    du[4] = - kSR_211 .* alph.^8 .* (u[3] .- 2 .* alph .* rP) .* u[1] .- kSR_322 .* alph.^12 .* (u[3] .- alph .* rP) .* u[2]
    du[4] += k322BH .* alph.^11 .* (M_pl ./ fa).^4 .* rP .* u[1].^2 .* u[2]
    
    # print(t, "\t", u[1] ./ Emax2, "\t", u[2] .* mu .* (GNew .* u[4]), "\t", u[3], "\t", u[4], "\n")
    
    du[1] *= mu ./ hbar .* 3.15e7
    du[2] *= mu ./ hbar .* 3.15e7
    du[3] *= mu ./ hbar .* 3.15e7
    du[4] *= mu.^2 .* (GNew .* u[4].^2) ./ hbar .* 3.15e7
    
    # print("Test \t", u[1] ./ du[1], "\t", u[2] ./ du[2], "\n")
    return
end



#### TESTING ZONE
# M_BH = 1e8
# aBH = 0.9
# massB = 5e-19
# f_a = 1e19
# tau_max = 1e8
# alpha_max_cut = 0.5
# super_rad_check(M_BH, aBH, massB, f_a, tau_max=tau_max, alpha_max_cut=alpha_max_cut, debug=true)
########################

