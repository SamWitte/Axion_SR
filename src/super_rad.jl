using Random
using OrdinaryDiffEq
using Statistics
using Distributions
using DelimitedFiles
using Interpolations
include("Constants.jl")
include("solve_sr_rates.jl")


function super_rad_check(M_BH, aBH, massB, f_a; spin=0, tau_max=1e4, alpha_max_cut=0.2, debug=false, solve_322=true, solve_n4=false, impose_low_cut=0.01, input_data="Masha", stop_on_a=0, eq_threshold=1e-4)
   
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
    
    final_spin = solve_system(massB, f_a, aBH, M_BH, tau_max, debug=debug, solve_322=solve_322, impose_low_cut=impose_low_cut, input_data=input_data, solve_n4=solve_n4, stop_on_a=stop_on_a, eq_threshold=eq_threshold)
    # print("Spin diff.. \t ", aBH, "\t", final_spin, "\t", alph, "\n")
    return final_spin
    
end

function emax_211(MBH, mu, aBH)
    alph = GNew .* MBH .* mu
    
    emax_N = 1 .- 8 .* alph.^2 .+ 8 .* alph.^3 .* aBH .- sqrt.(abs.(1 .- 16 .* alph.^2 .+ 32 .* aBH .* alph.^3 - 16 .* aBH.^2 .* alph.^4)) # abs just in case...
    emax_D = 8 .* (-alph.^3 .+ aBH .* alph.^4)
    return (emax_N ./ emax_D)
end
    
function solve_system(mu, fa, aBH, M_BH, t_max; n_times=100, debug=true, solve_322=true, impose_low_cut=0.01, return_all_info=false, input_data="Masha", solve_n4=false, eq_threshold=1e-4, stop_on_a=0)
    e_init = 1.0 ./ (GNew .* M_BH.^2 .* M_to_eV) # unitless
    if !solve_n4
        y0 = [e_init, e_init, aBH, M_BH]
    else
        y0 = [e_init, e_init, e_init, e_init, e_init, aBH, M_BH]
    end
    wait = 0


    u1_eq = false
    u1_fix = nothing
    u2_eq = false
    u2_fix = nothing
    u3_fix = nothing
    u4_fix = nothing
    u5_fix = nothing
    u6_fix = nothing
    u2_kill = false
    
    if input_data != "Doddy"
        e2_maxBN = 1024 * pi * (fa / M_pl).^2 ./ (9 * (GNew .* M_BH .* mu ).^3)
    else
        e2_maxBN = (5 .* 1e78 .* (2.0 .^4 ./ alph.^3) .* (MassBH ./ 10.0).^2 .* (10 .^log_f ./ M_pl).^2) .* e_init
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
    
    alph = GNew .* M_BH .* mu
    
    N_pts_interp = 50
    xtol_slv = 1e-30
    Ntot_slv = 2000
    iter_slv = 50
    
    alist, pts211 = compute_gridded(mu, M_BH, aBH, 2, 1, 1; Ntot=Ntot_slv, iter=iter_slv, xtol=xtol_slv, npts=N_pts_interp)
    itp_211 = LinearInterpolation(alist, log10.(pts211), extrapolation_bc=Line())
    SR211 = 10 .^itp_211(aBH)
    
    alist, pts322 = compute_gridded(mu, M_BH, aBH, 3, 2, 2; Ntot=Ntot_slv, iter=iter_slv, xtol=xtol_slv, npts=N_pts_interp)
    itp_322 = LinearInterpolation(alist, log10.(pts322), extrapolation_bc=Line())
    SR322 = 10 .^itp_322(aBH)
    
    if solve_n4
        alist, pts411 = compute_gridded(mu, M_BH, aBH, 4, 1, 1; Ntot=Ntot_slv, iter=iter_slv, xtol=xtol_slv, npts=N_pts_interp)
        itp_411 = LinearInterpolation(alist, log10.(pts411), extrapolation_bc=Line())
        SR411 = 10 .^itp_411(aBH)
        
        alist, pts422 = compute_gridded(mu, M_BH, aBH, 4, 2, 2; Ntot=Ntot_slv, iter=iter_slv, xtol=xtol_slv, npts=N_pts_interp)
        itp_422 = LinearInterpolation(alist, log10.(pts422), extrapolation_bc=Line())
        SR422 = 10 .^itp_422(aBH)
        
        alist, pts433 = compute_gridded(mu, M_BH, aBH, 4, 3, 3; Ntot=Ntot_slv, iter=iter_slv, xtol=xtol_slv, npts=N_pts_interp)
        itp_433 = LinearInterpolation(alist, log10.(pts433), extrapolation_bc=Line())
        SR433 = 10 .^itp_433(aBH)
    end

    function RHS_ax!(du, u, Mvars, t)
    
        # [e211, e322, aBH, MBH] or [e211, e322, e411, aBH, MBH]
        # mu [mass eV], fa[1/GeV]
        mu, fa, Emax2, solve_322, aBH_i, M_BH_i, impose_low_cut  = Mvars
        
        if !solve_n4
            spinI = 3
            massI = 4
        else
            spinI = 6
            massI = 7
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
            if u[4] < 0
                u[4] = 0.0
            end
            if u[5] < 0
                u[5] = 0.0
            end

        end

        OmegaH = u[spinI] ./ (2 .* (GNew .* u[massI]) .* (1 .+ sqrt.(1 .- u[spinI].^2)))
        

        if (OmegaH .< ergL(2, 1, 1, mu, u[massI], u[spinI]))
            SR211 = 0.0
        else
            SR211 = 10 .^itp_211(u[spinI])
        end
        
        if (2 .* OmegaH .< ergL(3, 2, 2, mu, u[massI], u[spinI]))
            SR322 = 0.0
        else
            SR322 = 10 .^itp_322(u[spinI])
        end
        
        if solve_n4
            if (u[3] > 1e-100)&&(OmegaH .> ergL(4, 1, 1, mu, u[massI], u[spinI]))
                SR411 = 10 .^itp_411(u[spinI])
            else
                SR411 = 0.0
            end
        
            if (u[4] > 1e-100)&&(OmegaH * 2 .> ergL(4, 2, 2, mu, u[massI], u[spinI]))
                SR422 = 10 .^itp_422(u[spinI])
            else
                SR422 = 0.0
            end
        
            if (u[5] > 1e-100)&&(OmegaH * 3 .> ergL(4, 3, 3, mu, u[massI], u[spinI]))
                SR433 = 10 .^itp_433(u[spinI])
            else
                SR433 = 0.0
            end
        end

        
        if u[1] .> Emax2
            SR211 *= 0.0
        end
        
        if input_data == "Doddy"
            k322BH = 0.0
            k2I_333 = 0.0
            kGW_22 = 0.0
            kGW_3t2 = 0.0
            kI_222 = 0.0
            kGW_33 = 0.0
        else
            k322BH = 4e-7  # k^322xBH_211x211
            k2I_333 = 1e-8 # k^211xInfinity_322x322
            kGW_22 = 1e-2 # k^GW_211x211
            kGW_3t2 = 5e-6 # k^GW_{322->211}
            kI_222 = 1.5e-8 # k^Infinity_(211)^3
            kGW_33 = 3e-8
        end
        
        kGW_23 = 0.0 # k^GW_211x322
        kI_223 = 0.0 # ??? k^Infinity_{(211)^2 x 322}
        kI_233 = 0.0 # ??? k^Infinity_{(211)x 322^2}
        kI_333 = 0.0 # ???
        

        if isnothing(u1_fix)
            du[1] = SR211 .* u[1] ./ mu
            du[1] += - 2 .* k322BH .* alph.^11 .* (M_pl ./ fa).^4 .* rP .* u[1].^2 .* u[2]
            du[1] += k2I_333 .* alph.^8 .* (M_pl ./ fa).^4 .* u[2].^2 .* u[1]
            du[1] += -2 .* kGW_22 .* alph.^14 .* u[1].^2
            du[1] += -3 .* kI_222 .* alph.^21 .* (M_pl ./ fa).^4 .* u[1].^3
            du[1] += kGW_3t2 .* alph.^10 .* u[1] .* u[2]
        else
            du[1] = 0.0
        end
        

        
        if isnothing(u2_fix)&&(solve_322==true)
            du[2] = SR322 .* u[2] ./ mu
            du[2] += k322BH .* alph.^11 .* (M_pl ./ fa).^4 .* rP .* u[1].^2 .* u[2]
            du[2] += -2 .* k2I_333 .* alph.^8 .* (M_pl ./ fa).^4 .* u[2].^2 .* u[1]
            du[2] += -2 .* kGW_33 .* alph.^18 .* u[2].^2
            du[2] += -kGW_23 .* alph.^16 .* u[1] .* u[2] .- kGW_3t2 .* alph.^10 .* u[1] .* u[2]
        else
            du[2] = 0.0
        end
        
        if solve_n4
            # 411
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
            
            ### 422
            
            du[4] = SR422 .* u[4] ./ mu
            
            du[4] += 1.5e-7 .* alph.^11 .* (M_pl ./ fa).^4 .* rP .* u[1].^2 .* u[4] # 211 x 211 -> 422 x BH
            du[1] += -2 * 1.5e-7 .* alph.^11 .* (M_pl ./ fa).^4 .* rP .* u[1].^2 .* u[4] # 211 x 211 -> 422 x BH
            
            du[4] += -5e-12 .* alph.^8 .* (M_pl ./ fa).^4 .* u[4] .* u[3] .* u[1] # 411 x 422 -> 211 x Inf
            du[3] += -5e-12 .* alph.^8 .* (M_pl ./ fa).^4 .* u[4] .* u[3] .* u[1] # 411 x 422 -> 211 x Inf
            du[1] +=  5e-12 .* alph.^8 .* (M_pl ./ fa).^4 .* u[4] .* u[3] .* u[1] # 411 x 422 -> 211 x Inf
            
            du[4] += 2.3e-7 .* alph.^7 .* (M_pl ./ fa).^4 .* rP .* u[3].^2 .* u[4] # 411 x 411 -> 422 x BH
            du[3] += - 2 * 2.3e-7 .* alph.^7 .* (M_pl ./ fa).^4 .* rP .* u[3].^2 .* u[4] # 411 x 411 -> 422 x BH
            
            
            du[4] += -1.1e-9 .* alph.^7 .* (M_pl ./ fa).^4 .* rP .* u[1] .* u[4] .* u[5]  # 211 x 422 -> 433 x 200 / BH
            du[1] += -1.1e-9 .* alph.^7 .* (M_pl ./ fa).^4 .* rP .* u[1] .* u[4] .* u[5]  # 211 x 422 -> 433 x 200 / BH
            
   
            ### 433
            du[5] = SR433 .* u[5] ./ mu
            
            du[5] += 1.1e-9 .* alph.^7 .* (M_pl ./ fa).^4 .* rP .* u[1] .* u[4] .* u[5] # 211 x 422 -> 433 x 200 [??]
             
            du[5] += -2 * 9.2e-11 .* alph.^8 .* (M_pl ./ fa).^4  .* u[1] .* u[5].^2 # 433 x 433 -> 211 x inf
            du[1] += 9.2e-11 .* alph.^8 .* (M_pl ./ fa).^4  .* u[1] .* u[5].^2
            
            du[5] += - 2.6e-9 .* alph.^8 .* (M_pl ./ fa).^4  .* u[1] .* u[5] .* u[2] # 322 x 433 -> 211 x inf
            du[1] +=   2.6e-9 .* alph.^8 .* (M_pl ./ fa).^4  .* u[1] .* u[5] .* u[2]
            du[2] += - 2.6e-9 .* alph.^8 .* (M_pl ./ fa).^4  .* u[1] .* u[5] .* u[2]
            
            du[5] += 9.1e-8 .* alph.^11 .* (M_pl ./ fa).^4  .* u[1]  .* u[2] .* u[5] # 211 x 322 -> 433 x BH
            du[1] += -9.1e-8 .* alph.^11 .* (M_pl ./ fa).^4  .* u[1]  .* u[2] .* u[5]
            du[2] += -9.1e-8 .* alph.^11 .* (M_pl ./ fa).^4  .* u[1]  .* u[2] .* u[5]
            
            
            # test =  -9.1e-8 .* alph.^11 .* (M_pl ./ fa).^4  .* u[1]  .* u[2] .* u[5]
            # test2 =  -1.1e-9 .* alph.^7 .* (M_pl ./ fa).^4 .* rP .* u[1] .* u[4] .* u[5]  # 211 x 422 -> 433 x 200 / BH
            # test3 = -2.5e-8 .* alph.^11 .* (M_pl ./ fa).^4 .* rP .* u[1] .* u[2] .* u[3]
            # print(test, "\t", test2, "\t", test3, "\t", du[1], "\n")
            
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
            du[spinI] = - SR211 .* u[1] ./ mu .- 2 .* SR322 .* u[2] ./ mu  .- SR411 .* u[3] ./ mu .- SR422 .* u[4] ./ mu .- SR433 .* u[5] ./ mu
            # print(SR211, "\t", SR322, "\t", SR411, "\t", SR422, "\t", SR433, "\n" )
            du[massI] = - SR211 .* u[1] ./ mu .-  SR322 .* u[2] ./ mu  .- SR411 .* u[3] ./ mu .- SR422 .* u[4] ./ mu .- SR433 .* u[5] ./ mu
            du[massI] += k322BH .* alph.^11 .* (M_pl ./ fa).^4 .* rP .* u[1].^2 .* u[2]
            du[massI] += 9.8e-11 .* alph.^11 .* (M_pl ./ fa).^4 .* rP .* u[2] .* u[3].^2
            du[massI] += 1.5e-7 .* alph.^11 .* (M_pl ./ fa).^4 .* rP .* u[1].^2 .* u[4]
            du[massI] += 9.1e-8 .* alph.^11 .* (M_pl ./ fa).^4  .* u[1]  .* u[2] .* u[5]
            du[massI] += 2.3e-7 .* alph.^7 .* (M_pl ./ fa).^4 .* rP .* u[3].^2 .* u[4] # 411 x 411 -> 422 x BH

            
            du[1] *= mu ./ hbar .* 3.15e7
            du[2] *= mu ./ hbar .* 3.15e7
            du[3] *= mu ./ hbar .* 3.15e7
            du[4] *= mu ./ hbar .* 3.15e7
            du[5] *= mu ./ hbar .* 3.15e7
            
            du[spinI] *= mu ./ hbar .* 3.15e7
            du[massI] *= (mu .* u[massI]) .* (mu .* GNew .* u[massI]) ./ hbar .* 3.15e7
            
            ## check bosenova
            if (u[1] > e2_maxBN)
                if du[1] > 0
                    du[1] = 0.0
                end
            end
            if (u[2] > e3_maxBN)
                if du[2] > 0
                    du[2] = 0.0
                end
            end
            if (u[3] > e4_maxBN)
                if du[3] > 0
                    du[3] = 0.0
                end
            end
            if (u[4] > e4_maxBN)
                if du[4] > 0
                    du[4] = 0.0
                end
            end
            if (u[5] > e4_maxBN)
                if du[5] > 0
                    du[5] = 0.0
                end
            end            
            
        end
        
        return
    end


    function check_timescale(u, t, integrator)
        
        
        if !solve_n4
            spinI = 3
            massI = 4
        else
            spinI = 6
            massI = 7
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
        
      
        
        
        if (u[1] > 0)
            t1 = abs.(u[1] ./ du[1])
        else
            t1 = 1e100
        end
        if (u[2] > 0)
            t2 = abs.(u[2] ./ du[2])
        else
            t2 = 1e100
        end
        if solve_n4
            if (integrator.u[3] > 0.0)
                t5 = abs.(u[3] ./ du[3])
            else
                t5 = 1e100
            end
            if (integrator.u[4] > 0.0)
                t6 = abs.(u[4] ./ du[4])
            else
                t6 = 1e100
            end
            if (integrator.u[5] > 0.0)
                t7 = abs.(u[5] ./ du[5])
            else
                t7 = 1e100
            end
        end
        
        if (u1_eq)||(u2_eq)      
            t1 = 1e100
            integrator.u[1] = u1_fix
            t2 = 1e100
            integrator.u[2] = u2_fix
            if u2_eq && solve_n4
                t5 = 1e100
                t6 = 1e100
                t7 = 1e100
                integrator.u[3] = u3_fix
                integrator.u[4] = u5_fix
                integrator.u[5] = u6_fix

            end
            
            du[spinI] = - SR211 .* u[1] ./ mu .-  2 .* SR322 .* u[2] ./ mu
            if solve_n4
                du[spinI] += - SR411 .* u[3] ./ mu - SR422 .* u[4] ./ mu - SR433 .* u[5] ./ mu
            end
        end
        t3 = abs.(u[spinI] ./ du[spinI])
        t4 = abs.(u[massI] ./ du[massI])
        
       
        
        if !solve_n4
            tcheck = minimum([t1 t2 t3])
        else
            # tcheck = minimum([t1 t2 t3 t5 t6 t7])
            tcheck = minimum([t1 t2 t5 t6 t7])
            if debug && (wait%10==0)
            #    print(integrator.dt, "\t", [t1 t2 t5 t6 t7], "\n\n")
            end
        end
        
        
        wait += 1
        
        
        if debug && (wait%10==0)
            if solve_n4
                print(t, "\t", u[1], "\t", u[2], "\t", u[3], "\t", u[4], "\t", u[5], "\t", u[6], "\n")
            else
                print(t, "\t", u[1], "\t", u[2], "\t", u[3], "\t", u[4], "\n")
            end
        end
        
        if (integrator.dt ./ tcheck .<= 1e-1) && (wait % 100 == 0) # can increase time step
            return true
        elseif (integrator.dt ./ tcheck .>= 1e-1)
            return true
        elseif (integrator.dt .<= 1e-10)
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
            spinI = 6
            massI = 7
        end
        du = get_du(integrator)

        
        if (integrator.u[1] > 0)
            t1 = abs.(integrator.u[1] ./ du[1])
        else
            t1 = 1e100
        end
        if (integrator.u[2] > 0)
            t2 = abs.(integrator.u[2] ./ du[2])
        else
            t2 = 1e100
        end
        if solve_n4
            if (integrator.u[3] > 0.0)
                t5 = abs.(integrator.u[3] ./ du[3])
            else
                t5 = 1e100
            end
            if (integrator.u[4] > 0.0)
                t6 = abs.(integrator.u[4] ./ du[4])
            else
                t6 = 1e100
            end
            if (integrator.u[5] > 0.0)
                t7 = abs.(integrator.u[5] ./ du[5])
            else
                t7 = 1e100
            end
        end

        if u1_eq || u2_eq
            
            
            t1 = 1e100
            integrator.u[1] = u1_fix
            t2 = 1e100
            integrator.u[2] = u2_fix
            if u2_eq
                t5 = 1e100
                t6 = 1e100
                t7 = 1e100
                integrator.u[3] = u3_fix
                integrator.u[4] = u5_fix
                integrator.u[5] = u6_fix
            end
           
            du[spinI] = - SR211 .* integrator.u[1] ./ mu .-  2 .* SR322 .* integrator.u[2] ./ mu
            if solve_n4
                du[spinI] += - SR411 .* integrator.u[3] ./ mu - SR422 .* integrator.u[4] ./ mu - SR433 .* integrator.u[5] ./ mu
            end
        end
        
        t3 = abs.(integrator.u[spinI] ./ du[spinI])
        if !solve_n4
            tcheck = minimum([t1 t2 t3])
        else
            tcheck = minimum([t1 t2 t5 t6 t7])
        end

        
        if (integrator.dt ./ tcheck .>= 1e-1)
            set_proposed_dt!(integrator, integrator.dt .* 0.9)
        elseif (integrator.dt .<= 1e-10)
            print("time step too small!! \n")
            terminate!(integrator)
        elseif (wait % 100 == 0)
            set_proposed_dt!(integrator, integrator.dt .* 1.02)
        end
    end
    
    function check_eq(u, t, integrator)
        
        if !solve_n4
            spinI = 3
            massI = 4
        else
            spinI = 6
            massI = 7
        end
        
        alph = GNew .* u[massI] .* mu
        du = get_du(integrator)
        # eq_threshold = 1e-4
        
        
        if u2_kill
            print("HERE ???? \n\n\n")
            terminate!(integrator)
        end
        if u2_eq
            return true
        end

        if u1_eq
            # check if n=4 can change equilibrium!
            if !solve_n4
                return true
            else
                # check if anything is perturbing equilibrium....
                rP = 1.0 .+ sqrt.(1 - u[spinI].^2)
                hold_n1 = -2.5e-8 .* alph.^11 .* (M_pl ./ fa).^4 .* rP .* u[1] .* u[2] .* u[3]
                hold_n1 += 3.8e-9 .* (alph.^8) .* (M_pl ./ fa).^4 .* u[1] .* u[2] .* u[3]
                hold_n1 += 9.2e-11 .* alph.^8 .* (M_pl ./ fa).^4  .* u[1] .* u[5].^2 .* mu .* integrator.dt ./ hbar .* 3.15e7
                hold_n1 += 2.6e-9 .* alph.^8 .* (M_pl ./ fa).^4  .* u[2] .* u[5] .* u[3]
                hold_n1 += -9.1e-8 .* alph.^11 .* (M_pl ./ fa).^4  .* u[1]  .* u[2] .* u[5]
                hold_n1 += - 2 * 1.5e-7 .* alph.^11 .* (M_pl ./ fa).^4 .* rP .* u[1].^2 .* u[4]
                hold_n1 += -1.1e-9 .* alph.^7 .* (M_pl ./ fa).^4 .* rP .* u[1] .* u[4] .* u[5]
                hold_n1 += 5e-12 .* alph.^8 .* (M_pl ./ fa).^4 .* u[4] .* u[3] .* u[1]
                hold_n1 *= mu .* integrator.dt ./ hbar .* 3.15e7
                
                 
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
    

        k2I_333 = 1e-8 # k^211xInfinity_322x322
        kGW_3t2 = 5e-6 # k^GW_{322->211}
        GR_211 = SR211 .* integrator.u[1] .+ k2I_333 .* alph.^8 .* (M_pl ./ fa).^4 .* integrator.u[2].^2 .* integrator.u[1] .* mu
        GR_211 += kGW_3t2 .* alph.^10 .* integrator.u[1] .* integrator.u[2] .* mu
        if solve_n4
            GR_211 += 2.6e-9 .* alph.^8 .* (M_pl ./ fa).^4  .* integrator.u[2] .* integrator.u[5] .* integrator.u[3] .* mu
            GR_211 += 9.2e-11 .* alph.^8 .* (M_pl ./ fa).^4  .* integrator.u[1] .* integrator.u[5].^2 .* mu
            GR_211 += 5e-12 .* alph.^8 .* (M_pl ./ fa).^4 .* integrator.u[4] .* integrator.u[3] .* integrator.u[1] .* mu
        end
        
            
            
        if solve_n4
            GR_411 = SR411 .* integrator.u[3]
            
            if (GR_411 > 0)&&(integrator.u[3] > 0)
                stable411 = abs.(du[3] ./ (GR_411 ./ hbar .* 3.15e7))
            else
                stable411 = 0.0
            end
            
            GR_422 = SR422 .* integrator.u[4] + 1.5e-7 .* alph.^11 .* (M_pl ./ fa).^4 .* rP .* integrator.u[1].^2 .* integrator.u[4] * mu # 211 x 211 -> 422 x BH
            GR_422 += 2.3e-7 .* alph.^7 .* (M_pl ./ fa).^4 .* rP .* integrator.u[3].^2 .* integrator.u[4] .* mu # 411 x 411 -> 422 x BH
            if (GR_422 > 0)&&(integrator.u[4] > 0)
                stable422 = abs.(du[4] ./ (GR_422 ./ hbar .* 3.15e7))
            else
                stable422 = 0.0
            end
            
            GR_433 = SR433 .* integrator.u[5] + 9.1e-8 .* alph.^11 .* (M_pl ./ fa).^4  .* integrator.u[1]  .* integrator.u[2] .* integrator.u[5] .* mu
            GR_433 += 1.1e-9 .* alph.^7 .* (M_pl ./ fa).^4 .* rP .* integrator.u[1] .* integrator.u[4] .* integrator.u[5] .* mu
            
            if (GR_433 > 0)&&(integrator.u[5] > 0)
                stable433 = abs.(du[5] ./ (GR_433 ./ hbar .* 3.15e7))
            else
                stable433 = 0.0
            end
        end
        
        if GR_211 > 0
            stable211 = abs.(du[1] ./ (GR_211 ./ hbar .* 3.15e7))
        else
            stable211 = 0.0
        end
        if GR_322 > 0
            stable322 = abs.(du[2] ./ (GR_322 ./ hbar .* 3.15e7))
        else
            stable322 = 0.0
        end
    
        
        if (stable211 .< eq_threshold)&&(stable322 .< eq_threshold)

            u1_eq = true
            u1_fix = integrator.u[1]
            u2_fix = integrator.u[2]
            u4_fix = integrator.u[massI]
            if solve_n4
                if (stable411 .< eq_threshold)&&(stable422 .< eq_threshold)&&(stable433 .< eq_threshold)
                    u2_eq = true
                    u3_fix = integrator.u[3]
                    u5_fix = integrator.u[4]
                    u6_fix = integrator.u[5]
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
            spinI = 6
            massI = 7
        end
        du = get_du(integrator)
        
        if u1_eq && u2_eq
            du[1] = 0.0
            du[2] = 0.0
            du[massI] = 0.0
            if !solve_n4
                set_u!(integrator, [u1_fix, u2_fix, integrator.u[spinI], u4_fix])
            else
                set_u!(integrator, [u1_fix, u2_fix, u3_fix, u5_fix, u6_fix,  integrator.u[spinI], u4_fix])
            end
        elseif u1_eq
            du[1] = 0.0
            du[2] = 0.0
            du[massI] = 0.0
            if !solve_n4
                set_u!(integrator, [u1_fix, u2_fix, integrator.u[spinI], u4_fix])
            else
                set_u!(integrator, [u1_fix, u2_fix, integrator.u[3], integrator.u[4], integrator.u[5], integrator.u[spinI], u4_fix])
            end
        end
        
    end
    

    function check_spin(u, t, integrator)
        if !solve_n4
            BHs_idx = 3
        else
            BHs_idx = 6
        end
        if u[BHs_idx] <= stop_on_a
            return true
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
            BHs_idx = 6
        end
        if integrator.u[BHs_idx] <= stop_on_a
            terminate!(integrator)
        end
        if integrator.u[BHs_idx] .> aBH
            integrator.u[BHs_idx] = aBH
        elseif integrator.u[BHs_idx] .< 0.0
            integrator.u[BHs_idx] = 0.0
        end
        set_proposed_dt!(integrator, integrator.dt .* 0.3)
    end
    
    function check_terminate_lvl(u, t, integrator)
        if solve_n4
            if ((u[3] < 1e-100)||(u[4] < 1e-100)||(u[5] < 1e-100)||(u[2] < 1e-100)||(u[1] < 1e-100))
                return true
            else
                return false
            end
        else
            return false
        end

    end
    function affect_terminate_lvl!(integrator)
        if integrator.u[1] < 1e-100
            integrator.u[1] = 0.0
            du = get_du(integrator)
            du[3] = 0.0
        end
        if integrator.u[2] < 1e-100
            integrator.u[2] = 0.0
            du = get_du(integrator)
            du[3] = 0.0
        end
        if integrator.u[3] < 1e-100
            integrator.u[3] = 0.0
            du = get_du(integrator)
            du[3] = 0.0
        end
        if integrator.u[4] < 1e-100
            integrator.u[4] = 0.0
            du = get_du(integrator)
            du[4] = 0.0
        end
        if integrator.u[5] < 1e-100
            integrator.u[5] = 0.0
            du = get_du(integrator)
            du[5] = 0.0
        end

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
    cback_term = DiscreteCallback(check_terminate_lvl, affect_terminate_lvl!, save_positions=(false, true))
    
    if solve_n4
    
        cbset = CallbackSet(cbackdt, cbackspin, cback_term)
        
        if (GNew .* M_BH .* mu ) < 0.1
            abstol = 1e-30
        else
            abstol = 1e-15
        end
        prob = ODEProblem(RHS_ax!, y0, tspan, Mvars, reltol=1e-3, abstol=abstol) #
        # sol = solve(prob, AutoTsit5(Rodas4()), dt=dt_guess, saveat=saveat, callback=cbset)
        
        
        # prob = ODEProblem(RHS_ax!, y0, tspan, Mvars, reltol=1e-5, abstol=1e-15)
        # prob = ODEProblem(RHS_ax!, y0, tspan, Mvars, reltol=1e-12, abstol=1e-20)
        sol = solve(prob, Rosenbrock23(), dt=dt_guess, saveat=saveat, callback=cbset, maxiters=5e5)
        # sol = solve(prob, KenCarp5(), dt=dt_guess, saveat=saveat, callback=cbset)
        # sol = solve(prob, Rodas4(), dt=dt_guess, saveat=saveat, callback=cbset)
        # sol = solve(prob, Euler(), dt=dt_guess, saveat=saveat, callback=cbset)
    else
        
        abstol = 1e-30
        
        cbset = CallbackSet(cback_equil, cbackdt, cbackspin)
        # prob = ODEProblem(RHS_ax!, y0, tspan, Mvars, reltol=1e-5, abstol=1e-15)
        # sol = solve(prob, Rodas4(), dt=dt_guess, saveat=saveat, callback=cbset)
        
        prob = ODEProblem(RHS_ax!, y0, tspan, Mvars, reltol=1e-3, abstol=abstol)
        # sol = solve(prob, Euler(), dt=dt_guess, saveat=saveat, callback=cbset)
        sol = solve(prob, Rosenbrock23(), dt=dt_guess, saveat=saveat, callback=cbset)
        # sol = solve(prob, Rodas4(), dt=dt_guess, saveat=saveat, callback=cbset)
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
    return spinBH[end]
 
end
