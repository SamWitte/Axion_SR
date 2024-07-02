using Random
using OrdinaryDiffEq
using Statistics
using Distributions
using DelimitedFiles
using Interpolations
include("Constants.jl")
include("solve_sr_rates.jl")
include("load_rates.jl")
using Printf

function super_rad_check(M_BH, aBH, massB, f_a; spin=0, tau_max=1e4, alpha_max_cut=10.0, debug=false, solve_322=true, solve_n4=false, solve_n5=false, impose_low_cut=0.01, input_data="Masha", stop_on_a=0, eq_threshold=1e-100, abstol=1e-30, non_rel=true)
   
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
    
    
    if input_data == "Masha"
        OmegaH = aBH ./ (2 .* (GNew .* M_BH) .* (1 .+ sqrt.(1 .- aBH.^2)))
        if (ergL(2, 1, 1, massB, M_BH, aBH) .>= OmegaH)&&(f_a .< 2.6e17)
            return aBH, M_BH
        end
    end
    
    final_spin, final_BH = solve_system(massB, f_a, aBH, M_BH, tau_max, debug=debug, solve_322=solve_322, impose_low_cut=impose_low_cut, input_data=input_data, solve_n4=solve_n4, solve_n5=solve_n5, stop_on_a=stop_on_a, eq_threshold=eq_threshold, abstol=abstol, non_rel=non_rel)
    
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

function solve_system(mu, fa, aBH, M_BH, t_max; n_times=10000, debug=true, solve_322=true, impose_low_cut=0.01, return_all_info=false, input_data="Masha", solve_n4=false, solve_n5=false, eq_threshold=1e-4, stop_on_a=0, abstol=1e-30, non_rel=true)

    if !solve_n4
        idx_lvl = 2 # number of states
    else
        if !solve_n5
            idx_lvl = 5
        else
            idx_lvl = 8
        end
    end
    spinI = idx_lvl + 1
    massI = spinI + 1
    
    
    e_init = 1.0 ./ (GNew .* M_BH.^2 .* M_to_eV) # unitless
    default_reltol = 5e-3
    y0 = []
    reltol = []
    for i in 1:idx_lvl
        append!(y0, e_init)
        append!(reltol, default_reltol)
    end
    append!(y0, aBH)
    append!(y0, M_BH)
    y0 = log.(y0)
    
    append!(reltol, default_reltol)
    append!(reltol, default_reltol)
    
    wait = 0 # tracker
    alph = GNew .* M_BH .* mu
    
    
    e2_maxBN = 1024 * pi * (fa / M_pl).^2 ./ (9 * alph.^3)
    e3_maxBN = e2_maxBN .* (3 ./ 2).^4
    e4_maxBN = e2_maxBN .* (4 ./ 2).^4
    e5_maxBN = e2_maxBN .* (5 ./ 2).^4
    if debug
        print("BN Threshold \t ", e2_maxBN, "\t", e3_maxBN, "\t", e4_maxBN, "\t", e5_maxBN, "\n")
    end
    
    
    Emax2 = 1.0
    OmegaH = aBH ./ (2 .* (GNew .* M_BH) .* (1 .+ sqrt.(1 .- aBH.^2)))
    if (OmegaH .> ergL(2, 1, 1, mu, M_BH, aBH))
        Emax2 = emax_211(M_BH, mu, aBH)
    end
    
    
    Mvars = [mu, fa, Emax2, solve_322, aBH, M_BH, impose_low_cut]
    tspan = (0.0, t_max)
    saveat = (tspan[2] .- tspan[1]) ./ n_times
    
    
    N_pts_interp = 10
    xtol_slv = 1e-15
    iter_slv = 50
    
    SR_rates = zeros(idx_lvl)
    m_list = [1 2 1 2 3 2 3 4] # [211, 322, 411, 422, 433, 522, 533, 544]
    bn_list = [e2_maxBN e3_maxBN e4_maxBN e4_maxBN e4_maxBN e5_maxBN e5_maxBN e5_maxBN]
    
    
    n = 2
    l = 1
    m = 1
    amin_guess = 8 * m * n.^2 .* alph .* (2 .* n.^2 .+ alph.^2) ./ (16 .* n.^4 .* alph.^2 .+ m.^2 .* (2 .* n.^2 .+ alph.^2).^2)
    alist, pts = compute_gridded(mu, M_BH, aBH, n, l, m; iter=iter_slv, xtol=xtol_slv, npts=N_pts_interp, amin=(amin_guess .* 0.99))
    condit = pts .> 1e-50
    if sum(condit) > 2
        itp_211 = LinearInterpolation(alist[condit], log10.(pts[condit]), extrapolation_bc=-100)
    else
        itp_211 = LinearInterpolation(alist, -100 .* ones(length(alist)), extrapolation_bc=-100)
    end
    SR_rates[1] = 10 .^itp_211(aBH)
    
    n = 3
    l = 2
    m = 2
    amin_guess = 8 * m * n.^2 .* alph .* (2 .* n.^2 .+ alph.^2) ./ (16 .* n.^4 .* alph.^2 .+ m.^2 .* (2 .* n.^2 .+ alph.^2).^2)
    alist, pts = compute_gridded(mu, M_BH, aBH, n, l, m; iter=iter_slv, xtol=xtol_slv, npts=N_pts_interp, amin=(amin_guess .* 0.99))
    condit = pts .> 1e-50
    if sum(condit) > 2
        itp_322 = LinearInterpolation(alist[condit], log10.(pts[condit]), extrapolation_bc=-100)
    else
        itp_322 = LinearInterpolation(alist, -100 .* ones(length(alist)), extrapolation_bc=-100)
    end
    SR_rates[2] = 10 .^itp_322(aBH)
    
    if solve_n4
        n = 4
        l = 1
        m = 1
        amin_guess = 8 * m * n.^2 .* alph .* (2 .* n.^2 .+ alph.^2) ./ (16 .* n.^4 .* alph.^2 .+ m.^2 .* (2 .* n.^2 .+ alph.^2).^2)
        alist, pts = compute_gridded(mu, M_BH, aBH, n, l, m; iter=iter_slv, xtol=xtol_slv, npts=N_pts_interp, amin=(amin_guess .* 0.99))
        condit = pts .> 1e-50
        if sum(condit) > 2
            itp_411 = LinearInterpolation(alist[condit], log10.(pts[condit]), extrapolation_bc=-100)
        else
            itp_411 = LinearInterpolation(alist, -100 .* ones(length(alist)), extrapolation_bc=-100)
        end
        SR_rates[3] = 10 .^itp_411(aBH)
        
        n = 4
        l = 2
        m = 2
        amin_guess = 8 * m * n.^2 .* alph .* (2 .* n.^2 .+ alph.^2) ./ (16 .* n.^4 .* alph.^2 .+ m.^2 .* (2 .* n.^2 .+ alph.^2).^2)
        alist, pts = compute_gridded(mu, M_BH, aBH, n, l, m; iter=iter_slv, xtol=xtol_slv, npts=N_pts_interp, amin=(amin_guess .* 0.99))
        condit = pts .> 1e-50
        if sum(condit) > 2
            itp_422 = LinearInterpolation(alist[condit], log10.(pts[condit]), extrapolation_bc=-100)
        else
            itp_422 = LinearInterpolation(alist, -100 .* ones(length(alist)), extrapolation_bc=-100)
        end
        SR_rates[4] = 10 .^itp_422(aBH)
        
        
        n = 4
        l = 3
        m = 3
        amin_guess = 8 * m * n.^2 .* alph .* (2 .* n.^2 .+ alph.^2) ./ (16 .* n.^4 .* alph.^2 .+ m.^2 .* (2 .* n.^2 .+ alph.^2).^2)
        alist, pts = compute_gridded(mu, M_BH, aBH, n, l, m; iter=iter_slv, xtol=xtol_slv, npts=N_pts_interp, amin=(amin_guess .* 0.99))
        condit = pts .> 1e-50
        if sum(condit) > 2
            itp_433 = LinearInterpolation(alist[condit], log10.(pts[condit]), extrapolation_bc=-100)
        else
            itp_433 = LinearInterpolation(alist, -100 .* ones(length(alist)), extrapolation_bc=-100)
        end
        SR_rates[5] = 10 .^itp_433(aBH)
        
        if solve_n5

            n = 5
            l = 2
            m = 2
            amin_guess = 8 * m * n.^2 .* alph .* (2 .* n.^2 .+ alph.^2) ./ (16 .* n.^4 .* alph.^2 .+ m.^2 .* (2 .* n.^2 .+ alph.^2).^2)
            alist, pts = compute_gridded(mu, M_BH, aBH, n, l, m; iter=iter_slv, xtol=xtol_slv, npts=N_pts_interp, amin=(amin_guess .* 0.99))
            condit = pts .> 1e-50
            if sum(condit) > 2
                itp_522 = LinearInterpolation(alist[condit], log10.(pts[condit]), extrapolation_bc=-100)
            else
                itp_522 = LinearInterpolation(alist, -100 .* ones(length(alist)), extrapolation_bc=-100)
            end
            SR_rates[6] = 10 .^itp_522(aBH)
        
            n = 5
            l = 3
            m = 3
            amin_guess = 8 * m * n.^2 .* alph .* (2 .* n.^2 .+ alph.^2) ./ (16 .* n.^4 .* alph.^2 .+ m.^2 .* (2 .* n.^2 .+ alph.^2).^2)
            alist, pts = compute_gridded(mu, M_BH, aBH, n, l, m; iter=iter_slv, xtol=xtol_slv, npts=N_pts_interp, amin=(amin_guess .* 0.99))
            condit = pts .> 1e-50
            if sum(condit) > 2
                itp_533 = LinearInterpolation(alist[condit], log10.(pts[condit]), extrapolation_bc=-100)
            else
                itp_533 = LinearInterpolation(alist, -100 .* ones(length(alist)), extrapolation_bc=-100)
            end
            SR_rates[7] = 10 .^itp_533(aBH)
        
            n = 5
            l = 4
            m = 4
            amin_guess = 8 * m * n.^2 .* alph .* (2 .* n.^2 .+ alph.^2) ./ (16 .* n.^4 .* alph.^2 .+ m.^2 .* (2 .* n.^2 .+ alph.^2).^2)
            alist, pts = compute_gridded(mu, M_BH, aBH, n, l, m; iter=iter_slv, xtol=xtol_slv, npts=N_pts_interp, amin=(amin_guess .* 0.99))
            condit = pts .> 1e-50
            if sum(condit) > 2
                itp_544 = LinearInterpolation(alist[condit], log10.(pts[condit]), extrapolation_bc=-100)
            else
                itp_544 = LinearInterpolation(alist, -100 .* ones(length(alist)), extrapolation_bc=-100)
            end
            SR_rates[8] = 10 .^itp_544(aBH)
        
        end
    end
    
    if debug
        print("SR rates @ prod \t", SR_rates, "\n")
    end
    
    rates = load_rate_coeffs(mu, M_BH, aBH, fa; non_rel=non_rel, input_data=input_data, solve_n4=solve_n4, solve_n5=solve_n5)
    turn_off = []
    for i in 1:idx_lvl
        append!(turn_off, false)
    end

    #### DEFINING FUNCTION FOR EVOLUTION ######
    function RHS_ax!(du, u, Mvars, t)
    
        # [e211, e322, aBH, MBH] or [e211, e322, e411, aBH, MBH]....
        # mu [mass eV], fa[1/GeV]
        u_real = exp.(u)
        
        mu, fa, Emax2, solve_322, aBH_i, M_BH_i, impose_low_cut  = Mvars

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
        
        SR_rates[1] = 10 .^itp_211(u_real[spinI])
        if (OmegaH .< ergL(2, 1, 1, mu, u_real[massI], u_real[spinI]))
            SR_rates[1] *= 0.0
        end

        
        SR_rates[2] = 10 .^itp_322(u_real[spinI])
        if (2 .* OmegaH .< ergL(3, 2, 2, mu, u_real[massI], u_real[spinI]))
            SR_rates[2] *= 0.0
        end
        
        if solve_n4
            SR_rates[3] = 10 .^itp_411(u_real[spinI])
            if (OmegaH .< ergL(4, 1, 1, mu, u_real[massI], u_real[spinI]))
                SR_rates[3] *= 0.0
            end
            
            SR_rates[4] = 10 .^itp_422(u_real[spinI])
            if (2 .* OmegaH .< ergL(4, 2, 2, mu, u_real[massI], u_real[spinI]))
                SR_rates[4] *= 0.0
            end
            
            SR_rates[5] = 10 .^itp_433(u_real[spinI])
            if (3 .* OmegaH .< ergL(4, 3, 3, mu, u_real[massI], u_real[spinI]))
                SR_rates[5] *= 0.0
            end
            
            if solve_n5
                SR_rates[6] = 10 .^itp_522(u_real[spinI])
                if (2 .* OmegaH .< ergL(5, 2, 2, mu, u_real[massI], u_real[spinI]))
                    SR_rates[6] *= 0.0
                end
                
                SR_rates[7] = 10 .^itp_533(u_real[spinI])
                if (3 .* OmegaH .< ergL(5, 3, 3, mu, u_real[massI], u_real[spinI]))
                    SR_rates[7] *= 0.0
                end
                
                SR_rates[8] = 10 .^itp_544(u_real[spinI])
                if (4 .* OmegaH .< ergL(5, 4, 4, mu, u_real[massI], u_real[spinI]))
                    SR_rates[8] *= 0.0
                end
            end
        end
            
        
        if u_real[1] .> Emax2
            SR_rates[1] *= 0.0
        end
        
    
    
        # SR terms
        du[spinI] = 0.0
        du[massI] = 0.0
        for i in 1:idx_lvl
            if u_real[i] < e_init
                du[i] = SR_rates[i] .* e_init ./ mu
                du[spinI] += - m_list[i] * SR_rates[i] .* e_init ./ mu
                du[massI] += - SR_rates[i] .* e_init ./ mu
            else
                du[i] = SR_rates[i] .* u_real[i] ./ mu
                du[spinI] += - m_list[i] * SR_rates[i] .*  u_real[i] ./ mu
                du[massI] += - SR_rates[i] .*  u_real[i] ./ mu
            end
        end
        
        # Scattering terms
        rate_keys = collect(keys(rates))
        for i in 1:length(rate_keys)
            idxV, sgn = key_to_indx(rate_keys[i]; solve_n4=solve_n4, solve_n5=solve_n5)
            u_term_tot = 1.0
            
            for j in 1:length(sgn)
                if (idxV[j] <= idx_lvl)&&(idxV[j] > 0)
                    if j == 3
                        u_term_tot *= (u_real[idxV[j]] .+ e_init)
                    else
                        u_term_tot *=  u_real[idxV[j]]
                    end
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
            if (sign_flip[i] == true)
                du[i] *= 0.0
            elseif ((abs.(u[i] .- log.(bn_list[i])) < 1e-2)||(u[i] > log.(bn_list[i])))&&(du[i] > 0)
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
    
    
    #### CALLBACK 1
    function check_oscillating(u, t, integrator)
    
        if length(integrator.sol.t) < 3000
            return false
        else
            for i in 1:idx_lvl
                if u[i] < log.(1e-60)
                    state = [integrator.sol.u[j][i] for j in 1:length(integrator.sol.u)]
                    short = round.(state[end-1000:end], digits=1)
                    cnt_extrm = 0
                    cnt_extrm += sum(abs.(diff(sign.(short[2:end] .- short[1:end-1])))) / 2
                    if cnt_extrm > 10
                        return true
                    end
                    
                end
            end
            return false
        end
    end
    
    function affect_oscillating!(integrator)
        for i in 1:idx_lvl
            state = [integrator.sol.u[j][i] for j in 1:length(integrator.sol.u)]
            
            short = round.(state[end-1000:end], digits=1)
            cnt_extrm = 0
            
            cnt_extrm += sum(abs.(diff(sign.(short[2:end] .- short[1:end-1])))) / 2
            if (cnt_extrm > 10)&&(integrator.u[i] < log.(1e-60))&&(turn_off[i] == false)
                turn_off[i] = true
                integrator.u[i] = log.(e_init)
#            else
#                print(i, "\t", turn_off[i], "\n")
#                turn_off[i] = false
            end
        end
    end
    
    sign_flip = []
    for i in 1:idx_lvl
        append!(sign_flip, false)
    end
    
    #### CALLBACK 2
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
        
        
        for i in 1:length(rate_keys)
            idxV, sgn = key_to_indx(rate_keys[i]; solve_n4=solve_n4, solve_n5=solve_n5)
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
        
        for i in 1:idx_lvl
            if (sign(all_contribs[i]) == sign(test[i]))
                sign_flip[i] = false
                reltol[i] = default_reltol
            else
                sign_flip[i] = true
                # print(i, " SGN FLIP \t", sign_flip, "\n")
                # reltol[i] = 1e-3
            end
        end
        
        integrator.opts.reltol =  reltol
        
        tlist = []
        # print(integrator.sol.u, "\n")
        for i in 1:idx_lvl
            condBN =  (abs.(u[i] .- log.(bn_list[i])) < 1e-2)
            
            if (u[i] > log.(e_init)) && condBN && (sign_flip[i]==false)
                append!(tlist, abs.(1.0 ./ du[i]))
            end
        end
            
        append!(tlist, 0.1 * 1.0 ./ du[spinI])
        tmin = minimum(abs.(tlist))
       
        # print("Tmin \t", integrator.dt, "\t", tlist, "\n")
        if (integrator.dt ./ tmin .>= 0.1)
            return true
        elseif (integrator.dt ./ tmin .<= 0.001)
            return true
        elseif (integrator.dt .<= 1e-10)
            return true
        elseif (integrator.dt ./ integrator.t .<= 1e-4)
            integrator.opts.dtmin = integrator.t * 1e-4
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
            if (integrator.u[i] > log.(e_init)) && condBN && (sign_flip[i]==false)
                append!(tlist, (1.0 ./ du[i]))
                append!(indx_list, i)
            end
        end
            
        append!(tlist, 0.1 * 1.0 ./ du[spinI])
        tmin = minimum(abs.(tlist))
        
        
        if (integrator.dt ./ integrator.t < 1e-6)
            for i in 1:idx_lvl
                if reltol[i] < 0.1
                    reltol[i] *= 2.0
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
        elseif (integrator.dt ./ tmin .>= 0.1)
            set_proposed_dt!(integrator, tmin .* 0.4)
        elseif (integrator.dt ./ tmin .<= 1e-3)&&(wait % 100 == 0)
            set_proposed_dt!(integrator, integrator.dt .* 1.01)
        elseif (integrator.dt ./ tmin .<= 1e-6)||(integrator.dt ./ integrator.t .<= 1e-4)
            for i in 1:idx_lvl
                if reltol[i] < 0.1
                    reltol[i] *= 2.0
                    integrator.opts.reltol =  reltol
                    
                else
                    if integrator.opts.abstol < 1e-10
                        integrator.opts.abstol *= 2.0
                    end
                end
            end
        
        elseif (integrator.dt .<= 1e-10)
            print("time step too small!! \n")
            terminate!(integrator)
        end
    end
    
    # Callback 3
    function check_spin(u, t, integrator)
        u_real = exp.(u)
        if debug && (wait%1==0)
            if solve_n4
                if !solve_n5
                    @printf("Time and Vals: \t %.7f  %.3e  %.3e  %.3e  %.3e  %.3e  %.3e  %.3e \n", t, u_real[1], u_real[2], u_real[3], u_real[4], u_real[5], u_real[6], u_real[7])
                    # print(integrator.dt ./ t, "\n")
                    if (t > 4200)&&(integrator.opts.reltol[1] < 0.1)
                        integrator.opts.reltol[1] *= 2.0
#                        print(integrator.opts.reltol, "\n")
#                        print(integrator.opts.abstol, "\n")
                    end
                else
                    @printf("Time and Vals: \t %.7f  %.3e  %.3e  %.3e  %.3e  %.3e  %.3e  %.3e  %.3e  %.3e  %.3e \n", t, u_real[1], u_real[2], u_real[3], u_real[4], u_real[5], u_real[6], u_real[7], u_real[8], u_real[9], u_real[10])
                end
                # du = get_du(integrator)
                # print(integrator.dt, "\t", du, "\t", "\n\n")
            else
                print(t, "\t", u_real[1], "\t", u_real[2], "\t", u_real[3], "\t", u_real[4], "\n")
            
            end
        end
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
    

    
    
    rP = 1.0 .+ sqrt.(1 - aBH.^2)
    t1 = (SR_rates[1] ./ hbar .* 3.15e7).^(-1)
    t2 = (SR_rates[2] ./ hbar .* 3.15e7).^(-1)
    if solve_n4
        t3 = (SR_rates[4] ./ hbar .* 3.15e7).^(-1)
        t4 = (SR_rates[5] ./ hbar .* 3.15e7).^(-1)
        tlist = [t1 t2 t3 t4]
    else
        tlist = [t1 t2]
    end
    dt_guess = minimum(tlist)
    if debug
        print("Time guess \t", dt_guess, "\n")
    end

    cback_osc = DiscreteCallback(check_oscillating, affect_oscillating!, save_positions=(false, true))
    cbackdt = DiscreteCallback(check_timescale, affect_timescale!, save_positions=(false, true))
    cbackspin = DiscreteCallback(check_spin, affect_spin!, save_positions=(false, true))
    cbset = CallbackSet(cbackspin, cbackdt, cback_osc)
    
    if solve_n4
        prob = ODEProblem(RHS_ax!, y0, tspan, Mvars, reltol=reltol, abstol=1e-10)
        # sol = solve(prob, Rosenbrock23(autodiff=true), dt=dt_guess, saveat=saveat, callback=cbset, maxiters=5e6, dtmin=(dt_guess / 1e5), force_dtmin=true)
        sol = solve(prob, Rosenbrock23(autodiff=false), dt=dt_guess, saveat=saveat, callback=cbset, maxiters=5e6, dtmin=(dt_guess / 1e5), force_dtmin=true)
        #sol = solve(prob, Euler(), dt=dt_guess, saveat=saveat, callback=cbset)
    else
        prob = ODEProblem(RHS_ax!, y0, tspan, Mvars,  reltol=reltol, abstol=1e-6)
        sol = solve(prob, Rosenbrock23(autodiff=false), dt=dt_guess, saveat=saveat, callback=cbset, maxiters=5e6)
    end
    
    
    
    state211 = [exp.(sol.u[i][1]) for i in 1:length(sol.u)]
    state322 = [exp.(sol.u[i][2]) for i in 1:length(sol.u)]
    if solve_n4
        state411 = [exp.(sol.u[i][3]) for i in 1:length(sol.u)]
        state422 = [exp.(sol.u[i][4]) for i in 1:length(sol.u)]
        state433 = [exp.(sol.u[i][5]) for i in 1:length(sol.u)]
        if solve_n5
            state522 = [exp.(sol.u[i][6]) for i in 1:length(sol.u)]
            state533 = [exp.(sol.u[i][7]) for i in 1:length(sol.u)]
            state544 = [exp.(sol.u[i][8]) for i in 1:length(sol.u)]
        end
    end
    spinBH = [exp.(sol.u[i][spinI]) for i in 1:length(sol.u)]
    MassB = [exp.(sol.u[i][massI]) for i in 1:length(sol.u)]
    

    if debug
        if !solve_n4
            print("Initial \t", state211[1], "\t", state322[1], "\t", spinBH[1], "\t", MassB[1], "\n")
            print("Final \t", state211[end], "\t", state322[end], "\t", spinBH[end], "\t", MassB[end], "\n")
        else
            if !solve_n5
                print("Initial \t", state211[1], "\t", state322[1] , "\t",  state411[1] , "\t", state422[1] , "\t", state433[1] , "\t", spinBH[1], "\t", MassB[1], "\n")
                print("Final \t", state211[end] , "\t", state322[end],"\t", state411[end],"\t", state422[end] , "\t",  state433[end] , "\t", spinBH[end], "\t", MassB[end], "\n")
            else
                print("Initial \t", state322[1], "\t", state422[1] , "\t", state433[1], "\t", state522[1], "\t", state533[1] , "\t", state544[1], "\t", spinBH[1], "\t", MassB[1], "\n")
                print("Final \t", state322[end], "\t", state422[end] , "\t", state433[end], "\t", state522[end], "\t", state533[end] , "\t", state544[end], "\t", spinBH[end], "\t", MassB[end], "\n")
            end
        end
    end

    if return_all_info
        if !solve_n4
            return sol.t, state211, state322, spinBH, MassB
        else
            if !solve_n5
                return sol.t, state211, state322, state411, state422, state433, spinBH, MassB
            else
                return sol.t, state211, state322, state411, state422, state433, state522, state533, state544, spinBH, MassB
            end
        end
    end
    
    if (sol.t[end] != t_max)&&(spinBH[end] > stop_on_a)
        print("Fail? Final time only \t", sol.t[end], "\n")
    end
    if isnan(spinBH[end])
        spinBH = spinBH[.!isnan.(spinBH)]
    end
    return spinBH[end], MassB[end]
 
end
