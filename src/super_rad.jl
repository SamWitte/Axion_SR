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

function super_rad_check(M_BH, aBH, massB, f_a; spin=0, tau_max=1e4, alpha_max_cut=100.0, debug=false, solve_322=true, solve_n4=false, solve_n5=false, impose_low_cut=0.01, input_data="Masha", stop_on_a=0, eq_threshold=1e-100, abstol=1e-30, non_rel=true, high_p=true, N_pts_interp=15, N_pts_interpL=10)
   
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
    print("Solving system... \n")
    if debug || (alph > 0.1)
        final_spin, final_BH = solve_system(massB, f_a, aBH, M_BH, tau_max, debug=debug, solve_322=solve_322, impose_low_cut=impose_low_cut, input_data=input_data, solve_n4=solve_n4, solve_n5=solve_n5, stop_on_a=stop_on_a, eq_threshold=eq_threshold, abstol=abstol, non_rel=non_rel, high_p=high_p, N_pts_interp=N_pts_interp, N_pts_interpL=N_pts_interpL)
    else
        # n = 4,5 don't alter evo for small alpha. accelerate.
        final_spin, final_BH = solve_system(massB, f_a, aBH, M_BH, tau_max, debug=debug, solve_322=solve_322, impose_low_cut=impose_low_cut, input_data=input_data, solve_n4=false, solve_n5=false, stop_on_a=stop_on_a, eq_threshold=eq_threshold, abstol=abstol, non_rel=non_rel, high_p=high_p, N_pts_interp=N_pts_interp, N_pts_interpL=N_pts_interpL)
    end
    
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

function solve_system(mu, fa, aBH, M_BH, t_max; n_times=10000, debug=false, solve_322=true, impose_low_cut=0.01, return_all_info=false, input_data="Masha", solve_n4=false, solve_n5=false, eq_threshold=1e-100, stop_on_a=0, abstol=1e-30, non_rel=true, max_m_2=false, high_p=true, N_pts_interp=10, N_pts_interpL=5, trace_term=false, trace_idx=1)
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
    
    if !solve_n4
        default_reltol = 1e-4
        idx_lvl = 3 # number of states
    else

        if !solve_n5
            idx_lvl = 6
        else
            idx_lvl = 10
        end
    end
    
    ### only for testing...
    default_reltol = 1e-5
    reltol_Thres = 1e-4
    
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
    
    def_spin_tol = 1e-4
    append!(reltol, def_spin_tol)
    append!(reltol, def_spin_tol)
    
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
    
    
    
    xtol_slv = 1e-15
    iter_slv = 50
    
    SR_rates = zeros(idx_lvl)
    m_list = [1 1 2 1 2 3 1 2 3 4] # [211, 311, 322, 411, 422, 433, 511, 522, 533, 544]
    bn_list = [e2_maxBN e3_maxBN e3_maxBN e4_maxBN e4_maxBN e4_maxBN e5_maxBN e5_maxBN e5_maxBN e5_maxBN]
    idx_fill = 1
    delt_a = 0.0001
    
    n = 2; l = 1; m = 1;
    amin_guess_211, run_high, a_list_high, out_high, run_low, a_list_low, out_low = pre_computed_sr_rates(n, l, m, alph, M_BH; n_high=100, n_low=100, delt_a=delt_a)
    if alph < 0.5
        # amin_guess_211 = 8 * m * n.^2 .* alph .* (2 .* n.^2 .+ alph.^2) ./ (16 .* n.^4 .* alph.^2 .+ m.^2 .* (2 .* n.^2 .+ alph.^2).^2)
        # alist, pts = compute_gridded(mu, M_BH, aBH, n, l, m; iter=iter_slv, xtol=xtol_slv, npts=N_pts_interp, amin=(amin_guess_211 .* 0.99))

        if !run_high
            a_list_high = [amin_guess_211, amin_guess_211 * 1.01]
            out_high = [1e-100, 1e-100]
        elseif length(a_list_high) == 1.0
            a_list_high = [alist[1], alist[1] * 1.01]
            out_high = [out_high[1], out_high[1]]
        end
    else
        amin_guess_211 = 1.0
        a_list_high = [amin_guess_211, amin_guess_211 * 1.01]
        out_high = [1e-100, 1e-100]
    end
    itp_211U = LinearInterpolation(a_list_high, log10.(out_high), extrapolation_bc=Interpolations.Line())
    if !run_low
        a_list_low = [amin_guess_211 .* 0.99, amin_guess_211]
        out_low = [1e-100, 1e-100]
    elseif length(a_list_low) == 1.0
        a_list_low = [a_list_low[1] .* 0.99, a_list_low[1]]
        out_low = [out_low[1], out_low[1]]
    end
    itp_211L = LinearInterpolation(a_list_low, log10.(out_low), extrapolation_bc=Interpolations.Line())
    abndryU_211 = a_list_high[1]
    abndryL_211 = a_list_low[end]
    function itp_211(aspin)
        if aspin .>= abndryU_211
            return 10.0 .^itp_211U(aspin)
        elseif aspin .<= abndryL_211
            return -10.0 .^itp_211L(aspin)
        else
            return 0.0
        end
    end
    SR_rates[idx_fill] = itp_211(aBH)
    idx_fill += 1
    
    n = 3; l = 1; m = 1;
    amin_guess_311, run_high, a_list_high, out_high, run_low, a_list_low, out_low = pre_computed_sr_rates(n, l, m, alph, M_BH; n_high=100, n_low=100, delt_a=delt_a)
    if alph < 0.5
        # amin_guess_211 = 8 * m * n.^2 .* alph .* (2 .* n.^2 .+ alph.^2) ./ (16 .* n.^4 .* alph.^2 .+ m.^2 .* (2 .* n.^2 .+ alph.^2).^2)
        # alist, pts = compute_gridded(mu, M_BH, aBH, n, l, m; iter=iter_slv, xtol=xtol_slv, npts=N_pts_interp, amin=(amin_guess_211 .* 0.99))

        if !run_high
            a_list_high = [amin_guess_311, amin_guess_311 * 1.01]
            out_high = [1e-100, 1e-100]
        elseif length(a_list_high) == 1.0
            a_list_high = [alist[1], alist[1] * 1.01]
            out_high = [out_high[1], out_high[1]]
        end
    else
        amin_guess_311 = 1.0
        a_list_high = [amin_guess_311, amin_guess_311 * 1.01]
        out_high = [1e-100, 1e-100]
    end
    itp_311U = LinearInterpolation(a_list_high, log10.(out_high), extrapolation_bc=Interpolations.Line())
    if !run_low
        a_list_low = [amin_guess_311 .* 0.99, amin_guess_311]
        out_low = [1e-100, 1e-100]
    elseif length(a_list_low) == 1.0
        a_list_low = [a_list_low[1] .* 0.99, a_list_low[1]]
        out_low = [out_low[1], out_low[1]]
    end
    itp_311L = LinearInterpolation(a_list_low, log10.(out_low), extrapolation_bc=Interpolations.Line())
    abndryU_311 = a_list_high[1]
    abndryL_311 = a_list_low[end]
    function itp_311(aspin)
        if aspin .>= abndryU_311
            return 10.0 .^itp_311U(aspin)
        elseif aspin .<= abndryL_311
            return -10.0 .^itp_311L(aspin)
        else
            return 0.0
        end
    end
    SR_rates[idx_fill] = itp_311(aBH)
    idx_fill += 1
    
    n = 3; l = 2; m = 2;
    amin_guess_322, run_high, a_list_high, out_high, run_low, a_list_low, out_low = pre_computed_sr_rates(n, l, m, alph, M_BH; n_high=1000, n_low=1000, delt_a=delt_a)
    
    if alph < 1.1
        # amin_guess_322 = 8 * m * n.^2 .* alph .* (2 .* n.^2 .+ alph.^2) ./ (16 .* n.^4 .* alph.^2 .+ m.^2 .* (2 .* n.^2 .+ alph.^2).^2)
        # alist, pts = compute_gridded(mu, M_BH, aBH, n, l, m; iter=iter_slv, xtol=xtol_slv, npts=N_pts_interp, amin=(amin_guess_322 .* 0.99))
        # cond = pts .> 0.0
        if !run_high
            a_list_high = [amin_guess_322, amin_guess_322 * 1.01]
            out_high = [1e-100, 1e-100]
        elseif length(a_list_high) == 1.0
            a_list_high = [alist[1], alist[1] * 1.01]
            out_high = [out_high[1], out_high[1]]
        end
    else
        amin_guess_211 = 1.0
        a_list_high = [amin_guess_322, amin_guess_322 * 1.01]
        out_high = [1e-100, 1e-100]
    end
    itp_322U = LinearInterpolation(a_list_high, log10.(out_high), extrapolation_bc=Interpolations.Line())
    if !run_low
        a_list_low = [amin_guess_322 .* 0.99, amin_guess_322]
        out_low = [1e-100, 1e-100]
    elseif length(a_list_low) == 1.0
        a_list_low = [a_list_low[1] .* 0.99, a_list_low[1]]
        out_low = [out_low[1], out_low[1]]
    end
    itp_322L = LinearInterpolation(a_list_low, log10.(out_low), extrapolation_bc=Interpolations.Line())
    
    abndryU_322 = a_list_high[1]
    abndryL_322 = a_list_low[end]
    function itp_322(aspin)
        if aspin .>= abndryU_322
            return 10.0 .^itp_322U(aspin)
        elseif aspin .<= abndryL_322
            return -10.0 .^itp_322L(aspin)
        else
            return (-10.0 .^itp_322L(aspin) + 10.0 .^itp_322U(aspin)) ./ 2.0 ## Need buffer region for numerical stability!
        end
    end
    SR_rates[idx_fill] = itp_322(aBH)
    idx_fill += 1
  
    amin_guess_433 = nothing
    if solve_n4
        n = 4; l = 1; m = 1;
        amin_guess_411, run_high, a_list_high, out_high, run_low, a_list_low, out_low = pre_computed_sr_rates(n, l, m, alph, M_BH; n_high=100, n_low=100, delt_a=delt_a)
        if alph < 0.5
           
            if !run_high
                a_list_high = [amin_guess_411, amin_guess_411 * 1.01]
                out_high = [1e-100, 1e-100]
            elseif length(a_list_high) == 1.0
                a_list_high = [alist[1], alist[1] * 1.01]
                out_high = [out_high[1], out_high[1]]
            end
        else
            amin_guess_411 = 1.0
            a_list_high = [amin_guess_411, amin_guess_411 * 1.01]
            out_high = [1e-100, 1e-100]
        end
        itp_411U = LinearInterpolation(a_list_high, log10.(out_high), extrapolation_bc=Interpolations.Line())
        if !run_low
            a_list_low = [amin_guess_411 .* 0.99, amin_guess_411]
            out_low = [1e-100, 1e-100]
        elseif length(a_list_low) == 1.0
            a_list_low = [a_list_low[1] .* 0.99, a_list_low[1]]
            out_low = [out_low[1], out_low[1]]
        end
        itp_411L = LinearInterpolation(a_list_low, log10.(out_low), extrapolation_bc=Interpolations.Line())
        abndryU_411 = a_list_high[1]
        abndryL_411 = a_list_low[end]
        function itp_411(aspin)
            if aspin .>= abndryU_411
                return 10.0 .^itp_411U(aspin)
            elseif aspin .<= abndryL_411
                return -10.0 .^itp_411L(aspin)
            else
                return 0.0 ## Need buffer region for numerical stability!
            end
        end
        SR_rates[idx_fill] = itp_411(aBH)
        idx_fill += 1
    
        
        n = 4; l = 2; m = 2;
        amin_guess_422, run_high, a_list_high, out_high, run_low, a_list_low, out_low = pre_computed_sr_rates(n, l, m, alph, M_BH; n_high=100, n_low=100, delt_a=delt_a)
        if alph < 1.1
           
            if !run_high
                a_list_high = [amin_guess_422, amin_guess_422 * 1.01]
                out_high = [1e-100, 1e-100]
            elseif length(a_list_high) == 1.0
                a_list_high = [alist[1], alist[1] * 1.01]
                out_high = [out_high[1], out_high[1]]
            end
        else
            amin_guess_422 = 1.0
            a_list_high = [amin_guess_422, amin_guess_422 * 1.01]
            out_high = [1e-100, 1e-100]
        end
        itp_422U = LinearInterpolation(a_list_high, log10.(out_high), extrapolation_bc=Interpolations.Line())
        if !run_low
            a_list_low = [amin_guess_422 .* 0.99, amin_guess_422]
            out_low = [1e-100, 1e-100]
        elseif length(a_list_low) == 1.0
            a_list_low = [a_list_low[1] .* 0.99, a_list_low[1]]
            out_low = [out_low[1], out_low[1]]
        end
        itp_422L = LinearInterpolation(a_list_low, log10.(out_low), extrapolation_bc=Interpolations.Line())
        abndryU_422 = a_list_high[1]
        abndryL_422 = a_list_low[end]
        function itp_422(aspin)
            if aspin .>= abndryU_422
                return 10 .^itp_422U(aspin)
            elseif aspin .<= abndryL_422
                return -10 .^itp_422L(aspin)
            else
                return 0.0 ## Need buffer region for numerical stability!
            end
        end
        SR_rates[idx_fill] = itp_422(aBH)
        idx_fill += 1
        
        n = 4; l = 3; m = 3;
        amin_guess_433, run_high, a_list_high, out_high, run_low, a_list_low, out_low = pre_computed_sr_rates(n, l, m, alph, M_BH; n_high=100, n_low=100, delt_a=delt_a)
        if alph < 1.7
           
            if !run_high
                a_list_high = [amin_guess_433, amin_guess_433 * 1.01]
                out_high = [1e-100, 1e-100]
            elseif length(a_list_high) == 1.0
                a_list_high = [alist[1], alist[1] * 1.01]
                out_high = [out_high[1], out_high[1]]
            end
        else
            amin_guess_433 = 1.0
            a_list_high = [amin_guess_433, amin_guess_433 * 1.01]
            out_high = [1e-100, 1e-100]
        end
        itp_433U = LinearInterpolation(a_list_high, log10.(out_high), extrapolation_bc=Interpolations.Line())
        if !run_low
            a_list_low = [amin_guess_433 .* 0.99, amin_guess_433]
            out_low = [1e-100, 1e-100]
        elseif length(a_list_low) == 1.0
            a_list_low = [a_list_low[1] .* 0.99, a_list_low[1]]
            out_low = [out_low[1], out_low[1]]
        end
        itp_433L = LinearInterpolation(a_list_low, log10.(out_low), extrapolation_bc=Interpolations.Line())
        abndryU_433 = a_list_high[1]
        abndryL_433 = a_list_low[end]
        function itp_433(aspin)
            if aspin .>= abndryU_433
                return 10 .^itp_433U(aspin)
            elseif aspin .<= abndryL_433
                return -10 .^itp_433L(aspin)
            else
                return 0.0 ## Need buffer region for numerical stability!
            end
        end
        SR_rates[idx_fill] = itp_433(aBH)
        idx_fill += 1
        
        if solve_n5
            n = 5; l = 1; m = 1;
            amin_guess_511, run_high, a_list_high, out_high, run_low, a_list_low, out_low = pre_computed_sr_rates(n, l, m, alph, M_BH; n_high=100, n_low=100, delt_a=delt_a)
            if alph < 0.5
               
                if !run_high
                    a_list_high = [amin_guess_511, amin_guess_511 * 1.01]
                    out_high = [1e-100, 1e-100]
                elseif length(a_list_high) == 1.0
                    a_list_high = [alist[1], alist[1] * 1.01]
                    out_high = [out_high[1], out_high[1]]
                end
            else
                amin_guess_511 = 1.0
                a_list_high = [amin_guess_511, amin_guess_511 * 1.01]
                out_high = [1e-100, 1e-100]
            end
            itp_511U = LinearInterpolation(a_list_high, log10.(out_high), extrapolation_bc=Interpolations.Line())
            if !run_low
                a_list_low = [amin_guess_511 .* 0.99, amin_guess_511]
                out_low = [1e-100, 1e-100]
            elseif length(a_list_low) == 1.0
                a_list_low = [a_list_low[1] .* 0.99, a_list_low[1]]
                out_low = [out_low[1], out_low[1]]
            end
            itp_511L = LinearInterpolation(a_list_low, log10.(out_low), extrapolation_bc=Interpolations.Line())
            abndryU_511 = a_list_high[1]
            abndryL_511 = a_list_low[end]
            function itp_511(aspin)
                if aspin .>= abndryU_511
                    return 10.0 .^itp_511U(aspin)
                elseif aspin .<= abndryL_511
                    return -10.0 .^itp_511L(aspin)
                else
                    return 0.0 ## Need buffer region for numerical stability!
                end
            end
            SR_rates[idx_fill] = itp_511(aBH)
            idx_fill += 1

            n = 5; l = 2; m = 2;
            amin_guess_522, run_high, a_list_high, out_high, run_low, a_list_low, out_low = pre_computed_sr_rates(n, l, m, alph, M_BH; n_high=100, n_low=100, delt_a=delt_a)
            if alph < 1.1
            
                if !run_high
                    a_list_high = [amin_guess_522, amin_guess_522 * 1.01]
                    out_high = [1e-100, 1e-100]
                elseif length(a_list_high) == 1.0
                    a_list_high = [alist[1], alist[1] * 1.01]
                    out_high = [out_high[1], out_high[1]]
                end
            else
                amin_guess_522 = 1.0
                a_list_high = [amin_guess_522, amin_guess_522 * 1.01]
                out_high = [1e-100, 1e-100]
            end
            itp_522U = LinearInterpolation(a_list_high, log10.(out_high), extrapolation_bc=Interpolations.Line())
            if !run_low
                a_list_low = [amin_guess_522 .* 0.99, amin_guess_522]
                out_low = [1e-100, 1e-100]
            elseif length(a_list_low) == 1.0
                a_list_low = [a_list_low[1] .* 0.99, a_list_low[1]]
                out_low = [out_low[1], out_low[1]]
            end
            itp_522L = LinearInterpolation(a_list_low, log10.(out_low), extrapolation_bc=Interpolations.Line())
            abndryU_522 = a_list_high[1]
            abndryL_522 = a_list_low[end]
            function itp_522(aspin)
                if aspin .>= abndryU_522
                    return 10 .^itp_522U(aspin)
                elseif aspin .<= abndryL_522
                    return -10 .^itp_522L(aspin)
                else
                    return 0.0 ## Need buffer region for numerical stability!
                end
            end
            SR_rates[idx_fill] = itp_522(aBH)
            idx_fill += 1
        
            n = 5; l = 3; m = 3;
            amin_guess_533, run_high, a_list_high, out_high, run_low, a_list_low, out_low = pre_computed_sr_rates(n, l, m, alph, M_BH; n_high=100, n_low=100, delt_a=delt_a)
            if alph < 1.7
            
                if !run_high
                    a_list_high = [amin_guess_533, amin_guess_533 * 1.01]
                    out_high = [1e-100, 1e-100]
                elseif length(a_list_high) == 1.0
                    a_list_high = [alist[1], alist[1] * 1.01]
                    out_high = [out_high[1], out_high[1]]
                end
            else
                amin_guess_533 = 1.0
                a_list_high = [amin_guess_533, amin_guess_533 * 1.01]
                out_high = [1e-100, 1e-100]
            end
            itp_533U = LinearInterpolation(a_list_high, log10.(out_high), extrapolation_bc=Interpolations.Line())
            if !run_low
                a_list_low = [amin_guess_533 .* 0.99, amin_guess_533]
                out_low = [1e-100, 1e-100]
            elseif length(a_list_low) == 1.0
                a_list_low = [a_list_low[1] .* 0.99, a_list_low[1]]
                out_low = [out_low[1], out_low[1]]
            end
            itp_533L = LinearInterpolation(a_list_low, log10.(out_low), extrapolation_bc=Interpolations.Line())
            abndryU_533 = a_list_high[1]
            abndryL_533 = a_list_low[end]
            function itp_533(aspin)
                if aspin .>= abndryU_533
                    return 10 .^itp_533U(aspin)
                elseif aspin .<= abndryL_533
                    return -10 .^itp_533L(aspin)
                else
                    return 0.0 ## Need buffer region for numerical stability!
                end
            end
            SR_rates[idx_fill] = itp_533(aBH)
            idx_fill += 1
            
            n = 5; l = 4; m = 4
            amin_guess_544, run_high, a_list_high, out_high, run_low, a_list_low, out_low = pre_computed_sr_rates(n, l, m, alph, M_BH; n_high=100, n_low=100, delt_a=delt_a)
            if alph < 2.2
            
                if !run_high
                    a_list_high = [amin_guess_544, amin_guess_544 * 1.01]
                    out_high = [1e-100, 1e-100]
                elseif length(a_list_high) == 1.0
                    a_list_high = [alist[1], alist[1] * 1.01]
                    out_high = [out_high[1], out_high[1]]
                end
            else
                amin_guess_544 = 1.0
                a_list_high = [amin_guess_544, amin_guess_544 * 1.01]
                out_high = [1e-100, 1e-100]
            end
            itp_544U = LinearInterpolation(a_list_high, log10.(out_high), extrapolation_bc=Interpolations.Line())
            if !run_low
                a_list_low = [amin_guess_544 .* 0.99, amin_guess_544]
                out_low = [1e-100, 1e-100]
            elseif length(a_list_low) == 1.0
                a_list_low = [a_list_low[1] .* 0.99, a_list_low[1]]
                out_low = [out_low[1], out_low[1]]
            end
            itp_544L = LinearInterpolation(a_list_low, log10.(out_low), extrapolation_bc=Interpolations.Line())
            abndryU_544 = a_list_high[1]
            abndryL_544 = a_list_low[end]
            function itp_544(aspin)
                if aspin .>= abndryU_544
                    return 10 .^itp_544U(aspin)
                elseif aspin .<= abndryL_544
                    return -10 .^itp_544L(aspin)
                else
                    return 0.0 ## Need buffer region for numerical stability!
                end
            end
            SR_rates[idx_fill] = itp_544(aBH)
            
        end
    end
    
    if debug
        print("SR rates @ prod \t", SR_rates, "\n")
    end
    
    rates = load_rate_coeffs(mu, M_BH, aBH, fa; non_rel=non_rel, input_data=input_data, solve_n4=solve_n4, solve_n5=solve_n5, amin_211=amin_guess_211, amin_322=amin_guess_322, amin_433=amin_guess_433)
    
    if debug
        println("Rates loaded...")
        println(rates)
    end
    
    ### for testing
    trace_terms_tot = 0
    if trace_term
        println("SR \t", SR_rates[trace_idx])
        rate_keys = collect(keys(rates))
        for i in 1:length(rate_keys)
            idxV, sgn = key_to_indx(rate_keys[i]; solve_n4=true, solve_n5=true)
            toV = sum(idxV .== trace_idx)
            # testGW = occursin("GW", rate_keys[i])
            testGW = false
            if (toV >= 1)&&(!testGW)
               print(rate_keys[i], "\t")
               trace_terms_tot += 1
            end
        end
    end
    ###
    
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
        
        
        SR_rates[1] = itp_211(u_real[spinI])
        SR_rates[2] = itp_311(u_real[spinI])
        SR_rates[3] = itp_322(u_real[spinI])

        
        if solve_n4
            SR_rates[4] = itp_411(u_real[spinI])
            SR_rates[5] = itp_422(u_real[spinI])
            SR_rates[6] = itp_433(u_real[spinI])
            
            
            if solve_n5
                SR_rates[7] = itp_511(u_real[spinI])
                SR_rates[8] = itp_522(u_real[spinI])
                SR_rates[9] = itp_533(u_real[spinI])
                SR_rates[10] = itp_544(u_real[spinI])

            end
        end
            
        
        if (u_real[1] .> Emax2)&&(SR_rates[1] > 0)
            SR_rates[1] *= 0.0
        end
        
    
        # SR terms
        du[spinI] = 0.0
        du[massI] = 0.0
        for i in 1:idx_lvl
            du[i] = SR_rates[i] .* u_real[i] ./ mu
            if (!max_m_2)||((i <= 5)||(i == 7)||(i == 8))
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
        
        if max_m_2
            du[massI] *= 0.0
        end
        
        
        
        return
    end
    
    state_check = 3000
    osc_thresh = 2
    #### CALLBACK 1
    function check_oscillating(u, t, integrator)
        
        if length(integrator.sol.t) <= state_check
            return false
        else
            for i in 1:idx_lvl
                
                state = [integrator.sol.u[j][i] for j in 1:length(integrator.sol.u)]
                short = round.(state[end-state_check:end], digits=1)
                cnt_extrm = 0
                cnt_extrm += sum(abs.(diff(sign.(short[2:end] .- short[1:end-1])))) / 2
                
                maxV = maximum(short)
                minV = minimum(short)
                medianV = median(short)
                
                cnd1 = (maxV - medianV) > osc_thresh
                cnd2 = (medianV - minV) > osc_thresh
                if (cnt_extrm > 10)&&cnd1&&cnd2
                    return true
                end
                
                
                
            end
            
        end
        
        return false
    end
    
    function affect_oscillating!(integrator)
        for i in 1:idx_lvl
            state = [integrator.sol.u[j][i] for j in 1:length(integrator.sol.u)]
            
            short = round.(state[end-state_check:end], digits=1)
            cnt_extrm = 0
            maxV = maximum(short)
            cnt_extrm += sum(abs.(diff(sign.(short[2:end] .- short[1:end-1])))) / 2
            maxV = maximum(short)
            minV = minimum(short)
            medianV = median(short)
            
            cnd1 = (maxV - medianV) > osc_thresh
            cnd2 = (medianV - minV) > osc_thresh
            
            if (cnt_extrm > 10)&&cnd2&&cnd2&&(turn_off[i] == false)
                turn_off[i] = true
#                 integrator.u[i] = log.(e_init)
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
    
    
    function check_stationary(u, t, integrator)
         if length(integrator.sol.t) < 3000
            return false
        end
        cnt_station = 0
        
        for i in 1:idx_lvl
            state = [integrator.sol.u[j][i] for j in 1:length(integrator.sol.u)]
            short = round.(state[end-1000:end], digits=1)
            
            cngs = sum(diff(short))
            if (cngs == 0)
                cnt_station += 1
            end
        end
        
        state = [integrator.sol.u[j][spinI] for j in 1:length(integrator.sol.u)]
        short = round.(state[end-1000:end], digits=2)
        cngs = sum(diff(short))
        
        state = [integrator.sol.u[j][massI] for j in 1:length(integrator.sol.u)]
        short = round.(state[end-1000:end], digits=2)
        cngs2 = sum(diff(short))
        
        check_SR = 0
        ## check if superrad is turning on
        for i in 1:idx_lvl
            t_SR = (abs.(SR_rates[i]) ./ hbar .* 3.15e7).^(-1)
            if (t > t_SR)&&(integrator.sol.t[end-2] < t_SR)
                for i in 1:idx_lvl
                    turn_off[i] = false
                end
                check_SR = 1
            end
        end
            
        if (cnt_station == idx_lvl)
            if (cngs == 0)&&(cngs2 == 0)&&(check_SR == 0)
                return true
            else
                for i in 1:idx_lvl
                    turn_off[i] = false
                end
                return false
            end
        else
            return false
        end
        
    end
    
    function affect_stationary!(integrator)
        for i in 1:idx_lvl
            turn_off[i] = true
        end
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
        
        for i in 1:idx_lvl
            if (u[i] .< log.(1e-70))&&(SR_rates[i] < 0)
                turn_off[i] = true
            end
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
                # reltol[i] = default_reltol
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
            
            if (u[i] > log.(e_init)) && condBN && (sign_flip[i]==false) && (du[i] != 0.0)
                append!(tlist, abs.(1.0 ./ du[i]))
            end
        end
        
        append!(tlist, def_spin_tol ./ du[spinI])
        tmin = minimum(abs.(tlist))
       
        # print("Tmin \t", integrator.dt, "\t", tlist, "\t", integrator.opts.reltol, "\t", integrator.opts.abstol, "\n")
        
        # print(du, "\n")
        if (integrator.dt ./ tmin .>= 0.1)
            return true
        elseif (integrator.dt ./ tmin .<= 0.001)
            return true
        elseif (integrator.dt .<= 1e-10)
            return true
        # elseif (integrator.dt ./ integrator.t .<= 1e-4)
        #     integrator.opts.dtmin = integrator.t * 1e-4
        #    return true
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
    
    # Callback 3
    function check_spin(u, t, integrator)
        wait += 1
        u_real = exp.(u)
        if debug && (wait%1==0)
            if solve_n4
                if !solve_n5
                    @printf("Time and Vals: \t %.7f  %.3e  %.3e  %.3e  %.3e  %.3e  %.3e  %.3e %.3e \n", t, u_real[1], u_real[2], u_real[3], u_real[4], u_real[5], u_real[6], u_real[7], u_real[8])
                    # println(integrator.opts.abstol, "\t", integrator.opts.reltol )
                    # print(integrator.dt ./ t, "\n")
#                    if (t > 4200)&&(integrator.opts.reltol[1] < 0.1)
#                        integrator.opts.reltol[1] *= 2.0
#                        print(integrator.opts.reltol, "\n")
#                        print(integrator.opts.abstol, "\n")
#                    end
                else
                    @printf("Time and Vals: \t %.7f  %.3e  %.3e  %.3e  %.3e  %.3e  %.3e  %.3e  %.3e  %.3e  %.3e %.3e  %.3e \n", t, u_real[1], u_real[2], u_real[3], u_real[4], u_real[5], u_real[6], u_real[7], u_real[8], u_real[9], u_real[10], u_real[11], u_real[12])
                    # print(integrator.opts.abstol, "\t ",integrator.opts.reltol, "\n")
                    # println(SR_rates)
                end
                # du = get_du(integrator)
                # print(integrator.dt, "\t", du, "\t", "\n\n")
            else
                print(t, "\t", u_real[1], "\t", u_real[2], "\t", u_real[3], "\t", u_real[4], "\t", u_real[5], "\n")
                # println(integrator.opts.abstol, "\t", integrator.opts.reltol )
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
    
    
    
    output_trace = []
    function sve_tr(u, t, integrator)
        out_row = zeros(trace_terms_tot + 1)
        out_row[1] = t
        cnt = 2
        for i in 1:length(rate_keys)
            idxV, sgn = key_to_indx(rate_keys[i]; solve_n4=solve_n4, solve_n5=solve_n5)
            # testGW = occursin("GW", rate_keys[i])
            testGW = false
            if (sum(idxV .== trace_idx) .> 0)&&!testGW
                u_term_tot = 1.0
                for j in 1:length(sgn)
                    
                    if (idxV[j] .<= idx_lvl)&&(idxV[j] .> 0)
                        u_term_tot *= exp.(u[idxV[j]])
                    end
                   
                end
                
                if (trace_idx == idxV[1])&&(idxV[1] == idxV[2])
                    out_row[cnt] += 2 * rates[rate_keys[i]] * u_term_tot
                else
                    out_row[cnt] += rates[rate_keys[i]] * u_term_tot
                end
                cnt += 1
            end
        end
        
        if length(output_trace) == 0
            output_trace = out_row'
        else
            output_trace = cat(output_trace, out_row', dims=1)
        end
        return false
    end
    function affect_tr!(integrator)
        nothing;
    end
    
    
    rP = 1.0 .+ sqrt.(1 - aBH.^2)
    
    tlist = []
    for i in 1:idx_lvl
        if SR_rates[i] .> 0
            append!(tlist, (abs.(SR_rates[i]) ./ hbar .* 3.15e7).^(-1))
        else
            append!(tlist, 1e6)
        end
    end
    dt_guess = minimum(tlist) ./ 5.0

    if debug
        print("Time guess \t", dt_guess, "\n")
    end
    
    
    max_real_time = 20.0 # Set the maximum allowed time for integration (in minutes)
    max_real_time *= 60 # convert to seconds
    start_time = Dates.now() # Initialize the start time
    
    # Define the callback function
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
    cback_station = DiscreteCallback(check_stationary, affect_stationary!, save_positions=(false, true))
    cback_osc = DiscreteCallback(check_oscillating, affect_oscillating!, save_positions=(false, true))
    cbackdt = DiscreteCallback(check_timescale, affect_timescale!, save_positions=(false, true))
    cbackspin = DiscreteCallback(check_spin, affect_spin!, save_positions=(false, true))
    if trace_term
        sve_trce = DiscreteCallback(sve_tr, affect_tr!, save_positions=(false, true))
    end
    
    
    if solve_n4
        if solve_n5
            # cbset = CallbackSet(cbackspin)
            # cbset = CallbackSet(cbackspin, cbackdt)
            cbset = CallbackSet(cbackspin, cbackdt, cback_station, cback_osc, callbackTIME)
            if trace_term
                cbset = CallbackSet(cbackspin, cbackdt, cback_station, cback_osc, sve_trce)
            end
        else
            cbset = CallbackSet(cbackspin, cbackdt, cback_station, cback_osc, callbackTIME)
            # cbset = CallbackSet(cbackspin, cbackdt)
        end
    else
        cbset = CallbackSet(cbackspin, cbackdt, cback_station, callbackTIME)
    end
    # cbset = CallbackSet(cbackspin, cbackdt, cback_station)
    # cbset = CallbackSet(cbackspin, cbackdt, cback_osc)
    
    if solve_n4
        prob = ODEProblem(RHS_ax!, y0, tspan, Mvars, reltol=reltol, abstol=1e-10)
        # prob = ODEProblem(RHS_ax!, y0, tspan, Mvars, reltol=1e-5, abstol=1e-10)
        # sol = solve(prob, Vern6(), dt=dt_guess, saveat=saveat, callback=cbset, maxiters=5e6, dtmin=(dt_guess / 1e5), force_dtmin=true)
        if high_p
            # cbset = CallbackSet(cbackspin, cbackdt)
            # cbset = CallbackSet(cbackspin, cbackdt, cback_station, callbackTIME)
            cbset = CallbackSet(cbackspin, cbackdt, callbackTIME)
        end
        
        # sol = solve(prob, Rosenbrock23(autodiff=false), dt=dt_guess, saveat=saveat, callback=cbset, maxiters=5e6, dtmin=(dt_guess / 1e5), force_dtmin=true)
        sol = solve(prob, TRBDF2(autodiff=false), dt=dt_guess, saveat=saveat, callback=cbset, maxiters=5e6)
        # sol = solve(prob, Vern6(), dt=dt_guess, saveat=saveat, callback=cbset, maxiters=5e6)
        # sol = solve(prob, QNDF(autodiff=false), dt=dt_guess, saveat=saveat, callback=cbset, maxiters=5e6, dtmin=(dt_guess / 1e5), force_dtmin=true)
        # sol = solve(prob, Euler(), dt=dt_guess / 1e2, saveat=saveat, callback=cbset)
    else
        prob = ODEProblem(RHS_ax!, y0, tspan, Mvars,  reltol=reltol, abstol=1e-10)
        sol = solve(prob, Rosenbrock23(autodiff=false), dt=dt_guess, saveat=saveat, callback=cbset, maxiters=5e6, dtmin=(dt_guess / 1e5),  force_dtmin=true)
    end
    
    
    
    state211 = [exp.(sol.u[i][1]) for i in 1:length(sol.u)]
    state311 = [exp.(sol.u[i][2]) for i in 1:length(sol.u)]
    state322 = [exp.(sol.u[i][3]) for i in 1:length(sol.u)]
    if solve_n4
        state411 = [exp.(sol.u[i][4]) for i in 1:length(sol.u)]
        state422 = [exp.(sol.u[i][5]) for i in 1:length(sol.u)]
        state433 = [exp.(sol.u[i][6]) for i in 1:length(sol.u)]
        if solve_n5
            state511 = [exp.(sol.u[i][7]) for i in 1:length(sol.u)]
            state522 = [exp.(sol.u[i][8]) for i in 1:length(sol.u)]
            state533 = [exp.(sol.u[i][9]) for i in 1:length(sol.u)]
            state544 = [exp.(sol.u[i][10]) for i in 1:length(sol.u)]
        end
    end
    spinBH = [exp.(sol.u[i][spinI]) for i in 1:length(sol.u)]
    MassB = [exp.(sol.u[i][massI]) for i in 1:length(sol.u)]

    if trace_term
        writedlm("test_store/Trace_term_$(trace_idx)_.dat", output_trace)
    end

    if debug
        if !solve_n4
            print("Initial \t", state211[1], "\t", state311[1], "\t", state322[1],  "\t", spinBH[1], "\t", MassB[1], "\n")
            print("Final \t", state211[end], "\t", state311[end], "\t", state322[end], "\t", spinBH[end], "\t", MassB[end], "\n")
        else
            if !solve_n5
                print("Initial \t", state211[1], "\t", state311[1] , "\t", state322[1] , "\t",  state411[1] , "\t", state422[1] , "\t", state433[1] , "\t", spinBH[1], "\t", MassB[1], "\n")
                print("Final \t", state211[end] , "\t", state311[end], "\t", state322[end],"\t", state411[end],"\t", state422[end] , "\t",  state433[end] , "\t", spinBH[end], "\t", MassB[end], "\n")
            else
                print("Initial \t", state322[1], "\t", state422[1] , "\t", state433[1], "\t", state522[1], "\t", state533[1] , "\t", state544[1], "\t", spinBH[1], "\t", MassB[1], "\n")
                print("Final \t", state322[end], "\t", state422[end] , "\t", state433[end], "\t", state522[end], "\t", state533[end] , "\t", state544[end], "\t", spinBH[end], "\t", MassB[end], "\n")
            end
        end
    end

    if return_all_info
        if !solve_n4
            return sol.t, state211, state311, state322, spinBH, MassB
        else
            if !solve_n5
                return sol.t, state211, state311, state322, state411, state422, state433, spinBH, MassB
            else
                return sol.t, state211, state311, state322, state411, state422, state433, state511, state522, state533, state544, spinBH, MassB
            end
        end
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
