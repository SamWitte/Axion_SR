include("solve_sr_rates.jl")
include("Core/constants.jl")
include("state_utils.jl")
using NPZ
using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
    

        "--run_analytic"
            arg_type = Bool
            default = false
            
        "--run_leaver"
            arg_type = Bool
            default = false

        "--solve_for_zeros"
            arg_type = Bool
            default = false
            
        "--solve_gridded"
            arg_type = Bool
            default = false
            
        "--nlmIn"
            arg_type = String
            default = "000"
            
    end
    return parse_args(s)
end


parsed_args = parse_commandline()

run_analytic = parsed_args["run_analytic"];
run_leaver = parsed_args["run_leaver"];
solve_for_zeros = parsed_args["solve_for_zeros"];
solve_gridded = parsed_args["solve_gridded"];
nlmIn = parsed_args["nlmIn"]

function main_gg(run_leaver, run_analytic, solve_for_zeros, solve_gridded)

    Ntot_safe = 5000
    debug = true
    xtol=1e-4
    ftol=1e-50
    iter=60

    Npoints = 110
    Iter = 30
    cvg_acc = 1e-10
    der_acc=1e-20
    Npts_r = 1000
    prec=200
    
    debug=true

    Lcheb=8
    aPts = 12
    alpha_pts = 10
    a_max = 0.998
    atemp = 10 .^LinRange(-3, log10.(0.99), aPts)
    


    nmax = 8
    mmax = 7
    loop_list = []
    for n in 1:nmax, l in 1:(n - 1), m in 1:l
        if m > mmax
            continue
        end
        if run_analytic
            ft1 = "Imag_zeroC"
            ft2 = "Imag_ergC_pos"
        else
            if run_leaver
                ft1 = "Imag_zero"
                ft2 = "Imag_erg_pos"
            else
                ft1 = "Imag_zeroC"
                ft2 = "Imag_ergC_pos"
            end
        end
        
        for (flag, suffix) in ((solve_for_zeros, ft1), (solve_gridded, ft2))
            if flag
                state_str = format_state_for_filename(n, l, m)
                file_out = "$(suffix)_$(state_str).$(suffix == ft1 ? "dat" : "npz")"
                full_path = "rate_sve/" * file_out

                if !isfile(full_path)
                    if nlmIn == "000"
                        push!(loop_list, [n, l, m])
                    else
                        # Parse nlmIn and compare quantum numbers
                        n_in, l_in, m_in = parse_state_string(nlmIn)
                        if (n == n_in) && (l == l_in) && (m == m_in)
                            push!(loop_list, [n, l, m])
                        end
                    end
                    break  # no need to check other flags once added
                end
            end
        end
    end


    if solve_for_zeros
        for nlm in loop_list
            println(nlm)
            n = nlm[1]; l = nlm[2]; m = nlm[3];
            
            alpha_max = a_max .* m ./ (2 .* (1 .+ sqrt.(1 .- a_max.^2))) .* 1.3
            alphList = LinRange(log10.(0.03), log10.(alpha_max), alpha_pts)
            state_str = format_state_for_filename(n, l, m)
            if run_leaver
                file_out = "Imag_zero_$(state_str).dat"
            else
                file_out = "Imag_zeroC_$(state_str).dat"
            end
           
            M=1;
            store_out = zeros(alpha_pts, 2)

            for i in 1:alpha_pts
                if debug
                    println("\n", "Alpha: ", 10 .^alphList[i])
                end
                a_min = 0.01
                a_guess = 0.5
                a_max_loop = a_max
        
                testF = find_im_part(10 .^ alphList[i] ./ (GNew .* M), M, 0.998, n, l, m; Ntot_force=Ntot_safe, return_both=false, for_s_rates=true)
              
                if testF .< 0
		    wR, testF = eigensys_Cheby(M, 0.998, 10 .^ alphList[i] ./ (GNew .* M), n, l, m, debug=false, return_wf=false, Npoints=Npoints, Iter=Iter, cvg_acc=cvg_acc, prec=prec, sfty_run=false, L=Lcheb)
		    if testF .< 0
 		        println("Stop running \t", testF)
                        continue
		    end
                end

                while_found = false
                cnt = 0
                while !while_found
                    if run_analytic
                        testF = sr_rates(n, l, m, 10 .^ alphList[i] ./ (GNew .* M), M, a_guess) .* GNew .* M
                    else
                        if run_leaver
                            testF = find_im_part(10 .^ alphList[i] ./ (GNew .* M), M, a_guess, n, l, m; Ntot_force=Ntot_safe, return_both=false, for_s_rates=true)
                        else
                            if a_guess > 0.95
                                npts_use = 120
                            else
                                npts_use = Npoints
                            end
                            wR, testF = eigensys_Cheby(M, a_guess, 10 .^ alphList[i] ./ (GNew .* M), n, l, m, debug=false, return_wf=false, Npoints=Npoints, Iter=Iter, cvg_acc=cvg_acc, prec=prec, sfty_run=false, L=Lcheb, der_acc=der_acc)
                        end
                    end
                    # print(cnt, "\t", a_guess, "\t", testF, "\n")
                    if testF > 0
                        a_max_loop = a_guess
                        a_guess = 0.5 .* (a_guess .+ a_min)
                    else
                        a_min = a_guess
                        a_guess = 0.5 .* (a_guess .+ a_max_loop)
                    end
                    
                    if (a_max_loop .- a_min) < 1e-4
                        while_found = true
                    end

                    if cnt > 30
                        while_found = true
                    end
                    cnt += 1
                end
                println(10 .^ alphList[i], "\t", a_guess)
                store_out[i, :] = [10 .^ alphList[i] a_guess]
            end
            
            store_out = store_out[store_out[:,1] .!== 0.0, :]
            writedlm("rate_sve/"*file_out, store_out)
            
        end
    end

    if solve_gridded
        a_min = 0.01
        for nlm in loop_list
            println(nlm)
            n = nlm[1]; l = nlm[2]; m = nlm[3];
           
            alpha_max = a_max .* m ./ (2 .* (1 .+ sqrt.(1 .- a_max.^2))) .* 1.3
            # alphList = LinRange(log10.(0.03), log10.(alpha_max), alpha_pts)
	    alphList = LinRange(log10.(0.3), log10.(alpha_max), alpha_pts)

            state_str = format_state_for_filename(n, l, m)
            if run_leaver
                file_out = "Imag_erg_pos_$(state_str).npz"
                file_out2 = "Imag_erg_neg_$(state_str).npz"
            else
                file_out = "Imag_ergC_pos_$(state_str).npz"
                file_out2 = "Imag_ergC_neg_$(state_str).npz"
            end

            M=1;
            store_P = zeros(alpha_pts, aPts, 3)
            store_M = zeros(alpha_pts, aPts, 3)

            # get sign flip a
            if run_leaver
                zerolist= readdlm("rate_sve/Imag_zero_$(state_str).dat")
            else
                zerolist= readdlm("rate_sve/Imag_zeroC_$(state_str).dat")
            end
            itp = LinearInterpolation(zerolist[:, 1], zerolist[:, 2], extrapolation_bc=Line())

            
            # alistP = a_max .- atemp
            # alistP = reverse(alistP)
            alistP = LinRange(0.01, a_max, aPts)
	    a_diff = alistP[2] .- alistP[1]
            
            for i in 1:alpha_pts
                a_mid = itp(10 .^ alphList[i])
                erg_store = nothing
                erg_store_prev = nothing
                cvg_thresh = 10.0  # max relative change in Re(ω) between spin steps before retrying (failures jump by orders of magnitude)
                for j in 1:length(alistP)
                    if alistP[j] > 0.925
                        npts_use = 120
                    else
                        npts_use = Npoints
                    end
                    if run_analytic
                        e_imgP = sr_rates(n, l, m, 10 .^ alphList[i] ./ (GNew .* M), M, alistP[j]) .* GNew .* M
                    else
                        if run_leaver
                            e_imgP = find_im_part(10 .^ alphList[i] ./ (GNew .* M), M, alistP[j], n, l, m; Ntot_force=Ntot_safe, return_both=false, for_s_rates=true)
                        else
                            if j == 1

                                wR, e_imgP = eigensys_Cheby(M, alistP[j], 10 .^ alphList[i] ./ (GNew .* M), n, l, m, debug=false, return_wf=false, Npoints=npts_use, Iter=Iter, cvg_acc=cvg_acc, prec=prec, sfty_run=true, L=Lcheb)
                                erg_store = wR + im .* e_imgP
                                erg_store_prev = nothing
                            else
                                # Linear extrapolation from last two good points if available, else last point
                                if !isnothing(erg_store_prev)
                                    nu_extrap = 2 .* erg_store .- erg_store_prev
                                else
                                    nu_extrap = erg_store
                                end

                                wR, e_imgP = eigensys_Cheby(M, alistP[j], 10 .^ alphList[i] ./ (GNew .* M), n, l, m, debug=false, return_wf=false, Npoints=npts_use, Iter=Iter, cvg_acc=cvg_acc, prec=prec, sfty_run=true, L=Lcheb, nu_guess=nu_extrap)

                                # Validate: if result is non-finite or Re(ω) jumps too far, retry with analytic guess
                                if !isfinite(wR) || !isfinite(e_imgP) || abs(wR - real(erg_store)) > cvg_thresh * abs(real(erg_store))
                                    println("Continuation failed at j=$j a=$(alistP[j]), retrying with analytic guess")
                                    wR, e_imgP = eigensys_Cheby(M, alistP[j], 10 .^ alphList[i] ./ (GNew .* M), n, l, m, debug=false, return_wf=false, Npoints=npts_use, Iter=Iter, cvg_acc=cvg_acc, prec=prec, sfty_run=true, L=Lcheb)
                                end

                                # Update history only if result is valid
                                if isfinite(wR) && isfinite(e_imgP)
                                    erg_store_prev = erg_store
                                    erg_store = wR + im .* e_imgP
                                end
                            end
                        end
                    end
                    
                    if e_imgP > 0
                        store_P[i, j, :] = [10 .^ alphList[i] alistP[j] e_imgP[1]]
                        store_M[i, j, :] = [10 .^ alphList[i] alistP[j] 1e-100]
                    else
                        store_M[i, j, :] = [10 .^ alphList[i] alistP[j] -e_imgP[1]]
                        store_P[i, j, :] = [10 .^ alphList[i] alistP[j] 1e-100]
                    end
                    println("Here \t", 10 .^ alphList[i], "\t", alistP[j], "\t", e_imgP[1])
                end
                    
            end
            
            npzwrite("rate_sve/"*file_out, store_P)
            npzwrite("rate_sve/"*file_out2, store_M)
        end

    end

end

main_gg(run_leaver, run_analytic, solve_for_zeros, solve_gridded)
