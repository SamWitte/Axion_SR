include("solve_sr_rates.jl")
using NPZ
using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
    
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

run_leaver = parsed_args["run_leaver"];
solve_for_zeros = parsed_args["solve_for_zeros"];
solve_gridded = parsed_args["solve_gridded"];
nlmIn = parsed_args["nlmIn"]

function main_gg(run_leaver, solve_for_zeros, solve_gridded)

    Ntot_safe = 5000
    debug = true
    xtol=1e-4
    ftol=1e-50
    iter=50

    Npoints = 60
    Iter = 20
    cvg_acc = 1e-3
    Npts_r = 1000
    prec=200
    
    debug=true


    aPts = 40
    alpha_pts = 60
    a_max = 0.998
    alist = LinRange(0.1, a_max, aPts)


    nmax = 8
    mmax = 5
    loop_list = []
    for n in 1:nmax, l in 1:(n - 1), m in 1:l
        if m > mmax
            continue
        end
        if run_leaver
            ft1 = "Imag_zero"
            ft2 = "Imag_erg_pos"
        else
            ft1 = "Imag_zeroC"
            ft2 = "Imag_ergC_pos"
        end
        
        for (flag, suffix) in ((solve_for_zeros, ft1), (solve_gridded, ft2))
            if flag
                file_out = "$(suffix)_$(n)$(l)$(m).$(suffix == ft1 ? "dat" : "npz")"
                full_path = "rate_sve/" * file_out
                
                if !isfile(full_path)
                    if nlmIn == "000"
                        push!(loop_list, [n, l, m])
                    else
                        if (suffix == ft2)&&(nlmIn == string(n)*string(l)*string(m))
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
            if run_leaver
                file_out = "Imag_zero_$(n)$(l)$(m).dat"
            else
                file_out = "Imag_zeroC_$(n)$(l)$(m).dat"
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
        
                if run_leaver
                    testF = find_im_part(10 .^ alphList[i] ./ (GNew .* M), M, 0.998, n, l, m; Ntot_force=Ntot_safe, return_both=false, for_s_rates=true)
                else
                    wR, testF = eigensys_Cheby(M, 0.998, 10 .^ alphList[i] ./ (GNew .* M), n, l, m, debug=false, return_wf=false, Npoints=Npoints, Iter=Iter, cvg_acc=cvg_acc, prec=prec)
                end
                if testF .< 0
                    continue
                end

                while_found = false
                cnt = 0
                while !while_found
                    
                    if run_leaver
                        testF = find_im_part(10 .^ alphList[i] ./ (GNew .* M), M, a_guess, n, l, m; Ntot_force=Ntot_safe, return_both=false, for_s_rates=true)
                    else
                        wR, testF = eigensys_Cheby(M, a_guess, 10 .^ alphList[i] ./ (GNew .* M), n, l, m, debug=false, return_wf=false, Npoints=Npoints, Iter=Iter, cvg_acc=cvg_acc, prec=prec)
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
                #println(10 .^ alphList[i], "\t", a_guess)
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
            alphList = LinRange(log10.(0.03), log10.(alpha_max), alpha_pts)
            
            if run_leaver
                file_out = "Imag_erg_pos_$(n)$(l)$(m).npz"
                file_out2 = "Imag_erg_neg_$(n)$(l)$(m).npz"
            else
                file_out = "Imag_ergC_pos_$(n)$(l)$(m).npz"
                file_out2 = "Imag_ergC_neg_$(n)$(l)$(m).npz"
            end
           
            M=1;
            store_P = zeros(alpha_pts, aPts, 3)
            store_M = zeros(alpha_pts, aPts, 3)

            # get sign flip a
            if run_leaver
                zerolist= readdlm("rate_sve/Imag_zero_$(n)$(l)$(m).dat")
            else
                zerolist= readdlm("rate_sve/Imag_zeroC_$(n)$(l)$(m).dat")
            end
            itp = LinearInterpolation(zerolist[:, 1], zerolist[:, 2], extrapolation_bc=Line())

            
            alistP = LinRange(a_min, a_max, aPts)
            for i in 1:alpha_pts
                a_mid = itp(10 .^ alphList[i])
                for j in 1:aPts
                    if run_leaver
                        e_imgP = find_im_part(10 .^ alphList[i] ./ (GNew .* M), M, alistP[j], n, l, m; Ntot_force=Ntot_safe, return_both=false, for_s_rates=true)
                    else
                        wR, e_imgP = eigensys_Cheby(M, alistP[j], 10 .^ alphList[i] ./ (GNew .* M), n, l, m, debug=false, return_wf=false, Npoints=Npoints, Iter=Iter, cvg_acc=cvg_acc, prec=prec)
                    end
                    
                    if e_imgP > 0
                        store_P[i, j, :] = [10 .^ alphList[i] alistP[j] e_imgP[1]]
                        store_M[i, j, :] = [10 .^ alphList[i] alistP[j] 1e-100]
                    else
                        store_M[i, j, :] = [10 .^ alphList[i] alistP[j] -e_imgP[1]]
                        store_P[i, j, :] = [10 .^ alphList[i] alistP[j] 1e-100]
                    end
                
                end
            end
            
            npzwrite("rate_sve/"*file_out, store_P)
            npzwrite("rate_sve/"*file_out2, store_M)
        end

    end

end

main_gg(run_leaver, solve_for_zeros, solve_gridded)
