include("solve_sr_rates.jl")
include("cheby_eigs.jl")

n_lvls = [2 3 4 5]
alph = LinRange(0.03, 1, 100)

l = 1
m = 1
M = 22.0
a = 0.9

Npts_Bnd=4000
rmaxT=100
Ntot_safe=5000

Npoints = 100
Iter = 20

for i in 1:length(n_lvls)
    print(n_lvls[i], "\n")
    outPr = []
    outPi = []
    
    outPr_2 = []
    outPi_2 = []
    
    for j in 1:length(alph)
        print("\t\t", alph[j], "\n\n")
    	mu = alph[j] ./ (GNew * M)
        if alph[j] < 0.1
            # wR, wI = find_im_part(mu, M, a, n_lvls[i], l, m; debug=false, iter=50, xtol=1e-20, ftol=1e-90, return_both=true, for_s_rates=true)
            
            rl, r1, erg = solve_radial(mu, M, a, n_lvls[i], l, m; rpts=Npts_Bnd, rmaxT=rmaxT, return_erg=true, Ntot_safe=Ntot_safe, iter=50, xtol=1e-20, ftol=1e-90, return_both=true, for_s_rates=true)
            wR = real.(erg)
            wI = imag.(erg)
        else
            itpR = LinearInterpolation(alph[1:j-1], outPr, extrapolation_bc=Line())
            itpI = LinearInterpolation(alph[1:j-1], outPi, extrapolation_bc=Line())
            
            guessV = itpR(alph[j]) .+ itpI(alph[j]) * im
            # print(alph[j], "\t", guessV, "\n")
            wR, wI = find_im_part(mu, M, a, n_lvls[i], l, m; debug=false, iter=50, xtol=1e-20, ftol=1e-90, return_both=false, for_s_rates=true, erg_Guess=guessV)
            
            rl, r1, erg = solve_radial(mu, M, a, n_lvls[i], l, m; rpts=Npts_Bnd, rmaxT=rmaxT, return_erg=true, Ntot_safe=Ntot_safe, iter=50, xtol=1e-20, ftol=1e-90, return_both=true, for_s_rates=true)
            wR = real.(erg)
            wI = imag.(erg)
        end
        append!(outPr, wR)
        append!(outPi, wI)
        println("First go \t ", wR, "\t", wI)
        
        
        
        wR, wI, x_values, y_values = eigensys_Cheby(M, a, alph[j] ./ (GNew .* M), n_lvls[i], l, m, debug=false, return_wf=true, Npoints=Npoints, Iter=Iter)
        println("Cheb go \t ", wR, "\t", wI)
        append!(outPr_2, wR)
        append!(outPi_2, wI)
        
        
        if j == 1
            writedlm("test_store/RadTest_LowA_1_n_$(n_lvls[i])_l_$(l)_m_$(m)_.dat", hcat(rl, Float64.(abs.(r1))))
            writedlm("test_store/RadTest_LowA_2_n_$(n_lvls[i])_l_$(l)_m_$(m)_.dat", hcat(x_values, Float64.(abs.(y_values))))
        elseif j == length(alph)
            writedlm("test_store/RadTest_HighA_1_n_$(n_lvls[i])_l_$(l)_m_$(m)_.dat", hcat(rl, Float64.(abs.(r1))))
            writedlm("test_store/RadTest_HighA_2_n_$(n_lvls[i])_l_$(l)_m_$(m)_.dat", hcat(x_values, Float64.(abs.(y_values))))
        end
        
    end
    writedlm("test_store/ErgEvol_n_$(n_lvls[i])_l_$(l)_m_$(m)_.dat", hcat(alph, Float64.(outPr)))
    writedlm("test_store/ErgEvol_I_n_$(n_lvls[i])_l_$(l)_m_$(m)_.dat", hcat(alph, Float64.(outPi)))
    
    writedlm("test_store/ErgEvol_Cheb_n_$(n_lvls[i])_l_$(l)_m_$(m)_.dat", hcat(alph, Float64.(outPr_2)))
    writedlm("test_store/ErgEvol_Cheb_I_n_$(n_lvls[i])_l_$(l)_m_$(m)_.dat", hcat(alph, Float64.(outPi_2)))
    
end
