include("solve_sr_rates.jl")

n_lvls = [1 2 3 4 5]
alph = LinRange(0.03, 1, 5000)

l = 0
m = 0
M = 22.0
a = 0.9

for i in 1:length(n_lvls)
    print(n_lvls[i], "\n")
    outPr = []
    outPi = []
    for j in 1:length(alph)
        print("\t\t", alph[j], "\n\n")
    	mu = alph[j] ./ (GNew * M)
        if alph[j] < 0.1
            wR, wI = find_im_part(mu, M, a, n_lvls[i], l, m; debug=true, iter=10000, xtol=1e-20, ftol=1e-90, return_both=true, for_s_rates=true)
        else
            itpR = LinearInterpolation(alph[1:j-1], outPr, extrapolation_bc=Line())
            itpI = LinearInterpolation(alph[1:j-1], outPi, extrapolation_bc=Line())
            
            guessV = itpR(alph[j]) .+ itpI(alph[j]) * im
            print(alph[j], "\t", guessV, "\n")
            wR, wI = find_im_part(mu, M, a, n_lvls[i], l, m; debug=true, iter=10000, xtol=1e-20, ftol=1e-90, return_both=true, for_s_rates=true, erg_Guess=guessV)
        end
	append!(outPr, wR)
	append!(outPi, wI)
    end
    writedlm("test_store/ErgEvol_n_$(n_lvls[i])_l_$(l)_m_$(m)_.dat", hcat(alph, Float64.(outPr)))
    writedlm("test_store/ErgEvol_I_n_$(n_lvls[i])_l_$(l)_m_$(m)_.dat", hcat(alph, Float64.(outPi)))
end
