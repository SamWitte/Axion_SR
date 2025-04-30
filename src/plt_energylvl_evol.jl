include("solve_sr_rates.jl")
include("cheby_eigs.jl")

# n_lvls = [2 3 4 5]
n_lvls = [3 4 5]
alph = LinRange(0.03, 1, 20)

l = 2
m = 2
M = 22.0
a = 0.9

Npts_Bnd=4000
rmaxT=100
Ntot_safe=7000

Npoints = 50
Iter = 20
cvg_acc = 1e-3
Npts_r = 1000
prec=200

for i in 1:length(n_lvls)
    print(n_lvls[i], "\n")
    outPr = []
    outPi = []
    
    outPr_2 = []
    outPi_2 = []
    
    for j in 1:length(alph)
        print("\t\t", alph[j], "\n\n")
    	mu = alph[j] ./ (GNew * M)
        
        # wR, wI = find_im_part(mu, M, a, n_lvls[i], l, m; debug=false, iter=100, xtol=1e-20, ftol=1e-90, return_both=true, for_s_rates=true)
        
        rl, r1, erg = solve_radial(mu, M, a, n_lvls[i], l, m; rpts=Npts_Bnd, rmaxT=rmaxT, return_erg=true, Ntot_safe=Ntot_safe, iter=100, xtol=1e-20, ftol=1e-90)
        wR = real.(erg)
        wI = imag.(erg)
        # rl, r1 = solve_radial(mu, M, a, n_lvls[i], l, m; rpts=Npts_Bnd, rmaxT=rmaxT, pre_compute_erg=(wR .+ im .* wI), Ntot_safe=Ntot_safe)
        println("First go \t ", wR, "\t", wI)
        append!(outPr, wR)
        append!(outPi, wI)
        
        
        
        
        wR, wI, x_values, y_values, nu = eigensys_Cheby(M, a, alph[j] ./ (GNew .* M), n_lvls[i], l, m, debug=false, return_wf=true, Npoints=Npoints, Iter=Iter, cvg_acc=cvg_acc, Npts_r=Npts_r, return_nu=true, prec=prec)
        rad_outT, routT = solve_radial(mu, M, a, n_lvls[i], l, m; rpts=Npts_Bnd, rmaxT=rmaxT, pre_compute_erg=(wR .+ im .* wI), Ntot_safe=Ntot_safe)
        # wRIG, wIIG, rad_outT, routT = eigensys_Cheby(M, a, alph[j] ./ (GNew .* M), n_lvls[i], l, m, debug=false, return_wf=true, Npoints=130, Iter=1, cvg_acc=cvg_acc, Npts_r=Npts_r, nu_guess=nu, prec=prec)
        
        println("Cheb go \t ", wR, "\t", wI)
        append!(outPr_2, wR)
        append!(outPi_2, wI)
        
        
        if j == 1
            println("Sving WF....")
            writedlm("test_store/RadTest_LowA_1_n_$(n_lvls[i])_l_$(l)_m_$(m)_.dat", hcat(rl, Float64.(abs.(r1))))
            writedlm("test_store/RadTest_LowA_2_n_$(n_lvls[i])_l_$(l)_m_$(m)_.dat", hcat(x_values, Float64.(abs.(y_values))))
            writedlm("test_store/RadTest_LowA_2v2_n_$(n_lvls[i])_l_$(l)_m_$(m)_.dat", hcat(rad_outT, Float64.(abs.(routT))))
        elseif j == 5
            println("Sving WF....")
            writedlm("test_store/RadTest_HighA_1_n_$(n_lvls[i])_l_$(l)_m_$(m)_.dat", hcat(rl, Float64.(abs.(r1))))
            writedlm("test_store/RadTest_HighA_2_n_$(n_lvls[i])_l_$(l)_m_$(m)_.dat", hcat(x_values, Float64.(abs.(y_values))))
            writedlm("test_store/RadTest_HighA_2v2_n_$(n_lvls[i])_l_$(l)_m_$(m)_.dat", hcat(rad_outT, Float64.(abs.(routT))))
        end
        
    end
    writedlm("test_store/ErgEvol_n_$(n_lvls[i])_l_$(l)_m_$(m)_.dat", hcat(alph, Float64.(outPr)))
    writedlm("test_store/ErgEvol_I_n_$(n_lvls[i])_l_$(l)_m_$(m)_.dat", hcat(alph, Float64.(outPi)))
    
    writedlm("test_store/ErgEvol_Cheb_n_$(n_lvls[i])_l_$(l)_m_$(m)_.dat", hcat(alph, Float64.(outPr_2)))
    writedlm("test_store/ErgEvol_Cheb_I_n_$(n_lvls[i])_l_$(l)_m_$(m)_.dat", hcat(alph, Float64.(outPi_2)))
    
end
