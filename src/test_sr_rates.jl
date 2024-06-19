include("solve_sr_rates.jl")
using DelimitedFiles

function main_plt()

    xtol = 1e-30
    ftol = 1e-30
    iter = 500
    # alist = range(0.5, 0.99, 7)
    alist = [0.5, 0.7, 0.9, 0.99]
    
    nn = 4
    ll = 3
    mm = 3
    #  amax = 4 * m * alpha ./ (m^2 .+ 4 * alpha.^2)
    
    # muL = range(0.01, 1.0, 50)
    

    for i in alist
        # muMax = i .* mm ./ (2 .+ 2 .* sqrt.(1 .- i.^2))
        muMax = (2 .* (1 + sqrt.(1.0 .- i^2)) .* nn.^2 .- sqrt.(8 * (1.0 .+ sqrt.(1.0 .- i.^2)) .* nn.^4 .- 2.0 .* i^2 * nn^2 * (mm.^2 .+ 2 * nn^2))) ./ (i .* mm)
        # print(muMax, "\n")
        muMin = 3e-2
        muL = muMax .- 10 .^ range(-4, log10.(muMax .* (1.0 .- muMin ./ muMax)), 20)
        
        lv211_d = []
        lv211_c = []
        
        for j in muL
            
            val = sr_rates(nn, ll, mm, j ./ (GNew), 1, i, impose_low_cut=1e-5) * (GNew)
            val2 = find_im_part(j ./ (GNew), 1, i, nn, ll, mm, iter=iter, xtol=xtol, ftol=ftol)
            if val > 0
                if length(lv211_d) > 0
                    lv211_d = [lv211_d; [j val ./ (j ./ (GNew))]]
                else
                    lv211_d = [j val ./ (j ./ (GNew))]
                end
            end
            
            if val2 > 0
                if length(lv211_c) > 0
                    lv211_c = [lv211_c; [j val2 ./ (j ./ (GNew))]]
                else
                    lv211_c = [j val2 ./ (j ./ (GNew))]
                end
            end
           

        end
        writedlm("test_store/$(nn)$(ll)$(mm)_det_spin_$(round(i, sigdigits=2))_test.dat", lv211_d)
        writedlm("test_store/$(nn)$(ll)$(mm)_new_spin_$(round(i, sigdigits=2))_test.dat", lv211_c)
       
#        writedlm("test_store/422_det_spin_$(round(i, sigdigits=2))_test.dat", lv422_d)
#        writedlm("test_store/422_new_spin_$(round(i, sigdigits=2))_test.dat", lv422_c)
#        writedlm("test_store/433_det_spin_$(round(i, sigdigits=2))_test.dat", lv433_d)
#        writedlm("test_store/433_new_spin_$(round(i, sigdigits=2))_test.dat", lv433_c)
        
        
#        writedlm("test_store/311_new_spin_$(round(i, sigdigits=2))_test_c.dat", lv311_c)
#        writedlm("test_store/321_new_spin_$(round(i, sigdigits=2))_test_c.dat", lv321_c)
#        writedlm("test_store/421_new_spin_$(round(i, sigdigits=2))_test_c.dat", lv421_c)
#        writedlm("test_store/431_new_spin_$(round(i, sigdigits=2))_test_c.dat", lv431_c)
#        writedlm("test_store/432_new_spin_$(round(i, sigdigits=2))_test_c.dat", lv432_c)
        
    end



end

@time main_plt()
