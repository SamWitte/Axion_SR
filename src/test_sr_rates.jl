include("solve_sr_rates.jl")
using DelimitedFiles

function main_plt()
    Ntot = 2000
    xtol = 1e-20
    iter = 50
    # alist = range(0.5, 0.99, 7)
    alist = [0.5, 0.7, 0.9, 0.99]
    muL = range(0.01, 1.0, 200)

    for i in alist
        lv211_d = []
        lv211_c = []

        lv411_d = []
        lv411_c = []

        lv322_d = []
        lv322_c = []
        
        for j in muL
            val = sr_rates(2, 1, 1, j ./ (GNew), 1, i, impose_low_cut=1e-5) * (GNew)
            val2 = find_im_part(j ./ (GNew), 1, i, 2, 1, 1, Ntot=Ntot, iter=iter, xtol=xtol)
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
            
            val = sr_rates(4, 1, 1, j ./ (GNew), 1, i, impose_low_cut=1e-5) * (GNew)
            val2 = find_im_part(j ./ (GNew), 1, i, 4, 1, 1, Ntot=Ntot, iter=iter, xtol=xtol)
            if val > 0
                if length(lv411_d) > 0
                    lv411_d = [lv411_d; [j val ./ (j ./ (GNew))]]
                else
                    lv411_d = [j val ./ (j ./ (GNew))]
                end
            end
            if val2 > 0
                if length(lv411_c) > 0
                    lv411_c = [lv411_c; [j val2 ./ (j ./ (GNew))]]
                else
                    lv411_c = [j val2 ./ (j ./ (GNew))]
                end
            end
            
            
            val = sr_rates(3, 2, 2, j ./ (GNew), 1, i, impose_low_cut=1e-5) * (GNew)
            val2 = find_im_part(j ./ (GNew), 1, i, 3, 2, 2, Ntot=Ntot, iter=iter, xtol=xtol)
            if val > 0
                if length(lv322_d) > 0
                    lv322_d = [lv322_d; [j val ./ (j ./ (GNew))]]
                else
                    lv322_d = [j val ./ (j ./ (GNew))]
                end
            end
            if val2 > 0
                if length(lv322_c) > 0
                    lv322_c = [lv322_c; [j val2 ./ (j ./ (GNew))]]
                else
                    lv322_c = [j val2 ./ (j ./ (GNew))]
                end
                
            end
            
        end
        writedlm("test_store/211_det_spin_$(round(i, sigdigits=2))_test.dat", lv211_d)
        writedlm("test_store/211_new_spin_$(round(i, sigdigits=2))_test.dat", lv211_c)
        writedlm("test_store/411_det_spin_$(round(i, sigdigits=2))_test.dat", lv411_d)
        writedlm("test_store/411_new_spin_$(round(i, sigdigits=2))_test.dat", lv411_c)
        writedlm("test_store/322_det_spin_$(round(i, sigdigits=2))_test.dat", lv322_d)
        writedlm("test_store/322_new_spin_$(round(i, sigdigits=2))_test.dat", lv322_c)
    end



end

main_plt()
