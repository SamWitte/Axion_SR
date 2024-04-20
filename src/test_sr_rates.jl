include("solve_sr_rates.jl")
using DelimitedFiles

function main_plt()
    Ntot = 2000
    xtol = 1e-30
    ftol = 1e-20
    iter = 100
    # alist = range(0.5, 0.99, 7)
    alist = [0.5, 0.7, 0.9, 0.99]
    # alist = [0.99]
    muL = range(0.01, 1.0, 200)

    for i in alist
        lv211_d = []
        lv211_c = []

        lv411_d = []
        lv411_c = []

        lv322_d = []
        lv322_c = []
        
        lv311_c = []
        lv321_c = []
        lv421_c = []
        lv431_c = []
        lv432_c = []
        lv433_c = []
        lv433_d = []
        lv422_c = []
        lv422_d = []
        
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

            
            # val2 = find_im_part(j ./ (GNew), 1, i, 3, 1, 1, Ntot=Ntot, iter=iter, xtol=xtol)
            val2 = sr_rates(3, 1, 1, j ./ (GNew), 1, i, impose_low_cut=1e-5) * (GNew)
            if val2 > 0
                if length(lv311_c) > 0
                    lv311_c = [lv311_c; [j val2 ./ (j ./ (GNew))]]
                else
                    lv311_c = [j val2 ./ (j ./ (GNew))]
                end
            end

            # val2 = find_im_part(j ./ (GNew), 1, i, 3, 2, 1, Ntot=Ntot, iter=iter, xtol=xtol)
            val2 = sr_rates(3, 2, 1, j ./ (GNew), 1, i, impose_low_cut=1e-5) * (GNew)
            if val2 > 0
                if length(lv321_c) > 0
                    lv321_c = [lv321_c; [j val2 ./ (j ./ (GNew))]]
                else
                    lv321_c = [j val2 ./ (j ./ (GNew))]
                end
            end

            # val2 = find_im_part(j ./ (GNew), 1, i, 4, 2, 1, Ntot=Ntot, iter=iter, xtol=xtol)
            val2 = sr_rates(4, 2, 1, j ./ (GNew), 1, i, impose_low_cut=1e-5) * (GNew)
            if val2 > 0
                if length(lv421_c) > 0
                    lv421_c = [lv421_c; [j val2 ./ (j ./ (GNew))]]
                else
                    lv421_c = [j val2 ./ (j ./ (GNew))]
                end
            end


            # val2 = find_im_part(j ./ (GNew), 1, i, 4, 3, 1, Ntot=Ntot, iter=iter, xtol=xtol)
            val2 = sr_rates(4, 3, 1, j ./ (GNew), 1, i, impose_low_cut=1e-5) * (GNew)
            if val2 > 0
                if length(lv431_c) > 0
                    lv431_c = [lv431_c; [j val2 ./ (j ./ (GNew))]]
                else
                    lv431_c = [j val2 ./ (j ./ (GNew))]
                end
            end

            # val2 = find_im_part(j ./ (GNew), 1, i, 4, 3, 2, Ntot=Ntot, iter=iter, xtol=xtol)
            val2 = sr_rates(4, 3, 2, j ./ (GNew), 1, i, impose_low_cut=1e-5) * (GNew)
            if val2 > 0
                if length(lv432_c) > 0
                    lv432_c = [lv432_c; [j val2 ./ (j ./ (GNew))]]
                else
                    lv432_c = [j val2 ./ (j ./ (GNew))]
                end
            end
            
            
            val2 = sr_rates(4, 3, 3, j ./ (GNew), 1, i, impose_low_cut=1e-5) * (GNew)
            if val2 > 0
                if length(lv433_d) > 0
                    lv433_d = [lv433_d; [j val2 ./ (j ./ (GNew))]]
                else
                    lv433_d = [j val2 ./ (j ./ (GNew))]
                end
            end
            val2 = find_im_part(j ./ (GNew), 1, i, 4, 3, 3, Ntot=Ntot, iter=iter, xtol=xtol, ftol=ftol)
            if val2 > 0
                if length(lv433_c) > 0
                    lv433_c = [lv433_c; [j val2 ./ (j ./ (GNew))]]
                else
                    lv433_c = [j val2 ./ (j ./ (GNew))]
                end
            end
            val2 = sr_rates(4, 2, 2, j ./ (GNew), 1, i, impose_low_cut=1e-5) * (GNew)
            if val2 > 0
                if length(lv422_d) > 0
                    lv422_d = [lv422_d; [j val2 ./ (j ./ (GNew))]]
                else
                    lv422_d = [j val2 ./ (j ./ (GNew))]
                end
            end
            val2 = find_im_part(j ./ (GNew), 1, i, 4, 2, 2, Ntot=Ntot, iter=iter, xtol=xtol, ftol=ftol)
            if val2 > 0
                if length(lv422_c) > 0
                    lv422_c = [lv422_c; [j val2 ./ (j ./ (GNew))]]
                else
                    lv422_c = [j val2 ./ (j ./ (GNew))]
                end
            end
            
        end
        writedlm("test_store/211_det_spin_$(round(i, sigdigits=2))_test.dat", lv211_d)
        writedlm("test_store/211_new_spin_$(round(i, sigdigits=2))_test.dat", lv211_c)
        writedlm("test_store/411_det_spin_$(round(i, sigdigits=2))_test.dat", lv411_d)
        writedlm("test_store/411_new_spin_$(round(i, sigdigits=2))_test.dat", lv411_c)
        writedlm("test_store/322_det_spin_$(round(i, sigdigits=2))_test.dat", lv322_d)
        writedlm("test_store/322_new_spin_$(round(i, sigdigits=2))_test.dat", lv322_c)
        
        writedlm("test_store/422_det_spin_$(round(i, sigdigits=2))_test.dat", lv422_d)
        writedlm("test_store/422_new_spin_$(round(i, sigdigits=2))_test.dat", lv422_c)
        writedlm("test_store/433_det_spin_$(round(i, sigdigits=2))_test.dat", lv433_d)
        writedlm("test_store/433_new_spin_$(round(i, sigdigits=2))_test.dat", lv433_c)
        
        
        writedlm("test_store/311_new_spin_$(round(i, sigdigits=2))_test_c.dat", lv311_c)
        writedlm("test_store/321_new_spin_$(round(i, sigdigits=2))_test_c.dat", lv321_c)
        writedlm("test_store/421_new_spin_$(round(i, sigdigits=2))_test_c.dat", lv421_c)
        writedlm("test_store/431_new_spin_$(round(i, sigdigits=2))_test_c.dat", lv431_c)
        writedlm("test_store/432_new_spin_$(round(i, sigdigits=2))_test_c.dat", lv432_c)
        
    end



end

main_plt()
