include("solve_sr_rates.jl")
using DelimitedFiles

mu = 1e-12
Mmin = 1
Mmax = 200
Num_samps = 50000

avals = rand(Num_samps)
Mvals = 10.0 .^ (rand(Num_samps) .* (log10.(Mmax) .- log10.(Mmin)) .+ log10.(Mmin) ) # rand(Num_samps) * (Mmax - Mmin) .+ Mmin

tscale = 1e10
ftag = "_tau_1e10_"

nmax=7
outarr = zeros(Num_samps, 2)
for i in 1:Num_samps
    age = 0.0
    alp = GNew .* Mvals[i] .* mu

    spun_down = false
    nn = 2
    atemp = avals[i]
    while !spun_down
        tau = 1.0 ./ (2 .*  sr_rates(nn, nn-1, nn-1, alp ./ (GNew .* Mvals[i]), Mvals[i], atemp) .* alp)  .* 6.58e-16 ./ 3.15e7 .* 180
#         tau = 1.0 ./ (2 .* find_im_part( mu , Mvals[i], atemp, nn, nn-1, nn-1) .* alp) .* 6.58e-16 ./ 3.15e7 .* 180
#        println(atemp, "\t", alp, "\t", nn, "\t", age, "\t", tau)
        if tau > 0.0
           age += tau
           if age > tscale
               spun_down = true
           else
               ahold = 4 .* alp .* (nn - 1.0) ./ ((nn .- 1.0).^2 .+ 4 .* alp.^2)
               if ahold .< atemp
                   atemp = ahold
               end
               nn += 1
           end
        else
            nn += 1
        end
        
        if nn > nmax
           spun_down = true
        end
    end
    outarr[i,:] .= [Mvals[i], atemp]
#    println("")
end

writedlm("test_store/Regge_Eg_"*ftag*".dat", outarr)
