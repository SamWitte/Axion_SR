using Random
using Distributions
using DelimitedFiles
using Suppressor
@suppress include("super_rad.jl")
include("tde_input.jl")
using Dates

fout = "_"
sve = false


f_a = 1e15
m_a = 1e-19
SpinBH = 0.9
MassBH = 4.0e8
tau_max = 5.0e7 # 5.0e7, 4.8e6

alpha_max_cut = 10.0
impose_low_cut = 1e-100
return_all_info = true
n_times = 100000
# n_times = 10000
eq_threshold=1e-100
stop_on_a = 0.0
abstol=1e-30
N_pts_interp=200
N_pts_interpL=200


Nmax = 5
cheby=true

non_rel = false
high_p = true

alph = GNew .* MassBH .* m_a
maxa = 4 .* alph ./ (1 .+ 4 .* alph.^2)


print("Alpha and max a spin down \t", alph, "\t", maxa, "\n")
println("fa ", f_a)
println("non rel? ", non_rel)


fname = "test_store/OutRun_"

if non_rel
    fname *= "NonRel_"
else
    fname *= "FullRel_"
end
fname *= "fa_$(f_a)_ma_$(m_a)_MBH_$(MassBH)_spin_$(SpinBH)"*fout*".dat"


timeT, StatesOut, idx_lvl, spin, massB = @time solve_system(m_a, f_a, SpinBH, MassBH, tau_max; impose_low_cut=impose_low_cut, return_all_info=return_all_info, n_times=n_times, eq_threshold=eq_threshold, abstol=abstol, non_rel=non_rel, debug=true, high_p=high_p, N_pts_interp=N_pts_interp, N_pts_interpL=N_pts_interpL, Nmax=Nmax, cheby=cheby)
 

if sve
    if !solve_n4
        outSve = zeros(length(timeT), 6)
        for i in 1:length(timeT)
            outSve[i, :] = [timeT[i] state211[i] state311[i] state322[i] spin[i] massB[i]]
        end
        writedlm(fname, outSve)
    else
        if !solve_n5
            outSve = zeros(length(timeT), 9)
            for i in 1:length(timeT)
                outSve[i, :] = [timeT[i] state211[i] state311[i] state322[i] state411[i] state422[i] state433[i] spin[i] massB[i]]
            end
            writedlm(fname, outSve)
        else
            outSve = zeros(length(timeT), 13)
            for i in 1:length(timeT)
                outSve[i, :] = [timeT[i] state211[i] state311[i] state322[i] state411[i] state422[i] state433[i] state511[i] state522[i] state533[i] state544[i] spin[i] massB[i]]
            end
            writedlm(fname, outSve)
        end
    end

end
