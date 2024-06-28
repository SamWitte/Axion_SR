using Random
using Distributions
using DelimitedFiles
include("super_rad.jl")
include("tde_input.jl")
using Dates

fout = "_"
sve = true


f_a = 1e14
m_a = 5.0e-12 # 2.0e-12
SpinBH = 0.99
MassBH = 22.2
tau_max = 4e6

alpha_max_cut = 10.0
impose_low_cut = 1e-3
return_all_info = true
solve_322 = true
n_times = 1000000
# n_times = 10000
input_data="Me"
eq_threshold=1e-100
stop_on_a = 0.0
abstol=1e-50

solve_n4 = true
non_rel = false

alph = GNew .* MassBH .* m_a
maxa = 4 .* alph ./ (1 .+ 4 .* alph.^2)
print(alph, "\t", maxa, "\n")

fname = "test_store/OutRun_"
if solve_n4
    fname *= "Lvln4_"
end
if non_rel
    fname *= "NonRel_"
else
    fname *= "FullRel_"
end
fname *= "fa_$(f_a)_ma_$(m_a)_MBH_$(MassBH)_spin_$(SpinBH)"*fout*".dat"


if !solve_n4
    time, state211, state322, spin, massB = @time solve_system(m_a, f_a, SpinBH, MassBH, tau_max; impose_low_cut=impose_low_cut, return_all_info=return_all_info, solve_322=solve_322, n_times=n_times, input_data=input_data, eq_threshold=eq_threshold, stop_on_a=stop_on_a, debug=true, abstol=abstol, non_rel=non_rel)
    # print(spin[end], "\n")
else
    time, state211, state322, state411, state422, state433, spin, massB = @time solve_system(m_a, f_a, SpinBH, MassBH, tau_max; impose_low_cut=impose_low_cut, return_all_info=return_all_info, solve_322=solve_322, n_times=n_times, input_data=input_data, solve_n4=true, eq_threshold=eq_threshold, abstol=abstol, non_rel=non_rel)
end

if sve
    if !solve_n4
        outSve = zeros(length(time), 5)
        for i in 1:length(time)
            outSve[i, :] = [time[i] state211[i] state322[i] spin[i] massB[i]]
        end
        writedlm(fname, outSve)
    else
        outSve = zeros(length(time), 8)
        for i in 1:length(time)
            outSve[i, :] = [time[i] state211[i] state322[i] state411[i] state422[i] state433[i] spin[i] massB[i]]
        end
        writedlm(fname, outSve)
    end

end
