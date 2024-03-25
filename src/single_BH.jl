using Random
using Distributions
using DelimitedFiles
include("super_rad.jl")
include("tde_input.jl")
using Dates

fout = ""
sve = false


f_a = 1e18
m_a = 1e-12
MassBH = 14.8
SpinBH = 0.9
alpha_max_cut = 1.0
impose_low_cut = 1e-3
tau_max = 5e6
return_all_info = true
solve_322 = true
n_times = 10000
input_data="Masha"

alph = GNew .* MassBH .* m_a
maxa = 4 .* alph ./ (1 .+ 4 .* alph.^2)
print(alph, "\t", maxa, "\n")

solve_n4 = true
if !solve_n4
    time, state211, state322, spin, massB = @time solve_system(m_a, f_a, SpinBH, MassBH, tau_max; impose_low_cut=impose_low_cut, return_all_info=return_all_info, solve_322=solve_322, n_times=n_times, input_data=input_data)
else
    time, state211, state322, state411, spin, massB = @time solve_system(m_a, f_a, SpinBH, MassBH, tau_max; impose_low_cut=impose_low_cut, return_all_info=return_all_info, solve_322=solve_322, n_times=n_times, input_data=input_data, solve_n4=true)
end

if sve
    outSve = zeros(length(time), 5)
    for i in 1:length(time)
        outSve[i, :] = [time[i] state211[i] state322[i] spin[i] massB[i]]
    end
    writedlm("test_store/OutRun_fa_$(f_a)_ma_$(m_a)_MBH_$(MassBH)_spin_$(SpinBH).dat", outSve)

end
