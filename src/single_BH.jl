using Random
using Distributions
using DelimitedFiles
include("super_rad.jl")
include("tde_input.jl")
using Dates

fout = ""
sve = true

f_a = 1e15
m_a = 3e-13
MassBH = 14.8
SpinBH = 0.99
alpha_max_cut = 1.0
impose_low_cut = 1e-3
tau_max = 5e6
return_all_info = true


alph = GNew .* MassBH .* m_a
print(alph, "\n")

time, state211, state322, spin, massB = solve_system(m_a, f_a, SpinBH, MassBH, tau_max; impose_low_cut=impose_low_cut, return_all_info=return_all_info)

if sve
    outSve = zeros(length(time), 5)
    for i in 1:length(time)
        outSve[i, :] = [time[i] state211[i] state322[i] spin[i] massB[i]]
    end
    writedlm("test_store/OutRun_fa_$(f_a)_ma_$(m_a)_MBH_$(MassBH)_spin_$(SpinBH).dat", outSve)

end
