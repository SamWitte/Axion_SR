using Random
using Distributions
using DelimitedFiles
include("super_rad.jl")
include("tde_input.jl")
using Dates






function run()
        
    
    f_guess = 1.0e14
    f_max = 1e16
    f_min = 1e10


    m_a = 2e-13
    SpinBH = 0.998
    SpinLim = 0.996
    MassBH = 22.2
    tau_max = 5e6

    alpha_max_cut = 10.0
    impose_low_cut = 1e-3
    return_all_info = false
    solve_322 = true
    n_times = 1000000
    # n_times = 10000
    input_data="Me"
    eq_threshold=1e-100
    stop_on_a = 0.9
    abstol=1e-25
    solve_n4 = true

    alph = GNew .* MassBH .* m_a
    maxa = 4 .* alph ./ (1 .+ 4 .* alph.^2)
    debug=false

    foundLim = false
    while !foundLim

        spin = solve_system(m_a, f_guess, SpinBH, MassBH, tau_max; impose_low_cut=impose_low_cut, return_all_info=return_all_info, solve_322=solve_322, n_times=n_times, input_data=input_data, eq_threshold=eq_threshold, stop_on_a=stop_on_a, debug=debug, abstol=abstol, solve_n4=solve_n4)
        
        print(f_guess, "\t", spin, "\n")
        if spin .< SpinLim
            f_max = f_guess
            f_guess = 0.5 * (f_guess + f_min)
        else
            f_min = f_guess
            f_guess = 0.5 * (f_guess + f_max)
        end
        
        if abs.(SpinLim - spin) .< 0.001
            foundLim = true
        end
        
    end

end


run()
