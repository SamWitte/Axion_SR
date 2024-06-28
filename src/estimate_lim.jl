using Random
using Distributions
using DelimitedFiles
include("super_rad.jl")
include("tde_input.jl")
using Dates






function run()
        
    
    f_guess = 1.0e14
    f_max = 1e18
    f_min = 1e12


    m_a = 2.7e-13
    SpinBH = 0.998
    SpinLim = 0.996
    MassBH = 22.2
    tau_max = 5e6
    athresh = 0.001

    alpha_max_cut = 10.0
    impose_low_cut = 1e-3
    return_all_info = false
    solve_322 = true
    n_times = 1000000
    # n_times = 10000
    input_data="Me"
    eq_threshold=1e-100
    stop_on_a = 0.7
    abstol=1e-25
    solve_n4 = true
    non_rel = true
    
    alph = GNew .* MassBH .* m_a
    maxa = 4 .* alph ./ (1 .+ 4 .* alph.^2)
    debug=false
    print("Mass \t ", m_a, "\n")

    foundLim = false
    while !foundLim

        spin = solve_system(m_a, f_guess, SpinBH, MassBH, tau_max; impose_low_cut=impose_low_cut, return_all_info=return_all_info, solve_322=solve_322, n_times=n_times, input_data=input_data, eq_threshold=eq_threshold, stop_on_a=stop_on_a, debug=debug, abstol=abstol, solve_n4=solve_n4, non_rel)
        
        print(f_guess, "\t", spin, "\n")
        if spin .< SpinLim
            f_max = f_guess
            f_guess = 10 .^ (0.5 * (log10.(f_guess) + log10.(f_min)))
        else
            f_min = f_guess
            f_guess = 10 .^ (0.5 * (log10.(f_guess) + log10.(f_max)))
        end
        
        if abs.(SpinLim - spin) .< athresh
            foundLim = true
        end
        
    end

end


run()
