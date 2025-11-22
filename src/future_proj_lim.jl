using Random
using Distributions
using DelimitedFiles
using Suppressor
@suppress include("super_rad.jl")
include("tde_input.jl")
using Dates
using ArgParse
using StatsBase

fname = "proj_optimistic_lim"

function run()
        
    # MassBH_list = 10 .^ LinRange(log10.(5.0), log10.(1e9), 10)
    
    massL = 10 .^ LinRange(log10.(1e-20), log10.(2e-11), 10)
    alphaL = LinRange(0.35, 0.41, 5)
    SpinLim = 0.97
    tau_max = 5e7
    

    athresh = 1e-4
    high_p = true
    alpha_max_cut = 10.0
    impose_low_cut = 1e-3
    return_all_info = false
    
    n_times = 1000000
    
    eq_threshold=1e-100
    stop_on_a = 0.9
    abstol=1e-30
    non_rel = false
    spinIn = 0.998
    Nmax = 3
    cheby = false
    
    
    debug=false
    out_massL = []
    for m_a in massL
        
        out_store = []
        for i in 1:length(alphaL)
            
            MassIn = alphaL[i] ./ (GNew .* m_a)
            if MassIn < 3.3
                MassIn = 3.3
            end
            if MassIn > 1e10
                MassIn = 1e10
            end
        
            f_guess = 1.0e14
            f_max = 1e18
            f_min = 1e10
            
            foundLim = false
            cnt = 0
            while !foundLim

                spin, massf = solve_system(m_a, f_guess, spinIn, MassIn, tau_max; impose_low_cut=impose_low_cut, return_all_info=return_all_info, n_times=n_times, eq_threshold=eq_threshold, stop_on_a=stop_on_a, debug=debug, abstol=abstol, non_rel=non_rel, high_p=high_p, Nmax=Nmax, cheby=cheby)
                
                # print("test \t ", f_guess, "\t", spin, "\n")
                if spin .< SpinLim
                    f_max = f_guess
                    # f_guess = 10 .^ (0.5 * (log10.(f_guess) + log10.(f_min)))
                    f_guess = 10 .^ rand(Uniform(log10.(f_min), log10.(f_guess)))
                else
                    f_min = f_guess
                    # f_guess = 10 .^ (0.5 * (log10.(f_guess) + log10.(f_max)))
                    f_guess = 10 .^ rand(Uniform(log10.(f_guess), log10.(f_max)))
                end
                
                if abs.(SpinLim - spin) .< athresh
                    foundLim = true
                    append!(out_store, f_guess)
                    # println(out_store)
                end
                
                if cnt > 20
                    foundLim = true
                    append!(out_store, f_guess)
                    # println(out_store)
                end
                
                if spin > 9.9e17
                    foundLim = true
                    append!(out_store, 1e18)
                end
                cnt+=1
                
            end
        end
        
        final_lim = minimum(out_store)
        print(m_a, "\t", final_lim, "\n")
        
        out_massL = cat(out_massL, [m_a final_lim], dims=1)
    end
    return out_massL
end

out = run()
println(out)
writedlm("test_store/"*fname*".dat", out)
