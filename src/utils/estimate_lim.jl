using Random
using Distributions
using DelimitedFiles
using Suppressor
@suppress include("super_rad.jl")
include("tde_input.jl")
using Dates
using ArgParse
using StatsBase

fname = "estim_lim_GRS_NR"
# fname = "estim_lim_Cyg_NR"
# fname = "estim_lim_Cyg_Rel"

function run()
        
    # massL = 10 .^ LinRange(log10.(1e-12), log10.(6e-12), 10)
    massL = 10 .^ LinRange(log10.(2e-12), log10.(6e-12), 7)
    


    runs = 10
    
    # m_a = 2.7e-13
    
#    SpinLim = 0.996
#    MassBH = 21.2
#    Mass_u = 0.10
#    tau_max = 4.8e6
    
    SpinLim = 0.95
    MassBH = 12.4
    Mass_u = 2.00
    tau_max = 5e7

    d_mass = Normal(MassBH, Mass_u)
    
    athresh = 1e-4
    high_p = true
    alpha_max_cut = 10.0
    impose_low_cut = 1e-3
    return_all_info = false
   
    n_times = 1000000
    # n_times = 10000
    
    eq_threshold=1e-100
    stop_on_a = 0.9
    abstol=1e-30
    
    non_rel = false
    print("Non rel \t ", non_rel, "\n")
    Nmax = 3
    cheby = false
    
    debug=false
    out_massL = []
    for m_a in massL
        alph = GNew .* MassBH .* m_a
        maxa = 4 .* alph ./ (1 .+ 4 .* alph.^2)
        # print("Mass \t ", m_a, "\n")

    
        out_store = []
        for i in 1:runs
            print(i, " ")
            
            MassIn = rand(d_mass,1)[1]
            # spinIn = rand(Uniform(SpinLim, 0.998))
            spinIn = rand(Uniform(0.99799, 0.998))
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
        
        final_lim = 10 .^ mean(log10.(out_store))
        # final_lim = 10 .^ percentile(log10.(out_store), 90)
        print(alph, "\t", m_a, "\t", final_lim, "\n")
        
        out_massL = cat(out_massL, [m_a final_lim], dims=1)
    end
    return out_massL
end

out = run()
println(out)
writedlm("test_store/"*fname*".dat", out)
