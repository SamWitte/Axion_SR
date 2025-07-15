using Random
using Distributions
import AffineInvariantMCMC
using Turing
using StatsPlots
using DelimitedFiles
using Suppressor
@suppress include("super_rad.jl")
using MCMCDiagnosticTools
using Dates
using KernelDensity

function profileL_func_minimize(data, mass_ax, Fname, Nsamples; fa_min=1e11, fa_max=1e18, tau_max=1e4, non_rel=true, Nmax=3, cheby=false, delt_M=0.05, thinning=1, numsamples_perwalker=2000, burnin=500)
    
    ## data format: [M_1, M_2, chi_1, chi_2] samples
    
    function llhood(x)
        ff = log_probability(x, data, mass_ax, fa_min=fa_min, fa_max=fa_max, tau_max=tau_max, non_rel=non_rel, Nmax=Nmax, cheby=cheby, Nsamples=Nsamples)
        println(x, "\t", ff, "\n")
        return ff
    end
    
    @model function scalar_model()
        θ ~ Uniform(log10.(fa_min), log10.(fa_max))              # Prior on θ
        Turing.@addlogprob! llhood(θ)      # Custom log-likelihood
    end
    
    chain = sample(scalar_model(), MH(), MCMCThreads(), numsamples_perwalker, Threads.nthreads())
    println(chain)
    println(describe(chain))
    
    posterior_samples = chain[burnin+1:end, :, :]
    
    samples = Array(chain)
    
    println(samples)
    println(size(samples))

    writedlm("output_mcmc/"*Fname*"_mcmc.dat", samples)
    
end
    
    

function log_probability(theta, data, mass_ax; fa_min=1e11, fa_max=1e20, tau_max=5e7, non_rel=false, Nmax=3, cheby=false, delt_M=0.05, Nsamples=3)
    
    return log_likelihood(theta, data, mass_ax, tau_max=tau_max, non_rel=non_rel, Nmax=Nmax, cheby=cheby, delt_M=delt_M, Nsamples=Nsamples)
end

function sample_spin(SpinBH_c, Spin_errD)

    p_or_neg = Int(round(rand()))
    SpinBH = nothing
    val_found = false
    if p_or_neg == 0
        if maxSpin > SpinBH_c
            SpinBH = rand(Uniform(SpinBH_c, maxSpin))
        else
            SpinBH = maxSpin
        end
    else
        d_spin = Normal(SpinBH_c, Spin_errD)
        while !val_found
            SpinBH = rand(d_spin,1)[1]
            if (SpinBH .<= SpinBH_c) && (SpinBH .>= 0.0)
                val_found = true
            end
        end
    end
    return SpinBH
end

function log_likelihood(theta, data, mass_ax; tau_max=1e4, non_rel=true, debug=false, Nmax=3, cheby=false, Nsamples=3, delt_M=0.05)
    
    log_f = theta
    sum_loglike = 0.0
    maxI = size(data)[1]
    idx_hold = rand(1:maxI, Nsamples)
    
    sampled_data = data[idx_hold, :]
    
    
    for i in 1:Nsamples
        M1 = sampled_data[i, 1]
        M2 = sampled_data[i, 2]
        
        kde_data = data[(data[:, 1] .>= M1 .* (1.0 .- delt_M)) .& (data[:, 1] .<= M1 .* (1.0 .+ delt_M)) .& (data[:, 2] .>= M2 .* (1.0 .- delt_M)) .& (data[:, 2] .<= M2 .* (1.0 .+ delt_M)), :]
        
        
        spin1_c = mean(kde_data[:, 3])
        spin2_c = mean(kde_data[:, 4])
        spin1_d = std(kde_data[:, 3])
        spin2_d = std(kde_data[:, 4])
        s1 = sample_spin(spin1_c, spin1_d)
        s2 = sample_spin(spin2_c, spin2_d)
        
        alph1 = GNew .* M1 .* mass_ax #
        alph2 = GNew .* M2 .* mass_ax #
        
        
        final_spin_1, final_mass_1 = super_rad_check(M1, s1, mass_ax, 10 .^log_f[1], tau_max=tau_max, non_rel=non_rel, Nmax=Nmax, cheby=cheby)
        final_spin_2, final_mass_2 = super_rad_check(M2, s2, mass_ax, 10 .^log_f[1], tau_max=tau_max, non_rel=non_rel, Nmax=Nmax, cheby=cheby)
        
        println("Init / Final BH 1: ", M1, " ", s1, " ", final_mass_1, " ", final_spin_1)
        println("Init / Final BH 2: ", M2, " ", s2, " ", final_mass_2, " ", final_spin_2)
        
        # Perform 2D kernel density estimation
        kdespins = kde((kde_data[:,3], kde_data[:, 4]))
        
        likelihood = pdf(kdespins, final_spin_1, final_spin_2)
        if likelihood .> 0.0
            sum_loglike += log.(likelihood)
        else
            sum_loglike += -Inf
        end
            
    end

    return sum_loglike
end




