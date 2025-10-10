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

function profileL_func_minimize(data, mass_ax, Fname, Nsamples; fa_min=1e11, fa_max=1e18, tau_max=1e4, non_rel=true, Nmax=3, cheby=false, delt_M=0.05, thinning=1, numsamples_perwalker=2000, burnin=500, use_kde=true, over_run=true, high_p=true, one_BH=false, high_spin_cut=nothing)
    
    ## data format: [M_1, M_2, chi_1, chi_2] samples
    
    function llhood(x)
        ff = log_probability(x, data, mass_ax, fa_min=fa_min, fa_max=fa_max, tau_max=tau_max, non_rel=non_rel, Nmax=Nmax, cheby=cheby, Nsamples=Nsamples, use_kde=use_kde, high_p=high_p, one_BH=one_BH, high_spin_cut=high_spin_cut)
        println(x, "\t", ff, "\n")
        return ff
    end
    
    @model function scalar_model()
        θ ~ Uniform(log10.(fa_min), log10.(fa_max))              # Prior on θ
        Turing.@addlogprob! llhood(θ)      # Custom log-likelihood
    end
    
    chain = sample(scalar_model(), MH(), MCMCThreads(), numsamples_perwalker, Threads.nthreads())
    # println(chain)
    # println(describe(chain))
    
    posterior_samples = chain[burnin+1:end, :, :]
    
    samples = Array(chain)
    
    fout = "output_mcmc/"*Fname*"_mcmc.dat"
    if isfile(fout) && over_run
        fin = readdlm(fout)
        fnew = cat(fin, samples, dims=1)
        writedlm(fout, fnew)
    else
        writedlm(fout, samples)
    end
    
end
    
    

function log_probability(theta, data, mass_ax; fa_min=1e11, fa_max=1e20, tau_max=5e7, non_rel=false, Nmax=3, cheby=false, delt_M=0.05, Nsamples=3, use_kde=true, high_p=true, one_BH=false, high_spin_cut=nothing)
    
    return log_likelihood(theta, data, mass_ax, tau_max=tau_max, non_rel=non_rel, Nmax=Nmax, cheby=cheby, delt_M=delt_M, Nsamples=Nsamples, use_kde=use_kde, high_p=high_p, one_BH=one_BH, high_spin_cut=high_spin_cut)
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

function log_likelihood(theta, data, mass_ax; tau_max=1e4, non_rel=true, debug=false, Nmax=3, cheby=false, Nsamples=3, delt_M=0.05, use_kde=true, high_p=true, one_BH=false, high_spin_cut=nothing)
    
    log_f = theta
    sum_loglike = 0.0
    maxI = size(data)[1]
    idx_hold = rand(1:maxI, Nsamples)
    
    sampled_data = data[idx_hold, :]
    
    ## note to self: check if strong correlations in mass/spin, could be something to think about putting Mf on cut
    for i in 1:Nsamples
        M1 = sampled_data[i, 1]
        M2 = sampled_data[i, 2]
        
        delt_M_use = delt_M
        kde_data = []
        N_min_kde = 50
        if length(data[:,1]) <= N_min_kde
            N_min_kde = length(data[:,1]) ./ 2
        end
        
        while length(kde_data) < N_min_kde
            kde_data = data[(data[:, 1] .>= M1 .* (1.0 .- delt_M_use)) .& (data[:, 1] .<= M1 .* (1.0 .+ delt_M_use)) .& (data[:, 2] .>= M2 .* (1.0 .- delt_M_use)) .& (data[:, 2] .<= M2 .* (1.0 .+ delt_M_use)), :]
            delt_M_use *= 1.1
        end
        
        spin1_c = mean(kde_data[:, 3])
        spin2_c = mean(kde_data[:, 4])
        spin1_d = std(kde_data[:, 3])
        spin2_d = std(kde_data[:, 4])
        if isnothing(high_spin_cut)
            s1 = sample_spin(spin1_c, spin1_d)
            s2 = sample_spin(spin2_c, spin2_d)
        else
            s1 = 1.0
            s2 = 1.0
            while s1 > high_spin_cut
                s1 = sample_spin(spin1_c, spin1_d)
            end
            while s2 > high_spin_cut
                s2 = sample_spin(spin2_c, spin2_d)
            end
        end
        
        alph1 = GNew .* M1 .* mass_ax #
        alph2 = GNew .* M2 .* mass_ax #
        
        println("fa \t", 10 .^log_f[1])
        final_spin_1, final_mass_1 = super_rad_check(M1, s1, mass_ax, 10 .^log_f[1], tau_max=tau_max, non_rel=non_rel, Nmax=Nmax, cheby=cheby, high_p=high_p)
        println("Init / Final BH 1: ", M1, " ", s1, " ", final_mass_1, " ", final_spin_1)
        if !one_BH
            final_spin_2, final_mass_2 = super_rad_check(M2, s2, mass_ax, 10 .^log_f[1], tau_max=tau_max, non_rel=non_rel, Nmax=Nmax, cheby=cheby, high_p=high_p)
            println("Init / Final BH 2: ", M2, " ", s2, " ", final_mass_2, " ", final_spin_2)
        else
            final_spin_2 = s2
        end
        
        
        if use_kde
            # Perform 2D kernel density estimation
            kdespins = kde((kde_data[:,3], kde_data[:, 4]))
            likelihood = pdf(kdespins, final_spin_1, final_spin_2)
        else
            likelihood = exp.(-(final_spin_1 .- spin1_c).^2 ./ (2 .* spin1_d.^2))
            if !one_BH
                likelihood .*= exp.(-(final_spin_2 .- spin2_c).^2 ./ (2 .* spin2_d.^2))
            end
        end
        if likelihood .> 0.0
            sum_loglike += log.(likelihood)
        else
            sum_loglike += -Inf
        end
            
    end

    return sum_loglike
end

function log_probability_s1(theta, data; tau_max=5e7, delt_M=0.05, Nsamples=3, use_kde=true, high_p=true)
    return log_likelihood_spinone(theta, data, tau_max=tau_max, delt_M=delt_M, Nsamples=Nsamples, use_kde=use_kde)
end


function profileL_func_minimize_spinone(data, Fname, Nsamples; minmass=1e-15, maxmass=1e-12, tau_max=1e4, delt_M=0.05, thinning=1, numsamples_perwalker=2000, burnin=500, use_kde=true, over_run=true, high_p=true)
    
    ## data format: [M_1, M_2, chi_1, chi_2] samples
    
    function llhood(x)
        ff = log_probability_s1(x, data, tau_max=tau_max, Nsamples=Nsamples, use_kde=use_kde)
        return ff
    end
    
    @model function scalar_model()
        θ ~ Uniform(log10.(minmass), log10.(maxmass))              # Prior on θ
        Turing.@addlogprob! llhood(θ)      # Custom log-likelihood
    end
    
    chain = sample(scalar_model(), MH(), MCMCThreads(), numsamples_perwalker, Threads.nthreads())
    # println(chain)
    # println(describe(chain))
    
    posterior_samples = chain[burnin+1:end, :, :]
    
    samples = Array(chain)
    
    fout = "output_mcmc/"*Fname*"_mcmc.dat"
    if isfile(fout) && over_run
        fin = readdlm(fout)
        fnew = cat(fin, samples, dims=1)
        writedlm(fout, fnew)
    else
        writedlm(fout, samples)
    end
    
end

function log_likelihood_spinone(theta, data; tau_max=1e4, debug=false, Nsamples=3, delt_M=0.05, use_kde=true)
    
    mass_ax = 10 .^theta
    sum_loglike = 0.0
    maxI = size(data)[1]
    idx_hold = rand(1:maxI, Nsamples)
    
    sampled_data = data[idx_hold, :]
    
    ## note to self: check if strong correlations in mass/spin, could be something to think about putting Mf on cut
    for i in 1:Nsamples
        M1 = sampled_data[i, 1]
        M2 = sampled_data[i, 2]
        
        delt_M_use = delt_M
        kde_data = []
        N_min_kde = 50
        if length(data[:,1]) <= N_min_kde
            N_min_kde = length(data[:,1]) ./ 2
        end
        
        while length(kde_data) < N_min_kde
            kde_data = data[(data[:, 1] .>= M1 .* (1.0 .- delt_M_use)) .& (data[:, 1] .<= M1 .* (1.0 .+ delt_M_use)) .& (data[:, 2] .>= M2 .* (1.0 .- delt_M_use)) .& (data[:, 2] .<= M2 .* (1.0 .+ delt_M_use)), :]
            delt_M_use *= 1.1
        end
        
        spin1_c = mean(kde_data[:, 3])
        spin2_c = mean(kde_data[:, 4])
        spin1_d = std(kde_data[:, 3])
        spin2_d = std(kde_data[:, 4])
        s1 = sample_spin(spin1_c, spin1_d)
        s2 = sample_spin(spin2_c, spin2_d)
        
        
        
        
        alph1 = GNew .* M1 .* mass_ax #
        alph2 = GNew .* M2 .* mass_ax #
        
        
        final_spin_1, final_mass_1 = super_rad_check(M1, s1, mass_ax, 1.0, spinone=true, tau_max=tau_max)
        final_spin_2, final_mass_2 = super_rad_check(M2, s2, mass_ax, 1.0, spinone=true, tau_max=tau_max)
        
        println("Init / Final BH 1: ", M1, " ", s1, " ", final_mass_1, " ", final_spin_1)
        println("Init / Final BH 2: ", M2, " ", s2, " ", final_mass_2, " ", final_spin_2)
        
        if use_kde
            # Perform 2D kernel density estimation
            kdespins = kde((kde_data[:,3], kde_data[:, 4]))
            likelihood = pdf(kdespins, final_spin_1, final_spin_2)

        else
            likelihood = exp.(-(final_spin_1 .- spin1_c).^2 ./ (2 .* spin1_d.^2)) .* exp.(-(final_spin_2 .- spin2_c).^2 ./ (2 .* spin2_d.^2))
        end
        if likelihood .> 0.0
            sum_loglike += log.(likelihood)
        else
            sum_loglike += -Inf
        end
            
    end

    return sum_loglike
end

