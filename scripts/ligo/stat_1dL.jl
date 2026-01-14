using Random
using Distributions
import AffineInvariantMCMC
using Turing
using StatsPlots
using DelimitedFiles
using Suppressor
@suppress include("../../src/super_rad.jl")
using MCMCDiagnosticTools
using Dates
using KernelDensity
using LinearAlgebra
using Statistics

"""
    multivariate_kde_pdf(data_matrix, eval_point; bandwidth_factor=1.0)

Compute the probability density at eval_point using a multivariate Gaussian KDE.

# Arguments
- `data_matrix`: N×D matrix where N is number of samples, D is dimensionality
- `eval_point`: D-dimensional vector at which to evaluate the density
- `bandwidth_factor`: Multiplier for the bandwidth (default 1.0 uses Scott's rule)

# Returns
- Probability density value at eval_point
"""
function multivariate_kde_pdf(data_matrix::Matrix{Float64}, eval_point::Vector{Float64}; bandwidth_factor=1.0)
    n, d = size(data_matrix)

    # Scott's rule: h = n^(-1/(d+4)) * σ
    # For multivariate, we use the covariance matrix
    data_cov = cov(data_matrix)

    # Compute bandwidth using Scott's rule
    h = n^(-1.0/(d+4)) * bandwidth_factor

    # Bandwidth covariance matrix
    H = h^2 * data_cov

    # Add small regularization to ensure positive definiteness
    H = H + 1e-6 * I

    # Evaluate KDE: sum of Gaussian kernels centered at each data point
    density = 0.0
    for i in 1:n
        x_i = data_matrix[i, :]
        diff = eval_point - x_i

        # Multivariate normal PDF: (2π)^(-d/2) |H|^(-1/2) exp(-0.5 * diff' * H^(-1) * diff)
        try
            mvn = MvNormal(x_i, H)
            density += pdf(mvn, eval_point)
        catch
            # If H is still singular, skip this kernel
            continue
        end
    end

    # Normalize by number of data points
    return density / n
end

function profileL_func_minimize(data, mass_ax, Fname, Nsamples; fa_min=1e11, fa_max=1e18, tau_max=1e4, non_rel=true, Nmax=3, cheby=true, delt_M=0.05, thinning=1, numsamples_perwalker=2000, burnin=500, use_kde=true, over_run=true, high_p=true, one_BH=false, high_spin_cut=nothing)
    
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
    
    fout = "../../src/output_mcmc/"*Fname*"_mcmc.dat"
    if isfile(fout) && over_run
        fin = readdlm(fout)
        fnew = cat(fin, samples, dims=1)
        writedlm(fout, fnew)
    else
        writedlm(fout, samples)
    end
    
end
    
    

function log_probability(theta, data, mass_ax; fa_min=1e11, fa_max=1e20, tau_max=5e7, non_rel=false, Nmax=3, cheby=true, delt_M=0.05, Nsamples=3, use_kde=true, high_p=true, one_BH=false, high_spin_cut=nothing)
    
    return log_likelihood(theta, data, mass_ax, tau_max=tau_max, non_rel=non_rel, Nmax=Nmax, cheby=cheby, delt_M=delt_M, Nsamples=Nsamples, use_kde=use_kde, high_p=high_p, one_BH=one_BH, high_spin_cut=high_spin_cut)
end

function sample_spin()
    # Sample spin uniformly between 0 and maxSpin (0.998)
    SpinBH = rand(Uniform(0.0, maxSpin))
    return SpinBH
end

function log_likelihood(theta, data, mass_ax; tau_max=1e4, non_rel=true, debug=false, Nmax=3, cheby=true, Nsamples=3, delt_M=0.05, use_kde=true, high_p=true, one_BH=false, high_spin_cut=nothing)

    log_f = theta
    sum_loglike = 0.0

    # Convert data to (M1, q, chi1, chi2) format where M1 >= M2 and q = M2/M1
    # Ensure M1 is the larger mass by swapping if needed
    data_transformed = copy(data)
    for i in 1:size(data_transformed, 1)
        if data_transformed[i, 2] > data_transformed[i, 1]
            # Swap masses and spins
            data_transformed[i, 1], data_transformed[i, 2] = data_transformed[i, 2], data_transformed[i, 1]
            data_transformed[i, 3], data_transformed[i, 4] = data_transformed[i, 4], data_transformed[i, 3]
        end
    end
    # Now convert M2 to q = M2/M1
    data_q = copy(data_transformed)
    data_q[:, 2] = data_q[:, 2] ./ data_q[:, 1]  # q = M2/M1

    # Now sample from the transformed data
    maxI = size(data_q)[1]
    idx_hold = rand(1:maxI, Nsamples)
    sampled_data = data_q[idx_hold, :]

    ## note to self: check if strong correlations in mass/spin, could be something to think about putting Mf on cut
    for i in 1:Nsamples
        # Now sampled_data is in (M1, q, chi1, chi2) format
        M1 = sampled_data[i, 1]
        q_sample = sampled_data[i, 2]
        M2 = q_sample * M1

        delt_M_use = delt_M
        kde_data = []
        N_min_kde = 50
        if length(data_q[:,1]) <= N_min_kde
            N_min_kde = length(data_q[:,1]) ./ 2
        end

        # Filter data based on M1 and q
        while length(kde_data) < N_min_kde
            kde_data = data_q[(data_q[:, 1] .>= M1 .* (1.0 .- delt_M_use)) .&
                              (data_q[:, 1] .<= M1 .* (1.0 .+ delt_M_use)) .&
                              (data_q[:, 2] .>= q_sample .* (1.0 .- delt_M_use)) .&
                              (data_q[:, 2] .<= q_sample .* (1.0 .+ delt_M_use)), :]
            delt_M_use *= 1.1
        end

        if isnothing(high_spin_cut)
            s1 = sample_spin()
            s2 = sample_spin()
        else
            s1 = 1.0
            s2 = 1.0
            while s1 > high_spin_cut
                s1 = sample_spin()
            end
            while s2 > high_spin_cut
                s2 = sample_spin()
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
            final_mass_2 = M2
        end

        # Ensure final M1 >= final M2, swap if needed for q calculation
        if final_mass_2 > final_mass_1
            final_M1 = final_mass_2
            final_M2 = final_mass_1
            final_chi1 = final_spin_2
            final_chi2 = final_spin_1
        else
            final_M1 = final_mass_1
            final_M2 = final_mass_2
            final_chi1 = final_spin_1
            final_chi2 = final_spin_2
        end
        final_q = final_M2 / final_M1

        if use_kde
            # Perform 4D multivariate kernel density estimation on (M1, q, chi1, chi2)
            eval_point = [final_M1, final_q, final_chi1, final_chi2]
            likelihood = multivariate_kde_pdf(kde_data, eval_point)
        else
            # Fallback to simple Gaussian if not using KDE
            M1_mean = mean(kde_data[:, 1])
            q_mean = mean(kde_data[:, 2])
            chi1_mean = mean(kde_data[:, 3])
            chi2_mean = mean(kde_data[:, 4])
            M1_std = std(kde_data[:, 1])
            q_std = std(kde_data[:, 2])
            chi1_std = std(kde_data[:, 3])
            chi2_std = std(kde_data[:, 4])

            likelihood = exp.(-(final_M1 .- M1_mean).^2 ./ (2 .* M1_std.^2))
            likelihood .*= exp.(-(final_q .- q_mean).^2 ./ (2 .* q_std.^2))
            likelihood .*= exp.(-(final_chi1 .- chi1_mean).^2 ./ (2 .* chi1_std.^2))
            if !one_BH
                likelihood .*= exp.(-(final_chi2 .- chi2_mean).^2 ./ (2 .* chi2_std.^2))
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
    
    fout = "../../src/output_mcmc/"*Fname*"_mcmc.dat"
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

    # Convert data to (M1, q, chi1, chi2) format where M1 >= M2 and q = M2/M1
    # Ensure M1 is the larger mass by swapping if needed
    data_transformed = copy(data)
    for i in 1:size(data_transformed, 1)
        if data_transformed[i, 2] > data_transformed[i, 1]
            # Swap masses and spins
            data_transformed[i, 1], data_transformed[i, 2] = data_transformed[i, 2], data_transformed[i, 1]
            data_transformed[i, 3], data_transformed[i, 4] = data_transformed[i, 4], data_transformed[i, 3]
        end
    end
    # Now convert M2 to q = M2/M1
    data_q = copy(data_transformed)
    data_q[:, 2] = data_q[:, 2] ./ data_q[:, 1]  # q = M2/M1

    # Now sample from the transformed data
    maxI = size(data_q)[1]
    idx_hold = rand(1:maxI, Nsamples)
    sampled_data = data_q[idx_hold, :]

    ## note to self: check if strong correlations in mass/spin, could be something to think about putting Mf on cut
    for i in 1:Nsamples
        # Now sampled_data is in (M1, q, chi1, chi2) format
        M1 = sampled_data[i, 1]
        q_sample = sampled_data[i, 2]
        M2 = q_sample * M1

        delt_M_use = delt_M
        kde_data = []
        N_min_kde = 50
        if length(data_q[:,1]) <= N_min_kde
            N_min_kde = length(data_q[:,1]) ./ 2
        end

        # Filter data based on M1 and q
        while length(kde_data) < N_min_kde
            kde_data = data_q[(data_q[:, 1] .>= M1 .* (1.0 .- delt_M_use)) .&
                              (data_q[:, 1] .<= M1 .* (1.0 .+ delt_M_use)) .&
                              (data_q[:, 2] .>= q_sample .* (1.0 .- delt_M_use)) .&
                              (data_q[:, 2] .<= q_sample .* (1.0 .+ delt_M_use)), :]
            delt_M_use *= 1.1
        end

        s1 = sample_spin()
        s2 = sample_spin()

        alph1 = GNew .* M1 .* mass_ax #
        alph2 = GNew .* M2 .* mass_ax #


        final_spin_1, final_mass_1 = super_rad_check(M1, s1, mass_ax, 1.0, spinone=true, tau_max=tau_max)
        final_spin_2, final_mass_2 = super_rad_check(M2, s2, mass_ax, 1.0, spinone=true, tau_max=tau_max)

        println("Init / Final BH 1: ", M1, " ", s1, " ", final_mass_1, " ", final_spin_1)
        println("Init / Final BH 2: ", M2, " ", s2, " ", final_mass_2, " ", final_spin_2)

        # Ensure final M1 >= final M2, swap if needed for q calculation
        if final_mass_2 > final_mass_1
            final_M1 = final_mass_2
            final_M2 = final_mass_1
            final_chi1 = final_spin_2
            final_chi2 = final_spin_1
        else
            final_M1 = final_mass_1
            final_M2 = final_mass_2
            final_chi1 = final_spin_1
            final_chi2 = final_spin_2
        end
        final_q = final_M2 / final_M1

        if use_kde
            # Perform 4D multivariate kernel density estimation on (M1, q, chi1, chi2)
            eval_point = [final_M1, final_q, final_chi1, final_chi2]
            likelihood = multivariate_kde_pdf(kde_data, eval_point)
        else
            # Fallback to simple Gaussian if not using KDE
            M1_mean = mean(kde_data[:, 1])
            q_mean = mean(kde_data[:, 2])
            chi1_mean = mean(kde_data[:, 3])
            chi2_mean = mean(kde_data[:, 4])
            M1_std = std(kde_data[:, 1])
            q_std = std(kde_data[:, 2])
            chi1_std = std(kde_data[:, 3])
            chi2_std = std(kde_data[:, 4])

            likelihood = exp.(-(final_M1 .- M1_mean).^2 ./ (2 .* M1_std.^2))
            likelihood .*= exp.(-(final_q .- q_mean).^2 ./ (2 .* q_std.^2))
            likelihood .*= exp.(-(final_chi1 .- chi1_mean).^2 ./ (2 .* chi1_std.^2))
            likelihood .*= exp.(-(final_chi2 .- chi2_mean).^2 ./ (2 .* chi2_std.^2))
        end
        if likelihood .> 0.0
            sum_loglike += log.(likelihood)
        else
            sum_loglike += -Inf
        end

    end

    return sum_loglike
end

