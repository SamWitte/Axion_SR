using Random
using Distributions
using DelimitedFiles
using Suppressor
@suppress include("../../src/super_rad.jl")
using Statistics
using LinearAlgebra
using Plots

"""
Test script to compute p(f_a | Data, M1, q) by marginalizing over spins.
"""

# Include the multivariate KDE function from stat_1dL.jl
function multivariate_kde_pdf(data_matrix::Matrix{Float64}, eval_point::Vector{Float64}; bandwidth_factor=1.0)
    n, d = size(data_matrix)
    data_cov = cov(data_matrix)
    h = n^(-1.0/(d+4)) * bandwidth_factor
    H = h^2 * data_cov
    H = H + 1e-6 * I

    density = 0.0
    for i in 1:n
        x_i = data_matrix[i, :]
        try
            mvn = MvNormal(x_i, H)
            density += pdf(mvn, eval_point)
        catch
            continue
        end
    end
    return density / n
end

function sample_spin()
    return rand(Uniform(0.0, maxSpin))
end

function compute_fa_likelihood(log_fa, M1, q, mass_ax, kde_data; tau_max=1e4, non_rel=true, Nmax=3, cheby=true, high_p=true, N_spin_samples=20)
    """
    Compute likelihood p(f_a | M1, q, Data) by marginalizing over spins.
    """

    M2 = q * M1
    fa = 10.0^log_fa

    sum_likelihood = 0.0

    for i in 1:N_spin_samples
        # Sample spins uniformly
        s1 = sample_spin()
        s2 = sample_spin()

        # Evolve forward
        final_spin_1, final_mass_1 = super_rad_check(M1, s1, mass_ax, fa, tau_max=tau_max, non_rel=non_rel, Nmax=Nmax, cheby=cheby, high_p=high_p)
        final_spin_2, final_mass_2 = super_rad_check(M2, s2, mass_ax, fa, tau_max=tau_max, non_rel=non_rel, Nmax=Nmax, cheby=cheby, high_p=high_p)

        # Ensure final M1 >= final M2
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

        # Evaluate KDE
        eval_point = [final_M1, final_q, final_chi1, final_chi2]
        likelihood = multivariate_kde_pdf(kde_data, eval_point)

        sum_likelihood += likelihood
    end

    # Average over spin samples (marginalization)
    return sum_likelihood / N_spin_samples
end

# Main script
function run_test(dataname; mass_ax=1e-12, fa_min=1e11, fa_max=1e18, N_fa_points=20, tau_max=1e4, non_rel=true, Nmax=3, cheby=true, high_p=true, N_spin_samples=20, N_mq_samples=3)

    println("Loading data...")
    data = readdlm("../../src/BH_data/$(dataname).dat")

    # Convert data to (M1, q, chi1, chi2) format
    println("Transforming data to (M1, q) format...")
    data_transformed = copy(data)
    for i in 1:size(data_transformed, 1)
        if data_transformed[i, 2] > data_transformed[i, 1]
            # Swap masses and spins
            data_transformed[i, 1], data_transformed[i, 2] = data_transformed[i, 2], data_transformed[i, 1]
            data_transformed[i, 3], data_transformed[i, 4] = data_transformed[i, 4], data_transformed[i, 3]
        end
    end
    kde_data = copy(data_transformed)
    kde_data[:, 2] = kde_data[:, 2] ./ kde_data[:, 1]  # q = M2/M1

    # Sample a few (M1, q) points from the data to test
    println("Sampling (M1, q) points from data...")
    sample_indices = rand(1:size(kde_data, 1), N_mq_samples)

# Function to run test with specific M1 and q values
function run_test_specific(dataname, M1_specific, q_specific; mass_ax=1e-12, fa_min=1e11, fa_max=1e18, N_fa_points=20, tau_max=1e4, non_rel=true, Nmax=3, cheby=true, high_p=true, N_spin_samples=20)

    println("Loading data...")
    data = readdlm("../../src/BH_data/$(dataname).dat")

    # Convert data to (M1, q, chi1, chi2) format
    println("Transforming data to (M1, q) format...")
    data_transformed = copy(data)
    for i in 1:size(data_transformed, 1)
        if data_transformed[i, 2] > data_transformed[i, 1]
            # Swap masses and spins
            data_transformed[i, 1], data_transformed[i, 2] = data_transformed[i, 2], data_transformed[i, 1]
            data_transformed[i, 3], data_transformed[i, 4] = data_transformed[i, 4], data_transformed[i, 3]
        end
    end
    kde_data = copy(data_transformed)
    kde_data[:, 2] = kde_data[:, 2] ./ kde_data[:, 1]  # q = M2/M1

    # Create f_a grid
    log_fa_grid = LinRange(log10(fa_min), log10(fa_max), N_fa_points)

    # For each (M1, q) sample, compute p(f_a | M1, q, Data)
    for (idx, sample_idx) in enumerate(sample_indices)
        M1_test = kde_data[sample_idx, 1]
        q_test = kde_data[sample_idx, 2]

        println("\n" * "="^60)
        println("Testing point $(idx): M1 = $(M1_test), q = $(q_test)")
        println("="^60)

        likelihoods = zeros(N_fa_points)

        for (i, log_fa) in enumerate(log_fa_grid)
            println("\nEvaluating f_a = $(10.0^log_fa) ($(i)/$(N_fa_points))...")
            likelihoods[i] = compute_fa_likelihood(log_fa, M1_test, q_test, mass_ax, kde_data,
                                                   tau_max=tau_max, non_rel=non_rel, Nmax=Nmax,
                                                   cheby=cheby, high_p=high_p, N_spin_samples=N_spin_samples)
            println("  Likelihood = $(likelihoods[i])")
        end

        # Normalize to get probability
        if sum(likelihoods) > 0
            probabilities = likelihoods ./ sum(likelihoods)
        else
            probabilities = likelihoods
        end

        # Plot
        p = plot(log_fa_grid, probabilities,
                xlabel="log₁₀(fₐ / GeV)",
                ylabel="p(fₐ | M₁=$(round(M1_test, digits=2)), q=$(round(q_test, digits=3)), Data)",
                title="m_axion = $(mass_ax) eV, marginalized over spins",
                legend=false,
                linewidth=2,
                marker=:circle)

        savefig(p, "test_fa_marginalization_M1_$(round(M1_test, digits=2))_q_$(round(q_test, digits=3)).png")
        println("\nPlot saved to test_fa_marginalization_M1_$(round(M1_test, digits=2))_q_$(round(q_test, digits=3)).png")

        # Save data
        output_data = hcat(log_fa_grid, likelihoods, probabilities)
        writedlm("test_fa_marginalization_M1_$(round(M1_test, digits=2))_q_$(round(q_test, digits=3)).dat", output_data)
        println("Data saved to test_fa_marginalization_M1_$(round(M1_test, digits=2))_q_$(round(q_test, digits=3)).dat")
    end

    println("\n" * "="^60)
    println("Test complete!")
    println("="^60)
end

# Run the test
# Example: run_test("LIGO_test", N_fa_points=10, N_spin_samples=10, N_mq_samples=2)
println("Script loaded. Run with:")
println("  run_test(\"LIGO_test\", N_fa_points=10, N_spin_samples=10, N_mq_samples=2)")
println("\nParameters:")
println("  dataname: name of data file (without .dat)")
println("  mass_ax: axion mass in eV (default: 1e-12)")
println("  fa_min, fa_max: f_a range (default: 1e11 to 1e18)")
println("  N_fa_points: number of f_a points to evaluate")
println("  N_spin_samples: number of spin samples for marginalization")
println("  N_mq_samples: number of (M1,q) points to test")
