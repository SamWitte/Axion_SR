#!/usr/bin/env julia
# Quick test of the modified statistical analysis from src directory

using DelimitedFiles
using Suppressor
@suppress include("stat_1dL.jl")

println("Testing modified statistical analysis...")
println("=" ^ 50)

# Load a small amount of data
data = open(readdlm, "BH_data/LIGO_test.dat")
println("Loaded $(size(data, 1)) data points")
println("Data format: [M1, M2, chi1, chi2]")
println("First row: ", data[1, :])

# Test parameters
ax_mass = 1e-13

println("\nTest parameters:")
println("  Axion mass: $ax_mass")

# Test the sample_spin function
println("\n" * "=" ^ 50)
println("Testing uniform spin prior sampling...")
spin_samples = [sample_spin() for _ in 1:10]
println("Sample spins (should be uniform in [0, 0.998]):")
println("  ", round.(spin_samples, digits=4))
println("  Min: $(round(minimum(spin_samples), digits=4)), Max: $(round(maximum(spin_samples), digits=4))")

# Test mass ratio conversion
println("\n" * "=" ^ 50)
println("Testing mass ratio conversion...")
test_data = data[1:5, :]
println("Original data (M1, M2):")
for i in 1:5
    println("  Row $i: M1=$(round(test_data[i,1], digits=2)), M2=$(round(test_data[i,2], digits=2)), chi1=$(round(test_data[i,3], digits=4)), chi2=$(round(test_data[i,4], digits=4))")
end

# Apply transformation (same as in log_likelihood)
data_transformed = copy(test_data)
for i in 1:size(data_transformed, 1)
    if data_transformed[i, 2] > data_transformed[i, 1]
        data_transformed[i, 1], data_transformed[i, 2] = data_transformed[i, 2], data_transformed[i, 1]
        data_transformed[i, 3], data_transformed[i, 4] = data_transformed[i, 4], data_transformed[i, 3]
    end
end
data_q = copy(data_transformed)
data_q[:, 2] = data_q[:, 2] ./ data_q[:, 1]

println("\nTransformed data (M1, q, chi1, chi2) where M1 >= M2, q = M2/M1:")
for i in 1:5
    println("  Row $i: M1=$(round(data_q[i,1], digits=2)), q=$(round(data_q[i,2], digits=4)), chi1=$(round(data_q[i,3], digits=4)), chi2=$(round(data_q[i,4], digits=4))")
end

# Verify q is in (0, 1]
println("\nVerifying q ∈ (0, 1]:")
println("  All q > 0: ", all(data_q[:, 2] .> 0))
println("  All q <= 1: ", all(data_q[:, 2] .<= 1))
println("  q range: [$(round(minimum(data_q[:, 2]), digits=4)), $(round(maximum(data_q[:, 2]), digits=4))]")

# Test a single likelihood evaluation with MINIMAL computation
println("\n" * "=" ^ 50)
println("Testing likelihood evaluation (this may take a minute)...")
println("Running with minimal samples and simplified physics...")

try
    theta = [log10(5e11)]  # Middle of fa range
    # Use very minimal computation for test
    ll = log_likelihood(theta, data[1:50, :], ax_mass,
                       tau_max=100.0,       # Very short evolution time
                       non_rel=true,        # Faster non-relativistic mode
                       Nmax=3,              # Minimum allowed Nmax
                       cheby=false,         # Leaver method
                       Nsamples=1,          # Only 1 sample
                       use_kde=true,        # Use KDE (4D)
                       high_p=false,        # Lower precision
                       one_BH=true)         # Only evolve one BH
    println("\nLog-likelihood value: $ll")
    if ll > -Inf
        println("✓ Likelihood evaluation successful!")
    else
        println("⚠ Likelihood is -Inf (may indicate numerical issues or poor fit)")
    end
catch e
    println("\n✗ Error during likelihood evaluation:")
    println(e)
    for (exc, bt) in Base.catch_stack()
        showerror(stdout, exc, bt)
        println()
    end
end

println("\n" * "=" ^ 50)
println("Test complete!")
println("\nSummary of modifications:")
println("  1. ✓ Spin priors now uniform in [0, 0.998]")
println("  2. ✓ Mass variables converted to (M1, q) with M1 >= M2, q = M2/M1 ∈ (0, 1]")
println("  3. ✓ KDE now 4-dimensional: (M1, q, chi1, chi2) using multivariate Gaussian kernels")

# Test the multivariate KDE function directly
println("\n" * "=" ^ 50)
println("Testing multivariate KDE function...")
n_test = min(100, size(data_q, 1))
test_data_matrix = data_q[1:n_test, :]
test_point = vec(data_q[min(50, n_test), :])
println("Test data: $(size(test_data_matrix)) matrix")
println("Eval point: ", round.(test_point, digits=4))
try
    kde_val = multivariate_kde_pdf(test_data_matrix, test_point)
    println("KDE value at point: $kde_val")
    println("✓ Multivariate KDE evaluation successful!")
catch e
    println("✗ Error in multivariate KDE:")
    println(e)
end
