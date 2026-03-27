#!/usr/bin/env julia
# Quick test of the modified statistical analysis

using DelimitedFiles
using Suppressor
@suppress include("stat_1dL.jl")

println("Testing modified statistical analysis...")
println("=" ^ 50)

# Load a small amount of data
data = open(readdlm, "../../src/BH_data/LIGO_test.dat")
println("Loaded $(size(data, 1)) data points")
println("Data format: [M1, M2, chi1, chi2]")
println("First row: ", data[1, :])

# Test parameters
ax_mass = 1e-13
fa_min = 1e11
fa_max = 1e12

println("\nTest parameters:")
println("  Axion mass: $ax_mass")
println("  fa range: [$fa_min, $fa_max]")

# Test the sample_spin function
println("\n" * "=" ^ 50)
println("Testing uniform spin prior sampling...")
spin_samples = [sample_spin() for _ in 1:10]
println("Sample spins (should be uniform in [0, 0.998]):")
println("  ", spin_samples)
println("  Min: $(minimum(spin_samples)), Max: $(maximum(spin_samples))")

# Test mass ratio conversion
println("\n" * "=" ^ 50)
println("Testing mass ratio conversion...")
test_data = data[1:5, :]
println("Original data (M1, M2):")
for i in 1:5
    println("  Row $i: M1=$(test_data[i,1]), M2=$(test_data[i,2]), chi1=$(test_data[i,3]), chi2=$(test_data[i,4])")
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

println("\nTransformed data (M1, q, chi1, chi2):")
for i in 1:5
    println("  Row $i: M1=$(data_q[i,1]), q=$(data_q[i,2]), chi1=$(data_q[i,3]), chi2=$(data_q[i,4])")
end

# Verify q is in (0, 1]
println("\nVerifying q ∈ (0, 1]:")
println("  All q > 0: ", all(data_q[:, 2] .> 0))
println("  All q <= 1: ", all(data_q[:, 2] .<= 1))

# Test a single likelihood evaluation
println("\n" * "=" ^ 50)
println("Testing likelihood evaluation...")
println("This will run a single likelihood calculation with minimal samples...")

try
    theta = [log10(5e11)]  # Middle of fa range
    ll = log_likelihood(theta, data[1:100, :], ax_mass,
                       tau_max=1e3, non_rel=true, Nmax=2, cheby=false,
                       Nsamples=1, use_kde=true, high_p=false, one_BH=true)
    println("\nLog-likelihood value: $ll")
    println("Test completed successfully!")
catch e
    println("\nError during likelihood evaluation:")
    println(e)
    for (exc, bt) in Base.catch_stack()
        showerror(stdout, exc, bt)
        println()
    end
end

println("\n" * "=" ^ 50)
println("Test complete!")
