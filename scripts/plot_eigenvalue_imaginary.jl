#!/usr/bin/env julia
"""
Compute and Save Eigenvalue Data for Full Range

For a 1 solar mass black hole at various spins, computes the imaginary components
of all eigenvalues and saves to HDF5 for plotting with matplotlib.

The dimensionless parameter α = μ * M * G_N is constrained to [0.03, 1].
This script computes the appropriate μ range based on this constraint.
"""

using HDF5
using DelimitedFiles
using Printf
using Statistics

# Include the necessary modules
src_dir = joinpath(@__DIR__, "..", "src")
include(joinpath(src_dir, "Core/constants.jl"))
include(joinpath(src_dir, "heunc.jl"))
include(joinpath(src_dir, "solve_sr_rates.jl"))

println("="^80)
println("EIGENVALUE DATA GENERATION (FULL RANGE)")
println("System: 1 solar mass black hole at various spins")
println("="^80)

# Physical parameters
M_BH = 1.0  # Solar masses (all calculations in these units)
@printf "\n[SETUP] Black hole mass: M = %.1f M☉\n" M_BH

# Spin parameters to explore
spins = [0.0, 0.5, 0.7, 0.9, 0.95, 0.99]
@printf "[SETUP] Spin parameters: %s\n" join(spins, ", ")

# Physical constraint: 0.03 ≤ α = μ * M * G_N ≤ 1
alpha_min = 0.03
alpha_max = 1.0

# Calculate mu range from alpha range
mu_min = alpha_min / (M_BH * GNew)
mu_max = alpha_max / (M_BH * GNew)

n_mu_points = 50
mu_values = 10 .^ (range(log10(mu_min), log10(mu_max), n_mu_points))
alpha_values = mu_values .* M_BH .* GNew

@printf "[SETUP] α range: [%.3f, %.3f]\n" alpha_min alpha_max
@printf "[SETUP] μ range: [%.3e, %.3e]\n" mu_min mu_max
@printf "[SETUP] Boson mass points: %d (log-spaced)\n" n_mu_points

# Quantum levels to consider: (n, l, m)
quantum_levels = [
    (2, 1, 1),
    (3, 1, 1),
    (3, 2, 1),
    (3, 2, 2),
    (4, 1, 1),
    (4, 3, 1),
    (4, 3, 3),
]
@printf "[SETUP] Quantum levels: %d\n" length(quantum_levels)

# Output directory for data
data_dir = joinpath(@__DIR__, "..", "data")
isdir(data_dir) || mkdir(data_dir)

# Data storage
eigenvalue_data = Dict()

total = length(quantum_levels) * length(spins) * length(mu_values)
println("\n[COMPUTING] Calculating eigenvalues ($total total computations)...")
println("-"^80)

comp_count = Ref(0)
for (n, l, m) in quantum_levels
    for a in spins
        for (i, mu) in enumerate(mu_values)
            comp_count[] += 1
            if comp_count[] % max(1, div(total, 20)) == 0
                @printf "  Progress: %d / %d (%.1f%%)\n" comp_count[] total (100 * comp_count[] / total)
            end

            try
                erg_r, erg = find_im_part(mu, M_BH, a, n, l, m; return_both=true, for_s_rates=true, Ntot_force=5000, debug=false)

                if !haskey(eigenvalue_data, (n, l, m, a))
                    eigenvalue_data[(n, l, m, a)] = (Float64[], Float64[], Float64[])
                end

                mu_arr, alpha_arr, imag_arr = eigenvalue_data[(n, l, m, a)]
                push!(mu_arr, mu)
                push!(alpha_arr, mu * M_BH * GNew)
                push!(imag_arr, imag(erg))

            catch e
                @warn "Failed: (n=$n, l=$l, m=$m, a=$a, μ=$mu)"
            end
        end
    end
end

println("\n[SAVING] Organizing and saving computed data...")

# Organize data by quantum level and save to HDF5
data_by_level = Dict()
h5_file = joinpath(data_dir, "eigenvalues_full.h5")

h5open(h5_file, "w") do file
    for ((n, l, m, a), (mus, alphas, imag_parts)) in eigenvalue_data
        if !haskey(data_by_level, (n, l, m))
            data_by_level[(n, l, m)] = Dict()
        end

        # Sort by mu
        sorted_indices = sortperm(mus)
        mus_sorted = mus[sorted_indices]
        alphas_sorted = alphas[sorted_indices]
        imag_sorted = imag_parts[sorted_indices]

        data_by_level[(n, l, m)][a] = (mus_sorted, alphas_sorted, imag_sorted)

        # Save to HDF5
        group_name = @sprintf "level_%d_%d_%d/spin_%.2f" n l m a
        file[group_name * "/mu"] = mus_sorted
        file[group_name * "/alpha"] = alphas_sorted
        file[group_name * "/imag_eigenvalue"] = imag_sorted
        file[group_name * "/M_BH"] = M_BH
        file[group_name * "/spin"] = a
    end
end

println("[SAVED] Data written to: $h5_file")

# Save a summary CSV for quick reference
summary_file = joinpath(data_dir, "eigenvalues_full_summary.csv")
open(summary_file, "w") do f
    write(f, "n,l,m,spin,n_points,im_min,im_max,im_mean\n")
    for (n, l, m) in quantum_levels
        if haskey(data_by_level, (n, l, m))
            for a in spins
                if haskey(data_by_level[(n, l, m)], a)
                    mus, alphas, imag_parts = data_by_level[(n, l, m)][a]
                    if length(imag_parts) > 0
                        @printf f "%d,%d,%d,%.2f,%d,%e,%e,%e\n" n l m a length(imag_parts) minimum(imag_parts) maximum(imag_parts) mean(imag_parts)
                    end
                end
            end
        end
    end
end

println("[SAVED] Summary written to: $summary_file")

println("\n" * "="^80)
println("DATA GENERATION COMPLETE")
println("="^80)
println("Output files:")
println("  - HDF5 data: $h5_file")
println("  - CSV summary: $summary_file")
println("\nNext step: Use plot_eigenvalues.py to create matplotlib plots")
println("  python scripts/plot_eigenvalues.py --input $h5_file --output plts/")
