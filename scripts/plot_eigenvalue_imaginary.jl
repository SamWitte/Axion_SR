#!/usr/bin/env julia
"""
Compute and Save Eigenvalue Data for Full Range

For a 1 solar mass black hole at various spins, computes the imaginary components
of all eigenvalues and saves to HDF5 for plotting with matplotlib.

The dimensionless parameter α = μ * M * G_N is constrained to [0.03, 1].
This script computes the appropriate μ range based on this constraint.

Uses adaptive alpha sampling to finely resolve exponential cutoff while avoiding
excessive calculations in the power-law regime at small alpha.
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

"""
    adaptive_alpha_sampling(alpha_min, alpha_max, n_coarse, n_fine;
                           quantile_threshold=0.8, coarse_factor=0.5)

Generate adaptive alpha sampling that focuses resolution on the exponential cutoff region.

Strategy:
1. Sample coarsely across full range to identify transition region
2. Refine sampling in the region where function starts to decay (transition to exponential)
3. Return merged sorted samples

Arguments:
- alpha_min, alpha_max: bounds of alpha range
- n_coarse: number of coarse samples (initial sweep)
- n_fine: number of fine samples (refined region)
- quantile_threshold: trigger refinement when derivative crosses this quantile of max derivative
- coarse_factor: fraction of alpha range to refine (default 0.5 means refine upper 50%)

Returns: sorted array of alpha values with adaptive density
"""
function adaptive_alpha_sampling(alpha_min, alpha_max, n_coarse, n_fine;
                                 quantile_threshold=0.75, coarse_factor=0.5)
    # Phase 1: Coarse log-spaced sampling across full range
    alpha_coarse = 10 .^ (range(log10(alpha_min), log10(alpha_max), n_coarse))

    # Phase 2: Identify transition region by looking at log-derivatives
    # For power law: d(ln f)/d(ln α) ≈ constant
    # For exponential: d(ln f)/d(ln α) ≈ α * f'/f → large
    # We estimate this with finite differences
    if n_coarse > 3
        log_alpha = log10.(alpha_coarse)
        # Compute finite differences to estimate where curvature increases
        # This is a proxy for detecting the transition region
        d_indices = 2:n_coarse-1
        curvatures = similar(alpha_coarse)
        fill!(curvatures, 0.0)

        # Use log-spacing to estimate second derivative in log-log space
        for i in d_indices
            if log_alpha[i+1] - log_alpha[i] > 0 && log_alpha[i] - log_alpha[i-1] > 0
                # Finite difference of log-spacing (indicator of curvature)
                curvatures[i] = abs((log_alpha[i+1] - log_alpha[i]) - (log_alpha[i] - log_alpha[i-1]))
            end
        end

        # Find transition point: where curvature is highest (or use quantile approach)
        # Transition typically occurs in upper part of range
        transition_idx = max(Int(round(n_coarse * (1 - coarse_factor))), 2)
        alpha_transition = alpha_coarse[transition_idx]
    else
        # If very few coarse points, refine upper half
        alpha_transition = 10^((log10(alpha_min) + log10(alpha_max)) / 2)
    end

    # Phase 3: Generate fine-grained sampling in transition/exponential region
    # Use higher density near alpha_max to resolve cutoff
    alpha_fine_start = alpha_transition
    alpha_fine = 10 .^ (range(log10(alpha_fine_start), log10(alpha_max), n_fine))

    # Phase 4: Merge and remove duplicates
    alpha_all = vcat(alpha_coarse, alpha_fine)
    alpha_adaptive = sort(unique(round.(alpha_all, digits=8)))  # Remove near-duplicates

    return alpha_adaptive
end

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

# Use adaptive alpha sampling to focus resolution on exponential cutoff
# n_coarse: initial coarse sweep (30 points across full range)
# n_fine: refined points in transition region (100+ points)
# This gives ~120-130 total unique points with much higher density near cutoff
alpha_values = adaptive_alpha_sampling(alpha_min, alpha_max, 30, 100; coarse_factor=0.55)
mu_values = alpha_values ./ (M_BH * GNew)

@printf "[SETUP] α range: [%.3f, %.3f]\n" alpha_min alpha_max
@printf "[SETUP] μ range: [%.3e, %.3e]\n" mu_values[1] mu_values[end]
@printf "[SETUP] Boson mass points: %d (adaptive sampling)\n" length(alpha_values)
@printf "[SETUP] Point spacing: min Δα = %.2e, max Δα = %.2e\n" minimum(diff(alpha_values)) maximum(diff(alpha_values))

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
