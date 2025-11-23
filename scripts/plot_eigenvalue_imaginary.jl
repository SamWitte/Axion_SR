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
    adaptive_alpha_sampling_simple(mu_min, mu_max, M_BH, a, n, l, m,
                                   alpha_coarse_spacing, alpha_fine_spacing)

Generate adaptive alpha sampling with manual spacing control.

Strategy:
1. Sparse coarse sampling (e.g., 0.05 spacing in log-alpha) across full range to probe for transition
2. Fine sampling (e.g., 0.001 spacing in alpha) in detected transition region
3. Minimal computation, user controls spacing directly

Arguments:
- mu_min, mu_max: bounds of mu range
- M_BH: black hole mass
- a: spin parameter
- n, l, m: quantum numbers
- alpha_coarse_spacing: spacing for coarse points (e.g., 0.05 in log10 scale)
- alpha_fine_spacing: spacing for fine points (e.g., 0.001 in linear alpha)

Returns: sorted array of alpha values with specified density
"""
function adaptive_alpha_sampling_simple(mu_min, mu_max, M_BH, a, n, l, m,
                                        alpha_coarse_spacing=0.05, alpha_fine_spacing=0.001)
    alpha_min = mu_min * M_BH * GNew
    alpha_max = mu_max * M_BH * GNew

    # Phase 1: Generate coarse log-spaced sampling
    log_alpha_coarse = range(log10(alpha_min), log10(alpha_max), step=alpha_coarse_spacing)
    alpha_coarse = 10 .^ log_alpha_coarse

    # Phase 2: Probe coarse points to find transition
    # Transition is where imaginary part starts to decay significantly
    max_imag = 0.0
    transition_idx = 1
    for (i, alpha) in enumerate(alpha_coarse)
        try
            mu = alpha / (M_BH * GNew)
            erg_r, erg = find_im_part(mu, M_BH, a, n, l, m; return_both=true, for_s_rates=true, Ntot_force=3000, debug=false)
            imag_val = imag(erg)
            if imag_val > max_imag
                max_imag = imag_val
                transition_idx = i
            end
        catch
        end
    end

    # Phase 3: Find where decay starts (roughly 80% of max, or within transition region)
    decay_start_idx = max(1, transition_idx - 2)
    alpha_decay_start = alpha_coarse[decay_start_idx]

    # Phase 4: Generate fine linear-spaced sampling in transition region
    # Use linear spacing in alpha (not log) for fine resolution
    alpha_fine = alpha_decay_start:alpha_fine_spacing:alpha_max
    alpha_fine = collect(alpha_fine)

    # Phase 5: Combine coarse + fine, keeping only those in proper range
    alpha_all = vcat(alpha_coarse, alpha_fine)
    alpha_adaptive = sort(unique(round.(alpha_all, digits=8)))

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
mu_min = alpha_min / (M_BH * GNew)
mu_max = alpha_max / (M_BH * GNew)

@printf "[SETUP] α range: [%.3f, %.3f]\n" alpha_min alpha_max
@printf "[SETUP] μ range: [%.3e, %.3e]\n" mu_min mu_max
@printf "[SETUP] Adaptive sampling strategy: probe transitions per (n,l,m,spin)\n"

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

println("\n[COMPUTING] Computing adaptive alpha sampling and eigenvalues...")
println("-"^80)

for (n, l, m) in quantum_levels
    for a in spins
        # Compute adaptive alpha sampling for this quantum state
        @printf "\nComputing adaptive sampling for (n,l,m,a) = (%d,%d,%d,%.2f)...\n" n l m a
        alpha_values = adaptive_alpha_sampling_simple(mu_min, mu_max, M_BH, a, n, l, m, 0.05, 0.002)
        mu_values = alpha_values ./ (M_BH * GNew)
        n_alpha = length(alpha_values)
        @printf "  Sampled %d alpha points\n" n_alpha
        if n_alpha > 1
            @printf "    Spacing: min = %.2e, max = %.2e\n" minimum(diff(alpha_values)) maximum(diff(alpha_values))
        end

        for (i, mu) in enumerate(mu_values)
            if i % max(1, div(n_alpha, 5)) == 0
                @printf "  State (%d,%d,%d,%.2f): %d / %d points\n" n l m a i n_alpha
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
