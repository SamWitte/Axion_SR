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
    adaptive_alpha_sampling_smart(mu_min, mu_max, M_BH, a, n, l, m, n_probe, n_fine)

Generate adaptive alpha sampling by probing imaginary part values.

Strategy:
1. Coarse probe across full range to find where imaginary part starts significant decay
2. Refine heavily in the sharp transition/cutoff region
3. Sparse sampling in flat power-law regime

The key insight: the imaginary part is relatively flat at small alpha (power law),
then shows a sharp exponential decay. We detect this transition by computing derivatives.

Arguments:
- mu_min, mu_max: bounds of mu range
- M_BH: black hole mass
- a: spin parameter
- n, l, m: quantum numbers
- n_probe: number of probe points for transition detection (~10-15)
- n_fine: number of fine points in sharp transition region (~100+)

Returns: sorted array of alpha values with adaptive density
"""
function adaptive_alpha_sampling_smart(mu_min, mu_max, M_BH, a, n, l, m, n_probe, n_fine)
    # Phase 1: Probe with sparse sampling to find transition point
    alpha_probe = 10 .^ (range(log10(mu_min * M_BH * GNew), log10(mu_max * M_BH * GNew), n_probe))

    # Evaluate imaginary part at probe points
    imag_vals = Float64[]
    for alpha in alpha_probe
        try
            mu = alpha / (M_BH * GNew)
            erg_r, erg = find_im_part(mu, M_BH, a, n, l, m; return_both=true, for_s_rates=true, Ntot_force=5000, debug=false)
            push!(imag_vals, imag(erg))
        catch
            push!(imag_vals, 0.0)  # Handle failures gracefully
        end
    end

    # Phase 2: Detect sharp transition by computing log-derivatives
    # Look for where |d(ln|imag|)/d(ln alpha)| is largest
    max_decay_idx = 1
    if n_probe > 2
        log_alpha_probe = log10.(alpha_probe)
        max_derivative = 0.0

        for i in 2:n_probe-1
            if imag_vals[i] > 1e-20 && imag_vals[i-1] > 1e-20 && imag_vals[i+1] > 1e-20
                # Log-log derivative: d(ln y)/d(ln x)
                dlog_y = log(imag_vals[i+1]) - log(imag_vals[i-1])
                dlog_x = log_alpha_probe[i+1] - log_alpha_probe[i-1]
                if dlog_x > 0
                    derivative = abs(dlog_y / dlog_x)
                    if derivative > max_derivative
                        max_derivative = derivative
                        max_decay_idx = i
                    end
                end
            end
        end
    end

    # Phase 3: Transition point is where max decay occurs
    alpha_transition = alpha_probe[max_decay_idx]

    # Only refine if transition is not at very end (need room to resolve)
    if max_decay_idx < n_probe - 1
        alpha_fine_start = alpha_probe[max(1, max_decay_idx - 1)]
        alpha_fine_end = alpha_probe[end]
    else
        # Fallback: refine in upper 30% of range
        alpha_fine_start = alpha_probe[max(1, Int(round(n_probe * 0.7)))]
        alpha_fine_end = alpha_probe[end]
    end

    # Phase 4: Generate fine sampling only in transition region
    alpha_fine = 10 .^ (range(log10(alpha_fine_start), log10(alpha_fine_end), n_fine))

    # Phase 5: Keep probe points + fine points (sparse elsewhere, dense at transition)
    alpha_all = vcat(alpha_probe, alpha_fine)
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
        alpha_values = adaptive_alpha_sampling_smart(mu_min, mu_max, M_BH, a, n, l, m, 15, 150)
        mu_values = alpha_values ./ (M_BH * GNew)
        n_alpha = length(alpha_values)
        @printf "  Sampled %d alpha points (spacing: %.2e to %.2e)\n" n_alpha minimum(diff(alpha_values)) maximum(diff(alpha_values))

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
