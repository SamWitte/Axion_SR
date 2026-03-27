#!/usr/bin/env julia
"""
Quick Test: Compute and Save Eigenvalue Data

Lighter version for testing - uses fewer quantum levels and boson mass points.
Computes imaginary eigenvalue components and saves to HDF5/CSV for plotting with matplotlib.

The dimensionless parameter α = μ * M * G_N is constrained to [0.03, 1].
This script computes the appropriate μ range based on this constraint.

Uses adaptive alpha sampling to finely resolve exponential cutoff while avoiding
excessive calculations in the power-law regime at small alpha.
"""

using HDF5
using Printf
using Statistics

# Include necessary modules
src_dir = joinpath(@__DIR__, "..", "src")
include(joinpath(src_dir, "Core/constants.jl"))
include(joinpath(src_dir, "heunc.jl"))
include(joinpath(src_dir, "solve_sr_rates.jl"))

"""
    adaptive_alpha_sampling_simple(mu_min, mu_max, M_BH, a, n, l, m,
                                   alpha_coarse_spacing, alpha_fine_spacing)

Generate adaptive alpha sampling by detecting zero crossings in derivatives.

Strategy:
1. Coarse run: sparse log-spaced sampling across full range to probe for zero-crossing point
2. Fine run: linear-spaced grid from 80% of zero-crossing alpha to larger values
3. Stop fine run when eigenvalue becomes negative

Arguments:
- mu_min, mu_max: bounds of mu range
- M_BH: black hole mass
- a: spin parameter
- n, l, m: quantum numbers
- alpha_coarse_spacing: spacing for coarse points (e.g., 0.05 in log10 scale)
- alpha_fine_spacing: spacing for fine points (e.g., 0.001 in linear alpha)

Returns: tuple (alpha_values, imag_values) - sorted alpha array and corresponding imaginary parts
"""
function adaptive_alpha_sampling_simple(mu_min, mu_max, M_BH, a, n, l, m,
                                        alpha_coarse_spacing=0.1, alpha_fine_spacing=0.002)
    alpha_min = mu_min * M_BH * GNew
    alpha_max = mu_max * M_BH * GNew

    # Phase 1: Generate coarse log-spaced sampling
    log_alpha_coarse = range(log10(alpha_min), log10(alpha_max), step=alpha_coarse_spacing)
    alpha_coarse = 10 .^ log_alpha_coarse

    # Phase 2: Probe coarse points to find where derivative of imaginary part is zero
    # This identifies the transition point from growth to decay
    coarse_imag = Float64[]
    alpha_erg_pairs = Tuple{Float64, Complex}[]  # Store (alpha, erg) pairs to reuse
    for alpha in alpha_coarse
        try
            mu = alpha / (M_BH * GNew)
            erg_r, erg = find_im_part(mu, M_BH, a, n, l, m; return_both=true, for_s_rates=true, Ntot_force=40000, debug=false, xtol=1e-90, prec=400)
            push!(coarse_imag, erg)
            push!(alpha_erg_pairs, (alpha, erg))
        catch
            push!(coarse_imag, NaN)
        end
    end

    # Phase 3: Find first point where derivative changes sign (zero crossing)
    zero_crossing_idx = 1
    for i in 2:length(coarse_imag)-1
        if !isnan(coarse_imag[i-1]) && !isnan(coarse_imag[i]) && !isnan(coarse_imag[i+1])
            deriv_before = coarse_imag[i] - coarse_imag[i-1]
            deriv_after = coarse_imag[i+1] - coarse_imag[i]
            # Check for zero crossing: derivative changes from positive to negative
            if deriv_before > 0 && deriv_after < 0
                println("xing \t", alpha_coarse[i])
                zero_crossing_idx = i
                break
            end
        end
    end

    # Phase 4: Define fine grid starting at 80% of the zero-crossing alpha value
    alpha_zero_crossing = alpha_coarse[zero_crossing_idx]
    alpha_fine_start = 0.8 * alpha_zero_crossing
    println("alpha fine start \t", alpha_fine_start)
    
    # Phase 5: Generate fine linear-spaced sampling, stopping when eigenvalue becomes negative
    alpha_erg_fine_pairs = Tuple{Float64, Complex}[]
    alpha = alpha_fine_start
    while alpha <= alpha_max
        try
            mu = alpha / (M_BH * GNew)
            erg_r, erg = find_im_part(mu, M_BH, a, n, l, m; return_both=true, for_s_rates=true, Ntot_force=40000, debug=false, xtol=1e-90, prec=400)
            imag_part = imag(erg)

            # Stop if eigenvalue becomes negative
            if imag_part < 0
                break
            end
            push!(alpha_erg_fine_pairs, (alpha, erg))
        catch
            # If computation fails, try to continue
        end
        alpha += alpha_fine_spacing
    end

    # Phase 6: Combine coarse + fine pairs, deduplicate, and sort
    all_pairs = vcat(alpha_erg_pairs, alpha_erg_fine_pairs)
    # Deduplicate by alpha values rounded to 8 decimals
    unique_pairs = Dict{Float64, Complex}()
    for (alpha, erg) in all_pairs
        alpha_key = round(alpha, digits=8)
        if !haskey(unique_pairs, alpha_key)
            unique_pairs[alpha_key] = erg
        end
    end

    # Sort by alpha
    sorted_alphas = sort(collect(keys(unique_pairs)))
    sorted_ergs = [unique_pairs[alpha] for alpha in sorted_alphas]

    return sorted_alphas, sorted_ergs
end

println("="^80)
println("QUICK EIGENVALUE DATA GENERATION TEST")
println("="^80)

M_BH = 1.0  # Solar masses
spins = [0.5, 0.9]  # Fewer spins for quick test

# Physical constraint: 0.03 ≤ α = μ * M * G_N ≤ 1
alpha_min = 0.05
alpha_max = 1.5
mu_min = alpha_min / (M_BH * GNew)
mu_max = alpha_max / (M_BH * GNew)

@printf "[SETUP] M_BH = %.1f M☉\n" M_BH
@printf "[SETUP] α range: [%.3f, %.3f]\n" alpha_min alpha_max
@printf "[SETUP] μ range: [%.3e, %.3e]\n" mu_min mu_max
@printf "[SETUP] Adaptive sampling strategy: probe transitions per (n,l,m,spin)\n"
@printf "[SETUP] Spins: %s\n" join(spins, ", ")

# Fewer quantum levels for quick test
quantum_levels = [
    (3, 2, 2),
]
@printf "[SETUP] Quantum levels: %d\n" length(quantum_levels)

data_dir = joinpath(@__DIR__, "..", "data")
isdir(data_dir) || mkdir(data_dir)

eigenvalue_data = Dict()

println("\nComputing adaptive alpha sampling and eigenvalues...")

eval_count = Ref(0)
for (n, l, m) in quantum_levels
    for a in spins
        # Compute adaptive alpha sampling for this quantum state
        @printf "\nComputing adaptive sampling for (n,l,m,a) = (%d,%d,%d,%.2f)...\n" n l m a
        alpha_values, erg_values = adaptive_alpha_sampling_simple(mu_min, mu_max, M_BH, a, n, l, m, 0.15, 0.005)
        println("Alpha values: ", alpha_values)
        print("\n\n")
        mu_values = alpha_values ./ (M_BH * GNew)
        n_alpha = length(alpha_values)
        @printf "  Sampled %d alpha points\n" n_alpha
        if n_alpha > 1
            @printf "    Spacing: min = %.2e, max = %.2e\n" minimum(diff(alpha_values)) maximum(diff(alpha_values))
        end

        if !haskey(eigenvalue_data, (n, l, m, a))
            eigenvalue_data[(n, l, m, a)] = (Float64[], Float64[], Float64[])
        end

        mu_arr, alpha_arr, imag_arr = eigenvalue_data[(n, l, m, a)]
        for (i, (alpha, erg)) in enumerate(zip(alpha_values, erg_values))
            eval_count[] += 1
            if i % max(1, div(n_alpha, 5)) == 0
                @printf "  State (%d,%d,%d,%.2f): %d / %d points\n" n l m a i n_alpha
            end

            mu = alpha / (M_BH * GNew)
            push!(mu_arr, mu)
            push!(alpha_arr, alpha)
            push!(imag_arr, erg)
        end
    end
end

println("\n[SAVING] Organizing and saving computed data...")

# Organize data by quantum level and save to HDF5
data_by_level = Dict()
h5_file = joinpath(data_dir, "eigenvalues_quick.h5")

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

# Also save a summary CSV for quick reference
summary_file = joinpath(data_dir, "eigenvalues_quick_summary.csv")
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
