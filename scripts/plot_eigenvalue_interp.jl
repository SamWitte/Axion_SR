#!/usr/bin/env julia
"""
Plot Eigenvalues from Interpolated Pre-Computed Data

Loads the pre-computed eigenvalue data from ../src/rate_sve/*.npz and uses
the interpolating functions in ../src/Numerics/rate_interpolation.jl to
evaluate the growth rate as a function of alpha at selected spin values.

Produces HDF5 output in the same format as plot_eigenvalue_quick.jl,
compatible with plot_eigenvalues.py for matplotlib plotting.
"""

using HDF5
using Printf
using NPZ

src_dir = joinpath(@__DIR__, "..", "src")
include(joinpath(src_dir, "Core", "constants.jl"))
include(joinpath(src_dir, "state_utils.jl"))
include(joinpath(src_dir, "Numerics", "rate_interpolation.jl"))
using .RateInterpolation

println("="^80)
println("EIGENVALUE PLOT FROM INTERPOLATED DATA")
println("="^80)

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
M_BH  = 10.0          # Solar masses (used only for unit conversion in rate_interpolation)
spins = [0.95]   # Spin values to evaluate

quantum_levels = [
    (2, 1, 1),
    (3, 2, 2),
    (4, 3, 3),
    (5, 4, 4),
    (6, 5, 5),
    (7, 6, 6),
    (8, 7, 7),
]

# Read the alpha grid directly from one of the data files so we stay on-grid
rate_sve_dir = joinpath(src_dir, "rate_sve")

function load_alpha_grid(n, l, m)
    state_str = format_state_for_filename(n, l, m)
    pos_file  = joinpath(rate_sve_dir, "Imag_ergC_pos_$(state_str).npz")
    neg_file  = joinpath(rate_sve_dir, "Imag_ergC_neg_$(state_str).npz")
    for fn in (pos_file, neg_file)
        if isfile(fn)
            data = npzread(fn)          # shape (n_alpha, n_spin, 3)
            return sort(unique(data[:, :, 1]))
        end
    end
    return nothing
end

# --------------------------------------------------------------------------
# Main loop
# --------------------------------------------------------------------------
eigenvalue_data = Dict()

for (n, l, m) in quantum_levels
    alpha_grid = load_alpha_grid(n, l, m)
    if isnothing(alpha_grid)
        @warn "No data files found for state ($n,$l,$m), skipping."
        continue
    end

    @printf "\nState (%d,%d,%d): %d alpha points in [%.3f, %.3f]\n" n l m length(alpha_grid) alpha_grid[1] alpha_grid[end]

    for a in spins
        mu_arr    = Float64[]
        alpha_arr = Float64[]
        rate_arr  = Float64[]

        for alph in alpha_grid
            try
                _, interp_func, _ = pre_computed_sr_rates_unified(n, l, m, alph, M_BH; cheby=true)
                rate = interp_func(a)
                push!(mu_arr,    alph / (M_BH * GNew))
                push!(alpha_arr, alph)
                push!(rate_arr,  rate)
            catch e
                @warn "  Interpolation failed at alpha=$(alph), spin=$(a): $e"
            end
        end

        @printf "  spin=%.2f: %d points evaluated\n" a length(alpha_arr)
        eigenvalue_data[(n, l, m, a)] = (mu_arr, alpha_arr, rate_arr)
    end
end

# --------------------------------------------------------------------------
# Save to HDF5 (same format as plot_eigenvalue_quick.jl)
# --------------------------------------------------------------------------
data_dir = joinpath(@__DIR__, "..", "data")
isdir(data_dir) || mkdir(data_dir)

h5_file = joinpath(data_dir, "eigenvalues_interp.h5")
h5open(h5_file, "w") do file
    for ((n, l, m, a), (mus, alphas, rates)) in eigenvalue_data
        sorted_idx     = sortperm(alphas)
        mus_s          = mus[sorted_idx]
        alphas_s       = alphas[sorted_idx]
        rates_s        = rates[sorted_idx]

        group_name = @sprintf "level_%d_%d_%d/spin_%.2f" n l m a
        file[group_name * "/mu"]              = mus_s
        file[group_name * "/alpha"]           = alphas_s
        file[group_name * "/imag_eigenvalue"] = rates_s
        file[group_name * "/M_BH"]            = M_BH
        file[group_name * "/spin"]            = a
    end
end

println("\n[SAVED] $h5_file")
println("\nNext step: python scripts/plot_eigenvalues.py --input $(h5_file) --output plts/")
println("="^80)
