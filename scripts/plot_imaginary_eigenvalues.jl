"""
    plot_imaginary_eigenvalues.jl

Compute imaginary parts of eigenvalues across a range of α values for a 1 solar mass BH.

For each spin a, computes Im(ω) as a function of α = μ * M * G_N from 0.05 to ~1.0
Uses adaptive resolution near the zero crossing.

Output: data/imaginary_eigenvalues_<level>.h5 (HDF5 format)
"""

using HDF5, DelimitedFiles, Suppressor

# Load necessary files
script_dir = dirname(@__FILE__)
proj_dir = dirname(script_dir)
cd(proj_dir)
include(joinpath(proj_dir, "src/Constants.jl"))
include(joinpath(proj_dir, "src/solve_sr_rates.jl"))

# Physical constants
M_BH = 1.0  # 1 solar mass

# Quantum numbers
levels = [
    (2, 1, 1),   # 211
    (2, 1, -1),  # 21-1
    (3, 1, 1),   # 311
    (3, 2, 2),   # 322
]

# Black hole spins
spins = [0.1, 0.3, 0.5, 0.7, 0.9, 0.99]

# Alpha range
alpha_min = 0.05
alpha_max = 1.0

# Create output directory
mkpath("data")

# ============================================================================
println("\n" * "="^70)
println("Computing Imaginary Eigenvalue Parts")
println("="^70)

for (n, l, m) in levels
    level_str = "$(n)$(l)$(abs(m))"

    println("\nLevel ($n, $l, $m) - Output: data/imaginary_eigenvalues_$(level_str).h5")

    # Open HDF5 file for writing
    h5open("data/imaginary_eigenvalues_$(level_str).h5", "w") do file
        for spin in spins
            println("  Spin a = $spin")

            # Create dataset with adaptive grid
            alpha_vals = Float64[]
            im_vals = Float64[]

            # Start with coarse grid
            alphas = collect(range(alpha_min, alpha_max, length=20))

            for alpha in alphas
                mu = alpha / (M_BH * GNew)

                try
                    # Compute with suppressed output
                    im_part = @suppress find_im_part(mu, M_BH, spin, n, l, m)

                    push!(alpha_vals, alpha)
                    push!(im_vals, im_part)

                    println("    α = $(round(alpha; digits=4)): Im(ω) = $(round(im_part; sigdigits=5))")

                    # If approaching zero, add extra points
                    if length(im_vals) > 1 && sign(im_vals[end]) != sign(im_vals[end-1]) && abs(im_vals[end-1]) > 1e-15
                        println("      >> Sign flip detected - adding refined points")

                        alpha_prev = alpha_vals[end-1]
                        refined_alphas = collect(range(alpha_prev, alpha, length=8))[2:end-1]

                        for alpha_r in refined_alphas
                            mu_r = alpha_r / (M_BH * GNew)
                            im_r = @suppress find_im_part(mu_r, M_BH, spin, n, l, m)
                            push!(alpha_vals, alpha_r)
                            push!(im_vals, im_r)
                            println("      α = $(round(alpha_r; digits=5)): Im(ω) = $(round(im_r; sigdigits=5))")
                        end
                    end

                catch e
                    println("    α = $(round(alpha; digits=4)): Failed - $(typeof(e).__name__)")
                    break
                end
            end

            # Sort and store in HDF5
            if length(alpha_vals) > 0
                sort_idx = sortperm(alpha_vals)
                alpha_sorted = alpha_vals[sort_idx]
                im_sorted = im_vals[sort_idx]

                # Create group for this spin
                grp = create_group(file, "spin_$(spin)")
                grp["alpha"] = alpha_sorted
                grp["im_part"] = im_sorted
                grp["n"] = n
                grp["l"] = l
                grp["m"] = m

                println("    >> Stored $(length(alpha_sorted)) points")
            end
        end
    end

    println("  ✓ Saved to data/imaginary_eigenvalues_$(level_str).h5")
end

println("\n" * "="^70)
println("Computation Complete!")
println("="^70 * "\n")
