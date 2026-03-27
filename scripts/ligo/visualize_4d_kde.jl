using DelimitedFiles
using Plots
using KernelDensity
using Statistics
using LinearAlgebra
using Distributions

# Load the data
data = readdlm("/Users/samuelwitte/Dropbox/Axion_SR/src/BH_data/GW231123_NRSur.dat")
println("Data shape: ", size(data))
println("Data columns: M1, M2, chi1, chi2")

# Transform to (M1, q, chi1, chi2) format
data_transformed = copy(data)
for i in 1:size(data_transformed, 1)
    if data_transformed[i, 2] > data_transformed[i, 1]
        # Swap masses and spins to ensure M1 >= M2
        data_transformed[i, 1], data_transformed[i, 2] = data_transformed[i, 2], data_transformed[i, 1]
        data_transformed[i, 3], data_transformed[i, 4] = data_transformed[i, 4], data_transformed[i, 3]
    end
end

# Convert M2 to q = M2/M1
kde_data = copy(data_transformed)
kde_data[:, 2] = kde_data[:, 2] ./ kde_data[:, 1]

# Extract columns
M1 = kde_data[:, 1]
q = kde_data[:, 2]
chi1 = kde_data[:, 3]
chi2 = kde_data[:, 4]

println("\nData statistics:")
println("M1: min=$(minimum(M1)), max=$(maximum(M1)), mean=$(mean(M1))")
println("q:  min=$(minimum(q)), max=$(maximum(q)), mean=$(mean(q))")
println("chi1: min=$(minimum(chi1)), max=$(maximum(chi1)), mean=$(mean(chi1))")
println("chi2: min=$(minimum(chi2)), max=$(maximum(chi2)), mean=$(mean(chi2))")

# Function to compute 2D KDE on a grid
function compute_2d_kde(x, y; n_points=100, bandwidth_factor=1.0)
    # Create grid
    x_range = range(minimum(x) - 0.1*(maximum(x)-minimum(x)),
                    maximum(x) + 0.1*(maximum(x)-minimum(x)), length=n_points)
    y_range = range(minimum(y) - 0.1*(maximum(y)-minimum(y)),
                    maximum(y) + 0.1*(maximum(y)-minimum(y)), length=n_points)

    # Use the KernelDensity package for 2D KDE
    k = kde((x, y))

    # Interpolate onto our grid
    z = zeros(n_points, n_points)
    for (i, xi) in enumerate(x_range)
        for (j, yj) in enumerate(y_range)
            z[j, i] = pdf(k, xi, yj)
        end
    end

    return x_range, y_range, z
end

# Variable names and data
variables = ["M₁ [M☉]", "q = M₂/M₁", "χ₁", "χ₂"]
var_data = [M1, q, chi1, chi2]

# Create all 6 pairwise 2D density plots
pairs = [(1,2), (1,3), (1,4), (2,3), (2,4), (3,4)]

# Create a 2x3 subplot layout
plots_list = []

for (i, j) in pairs
    x_range, y_range, z = compute_2d_kde(var_data[i], var_data[j])

    # Use log scale for better visualization
    z_pos = max.(z, 1e-10)  # Ensure all values are positive
    z_log = log10.(z_pos)

    # Clip to reasonable range
    z_min = maximum(z_log) - 4  # 4 orders of magnitude below peak
    z_log = clamp.(z_log, z_min, maximum(z_log))

    p = heatmap(collect(x_range), collect(y_range), z_log,
                xlabel=variables[i], ylabel=variables[j],
                title="$(variables[i]) vs $(variables[j])",
                c=:inferno, colorbar=true,
                size=(400, 350))

    push!(plots_list, p)
end

# Combine into single figure
combined_plot = plot(plots_list..., layout=(2,3), size=(1200, 800),
                     plot_title="4D KDE: 2D Marginal Density Projections (log scale)\n(M₁, q, χ₁, χ₂)",
                     margin=5Plots.mm)

savefig(combined_plot, "/Users/samuelwitte/Dropbox/Axion_SR/scripts/ligo/kde_2d_projections.png")
println("\nSaved combined plot to kde_2d_projections.png")

# Also create individual larger plots for better detail
println("\nCreating individual detailed plots...")

for (idx, (i, j)) in enumerate(pairs)
    x_range, y_range, z = compute_2d_kde(var_data[i], var_data[j], n_points=150)

    # Use log scale for better visualization
    z_pos = max.(z, 1e-10)  # Ensure all values are positive
    z_log = log10.(z_pos)
    z_min = maximum(z_log) - 4
    z_log = clamp.(z_log, z_min, maximum(z_log))

    p = heatmap(collect(x_range), collect(y_range), z_log,
                xlabel=variables[i], ylabel=variables[j],
                title="2D KDE Density: $(variables[i]) vs $(variables[j])",
                c=:inferno, colorbar_title="log₁₀(Density)",
                size=(600, 500))

    # Add contour lines on the original density
    contour!(p, collect(x_range), collect(y_range), z_pos,
             levels=6, color=:white, alpha=0.6, linewidth=1.5, label="")

    fname = "/Users/samuelwitte/Dropbox/Axion_SR/scripts/ligo/kde_2d_$(i)_$(j).png"
    savefig(p, fname)
    println("Saved $(variables[i]) vs $(variables[j]) to kde_2d_$(i)_$(j).png")
end

println("\nDone! All plots saved.")
