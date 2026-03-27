#!/usr/bin/env julia
"""
Plot eigenvalues vs alpha for a fixed spin value.

Loads all heunc rate files (files with "C" in the name) and plots:
- Positive eigenvalues as a function of alpha (superradiant growth rates)
- Negative eigenvalues as a function of alpha (decay rates)

Usage:
    julia plot_eigenvalues_vs_alpha.jl [--spin=0.9] [--output=plts/]
"""

using NPZ
using Interpolations
using Plots
using LaTeXStrings
using ArgParse

# Path to rate files
const RATE_SVE_DIR = joinpath(@__DIR__, "..", "src", "rate_sve")

"""
    parse_state_from_filename(filename::String) -> (n, l, m)

Parse quantum numbers from a rate file name like "Imag_ergC_pos_211.npz"
"""
function parse_state_from_filename(filename::String)
    # Extract the state string (e.g., "211" from "Imag_ergC_pos_211.npz")
    m = match(r"Imag_ergC_(?:pos|neg)_(\d+)\.npz", filename)
    if m === nothing
        return nothing
    end
    state_str = m.captures[1]
    if length(state_str) != 3
        return nothing
    end
    n = parse(Int, state_str[1])
    l = parse(Int, state_str[2])
    m_val = parse(Int, state_str[3])
    return (n, l, m_val)
end

"""
    load_rate_data(filename::String) -> (alpha_vals, spin_vals, eigenvalues)

Load rate data from an NPZ file.
Returns the alpha values, spin values, and eigenvalue matrix.
"""
function load_rate_data(filename::String)
    data = npzread(filename)
    # Data is (n_alpha, n_spin, 3) where [:,:,1]=alpha, [:,:,2]=spin, [:,:,3]=eigenvalue
    alpha_vals = unique(data[:, :, 1])
    spin_vals = unique(data[:, :, 2])
    eigenvalues = data[:, :, 3]
    return alpha_vals, spin_vals, eigenvalues
end

# Number of interpolation points for smooth curves
const N_ALPHA_POINTS = 500

"""
    get_eigenvalue_vs_alpha(filename::String, fixed_spin::Float64) -> (alpha_vals, eigenvalues)

Load rate data and interpolate to get eigenvalue vs alpha at a fixed spin.
Uses a fine grid of alpha values for smooth curves.
Interpolates in log10 space for proper behavior across orders of magnitude.
"""
function get_eigenvalue_vs_alpha(filename::String, fixed_spin::Float64)
    alpha_vals_raw, spin_vals, eigenvalues = load_rate_data(filename)

    # Replace zeros/negatives with small value before taking log
    eigenvalues_safe = copy(eigenvalues)
    eigenvalues_safe[eigenvalues_safe .<= 0] .= 1e-100

    # Interpolate in log10 space (as done in rate_computation.jl)
    log_eigenvalues = log10.(eigenvalues_safe)
    itp = LinearInterpolation((alpha_vals_raw, spin_vals), log_eigenvalues, extrapolation_bc=Line())

    # Create fine alpha grid for smooth curves
    alpha_min, alpha_max = extrema(alpha_vals_raw)
    alpha_vals = collect(LinRange(alpha_min, alpha_max, N_ALPHA_POINTS))

    # Evaluate at fixed spin for all alpha values, then convert back from log
    eig_at_spin = [10.0^itp(a, fixed_spin) for a in alpha_vals]

    return alpha_vals, eig_at_spin
end

"""
    find_all_rate_files() -> Dict

Find all heunc rate files and organize them by quantum state and sign.
Returns Dict mapping (n,l,m) -> Dict("pos" => filename, "neg" => filename)
"""
function find_all_rate_files()
    files = readdir(RATE_SVE_DIR)

    # Filter for heunc files (contain "C" in name)
    heunc_files = filter(f -> occursin("Imag_ergC_", f) && endswith(f, ".npz"), files)

    # Organize by state and sign
    rate_files = Dict{Tuple{Int,Int,Int}, Dict{String,String}}()

    for f in heunc_files
        state = parse_state_from_filename(f)
        if state === nothing
            continue
        end

        sign = occursin("_pos_", f) ? "pos" : "neg"

        if !haskey(rate_files, state)
            rate_files[state] = Dict{String,String}()
        end
        rate_files[state][sign] = joinpath(RATE_SVE_DIR, f)
    end

    return rate_files
end

# Plot style constants
const ALPHA_MIN = 0.1           # Minimum alpha to display
const EIGENVALUE_MIN = 1e-15    # Minimum eigenvalue to display

"""
    create_eigenvalue_plots(fixed_spin::Float64, output_dir::String)

Create plots of eigenvalues vs alpha for all quantum states at fixed spin.
"""
function create_eigenvalue_plots(fixed_spin::Float64, output_dir::String)
    mkpath(output_dir)

    rate_files = find_all_rate_files()
    println("Found $(length(rate_files)) quantum states")

    # Collect data for positive and negative eigenvalues
    # Apply filters: alpha >= ALPHA_MIN, eigenvalue >= EIGENVALUE_MIN
    pos_data = Dict{Tuple{Int,Int,Int}, Tuple{Vector{Float64}, Vector{Float64}}}()
    neg_data = Dict{Tuple{Int,Int,Int}, Tuple{Vector{Float64}, Vector{Float64}}}()

    for (state, files) in rate_files
        n, l, m = state

        # Process positive eigenvalues
        if haskey(files, "pos")
            try
                alpha_vals, eig_vals = get_eigenvalue_vs_alpha(files["pos"], fixed_spin)
                # Filter: alpha >= ALPHA_MIN, eigenvalue >= EIGENVALUE_MIN
                valid_mask = (alpha_vals .>= ALPHA_MIN) .& (eig_vals .>= EIGENVALUE_MIN)
                if any(valid_mask)
                    pos_data[state] = (alpha_vals[valid_mask], eig_vals[valid_mask])
                end
            catch e
                println("Warning: Could not process positive file for state $state: $e")
            end
        end

        # Process negative eigenvalues
        if haskey(files, "neg")
            try
                alpha_vals, eig_vals = get_eigenvalue_vs_alpha(files["neg"], fixed_spin)
                # Filter: alpha >= ALPHA_MIN, eigenvalue >= EIGENVALUE_MIN
                valid_mask = (alpha_vals .>= ALPHA_MIN) .& (eig_vals .>= EIGENVALUE_MIN)
                if any(valid_mask)
                    neg_data[state] = (alpha_vals[valid_mask], eig_vals[valid_mask])
                end
            catch e
                println("Warning: Could not process negative file for state $state: $e")
            end
        end
    end

    println("States with visible positive eigenvalues: $(length(pos_data))")
    println("States with visible negative eigenvalues: $(length(neg_data))")

    # Create distinct color palette for the states that will actually be plotted
    # Use a qualitative palette that's easier to distinguish
    function get_colors(n)
        if n <= 10
            return palette(:tab10, n)
        elseif n <= 20
            return palette(:tab20, n)
        else
            return palette(:rainbow, n)
        end
    end

    # Plot positive eigenvalues
    println("\nCreating positive eigenvalue plot...")
    sorted_states_pos = sort(collect(keys(pos_data)))
    colors_pos = get_colors(length(sorted_states_pos))

    p_pos = plot(
        xlabel=L"\alpha",
        ylabel=L"\Gamma",
        title="Positive Eigenvalues (a = $(fixed_spin))",
        xscale=:identity,  # Linear x-axis
        yscale=:log10,
        legend=:outerright,
        legendfontsize=8,
        labelfontsize=14,
        tickfontsize=11,
        size=(1000, 600),
        dpi=150,
        framestyle=:box,
        minorticks=true,
        tick_direction=:in,
        xlims=(ALPHA_MIN, :auto),
        ylims=(EIGENVALUE_MIN, :auto),
        xtickfontsize=12,
        ytickfontsize=12,
        guidefontsize=14
    )

    for (i, state) in enumerate(sorted_states_pos)
        n, l, m = state
        alpha_vals, eig_vals = pos_data[state]
        plot!(p_pos, alpha_vals, eig_vals,
              label="($n,$l,$m)",
              color=colors_pos[i],
              linewidth=2)
    end

    # Adjust tick lengths
    plot!(p_pos, ticklength=8)

    pos_filename = joinpath(output_dir, "eigenvalues_positive_a$(fixed_spin).png")
    savefig(p_pos, pos_filename)
    println("  Saved: $pos_filename")

    # Plot negative eigenvalues (decay rates - stored as magnitudes)
    println("Creating negative eigenvalue plot...")
    sorted_states_neg = sort(collect(keys(neg_data)))
    colors_neg = get_colors(length(sorted_states_neg))

    p_neg = plot(
        xlabel=L"\alpha",
        ylabel=L"|\Gamma|",
        title="Negative Eigenvalues (Decay Rates, a = $(fixed_spin))",
        xscale=:identity,  # Linear x-axis
        yscale=:log10,
        legend=:outerright,
        legendfontsize=8,
        labelfontsize=14,
        tickfontsize=11,
        size=(1000, 600),
        dpi=150,
        framestyle=:box,
        minorticks=true,
        tick_direction=:in,
        xlims=(ALPHA_MIN, :auto),
        ylims=(EIGENVALUE_MIN, :auto),
        xtickfontsize=12,
        ytickfontsize=12,
        guidefontsize=14
    )

    for (i, state) in enumerate(sorted_states_neg)
        n, l, m = state
        alpha_vals, eig_vals = neg_data[state]
        plot!(p_neg, alpha_vals, eig_vals,
              label="($n,$l,$m)",
              color=colors_neg[i],
              linewidth=2)
    end

    # Adjust tick lengths
    plot!(p_neg, ticklength=8)

    neg_filename = joinpath(output_dir, "eigenvalues_negative_a$(fixed_spin).png")
    savefig(p_neg, neg_filename)
    println("  Saved: $neg_filename")

    println("\nDone! Created 2 plots in $output_dir")

    return p_pos, p_neg
end

function parse_commandline()
    s = ArgParseSettings(description="Plot eigenvalues vs alpha at fixed spin")

    @add_arg_table! s begin
        "--spin"
            help = "Fixed black hole spin value"
            arg_type = Float64
            default = 0.9
        "--output", "-o"
            help = "Output directory for plots"
            arg_type = String
            default = "plts/"
    end

    return parse_args(s)
end

function main()
    args = parse_commandline()

    fixed_spin = args["spin"]
    output_dir = args["output"]

    println("="^60)
    println("EIGENVALUE VS ALPHA PLOTTING SCRIPT")
    println("="^60)
    println("Fixed spin: $fixed_spin")
    println("Output directory: $output_dir")
    println()

    create_eigenvalue_plots(fixed_spin, output_dir)

    println("\n" * "="^60)
    println("COMPLETE")
    println("="^60)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
