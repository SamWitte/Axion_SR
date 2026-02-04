"""
    rate_interpolation.jl

Improved rate interpolation module that eliminates flat regions near the transition.

The key insight: the data files contain artificial floor values (1e-100) near the
transition where actual rates are small. This module:
1. Filters out floored data points
2. Uses only genuine rate measurements
3. Enforces rate = 0 at the known zero-crossing point (a_mid)
4. Builds a smooth interpolation through valid points + the zero crossing
"""

module RateInterpolation

using Interpolations
using NPZ
using DelimitedFiles

# Import constants and utilities
include(joinpath(@__DIR__, "..", "Core/constants.jl"))
include(joinpath(@__DIR__, "..", "state_utils.jl"))

export compute_sr_rates_smooth, pre_computed_sr_rates_unified
export symlog, inverse_symlog

# Floor value used in the data files - anything at or below this is artificial
const RATE_FLOOR = 1e-90

"""
    symlog(x, linthresh)

Symmetric logarithm transformation for smooth sign changes.
"""
function symlog(x::Real, linthresh::Real=1e-30)
    return sign(x) * log10(1.0 + abs(x) / linthresh)
end

function inverse_symlog(y::Real, linthresh::Real=1e-30)
    return sign(y) * linthresh * (10.0^abs(y) - 1.0)
end

symlog(x::AbstractArray, linthresh::Real=1e-30) = symlog.(x, linthresh)
inverse_symlog(y::AbstractArray, linthresh::Real=1e-30) = inverse_symlog.(y, linthresh)

"""
    load_zero_crossing(state_str, rate_sve_dir; cheby=true)

Load the zero-crossing data: alpha -> spin where rate = 0
"""
function load_zero_crossing(state_str::String, rate_sve_dir::String; cheby::Bool=true)
    prefix = cheby ? "Imag_zeroC" : "Imag_zero"
    fn = joinpath(rate_sve_dir, "$(prefix)_$(state_str).dat")

    if !isfile(fn)
        return nothing
    end

    zerolist = readdlm(fn)
    itp = LinearInterpolation(zerolist[:, 1], zerolist[:, 2], extrapolation_bc=Line())
    return itp
end

"""
    extract_valid_rates_at_alpha(data_3d, alph, M; floor_threshold=RATE_FLOOR)

Extract rate vs spin at a specific alpha, filtering out floored values.
Returns (spin_array, rate_array) with only valid (non-floored) data.
"""
function extract_valid_rates_at_alpha(data_3d, alph, M; floor_threshold=RATE_FLOOR)
    alpha_vals = unique(data_3d[:, :, 1])
    spin_vals = unique(data_3d[:, :, 2])
    n_alpha = length(alpha_vals)
    n_spin = length(spin_vals)

    # Reshape rate data: (n_alpha, n_spin)
    rates = reshape(data_3d[:, :, 3], n_alpha, n_spin)

    # Apply scale factor and normalize
    rates .*= 2.0
    rates ./= (GNew * M)

    # Create 2D interpolator to get rates at our specific alpha
    itp_2d = LinearInterpolation((alpha_vals, spin_vals), rates, extrapolation_bc=Line())

    # Extract 1D slice at our alpha
    rates_at_alpha = [itp_2d(alph, a) for a in spin_vals]

    # Filter out floored values
    valid_mask = abs.(rates_at_alpha) .> floor_threshold

    return spin_vals[valid_mask], rates_at_alpha[valid_mask]
end

"""
    pre_computed_sr_rates_unified(n, l, m, alph, M; cheby=true)

Build a smooth interpolation that:
1. Uses valid (non-floored) positive rates from the positive file
2. Uses valid (non-floored) negative rates from the negative file
3. Lets the interpolation naturally cross zero where the data dictates
4. No artificial forcing - just smooth interpolation through real data
"""
function pre_computed_sr_rates_unified(n, l, m, alph, M; cheby::Bool=true)
    state_str = format_state_for_filename(n, l, m)
    rate_sve_dir = joinpath(@__DIR__, "..", "rate_sve")

    # Get zero-crossing point (for reference, but we won't force it)
    zero_itp = load_zero_crossing(state_str, rate_sve_dir; cheby=cheby)
    if isnothing(zero_itp)
        a_mid = 0.5
    else
        a_mid = zero_itp(alph)
        a_mid = clamp(a_mid, minSpin, maxSpin)
    end

    # Load data files
    prefix = cheby ? "Imag_ergC" : "Imag_erg"
    pos_file = joinpath(rate_sve_dir, "$(prefix)_pos_$(state_str).npz")
    neg_file = joinpath(rate_sve_dir, "$(prefix)_neg_$(state_str).npz")

    # Collect all valid data points
    all_spins = Float64[]
    all_rates = Float64[]

    # Add positive rate data (these are growth rates, keep as positive)
    if isfile(pos_file)
        pos_data = npzread(pos_file)
        spins_pos, rates_pos = extract_valid_rates_at_alpha(pos_data, alph, M)

        for (a, r) in zip(spins_pos, rates_pos)
            if r > 0  # Valid non-floored data
                push!(all_spins, a)
                push!(all_rates, r)  # Positive rates
            end
        end
    end

    # Add negative rate data (these are damping rates, make negative)
    if isfile(neg_file)
        neg_data = npzread(neg_file)
        spins_neg, rates_neg = extract_valid_rates_at_alpha(neg_data, alph, M)

        for (a, r) in zip(spins_neg, rates_neg)
            if r > 0  # Valid non-floored data (stored as positive)
                push!(all_spins, a)
                push!(all_rates, -r)  # Negative (damping) rates
            end
        end
    end

    # Sort by spin
    sort_idx = sortperm(all_spins)
    all_spins = all_spins[sort_idx]
    all_rates = all_rates[sort_idx]

    # Handle edge cases
    if length(all_spins) < 3
        interp_func = aspin -> 0.0
        return a_mid, interp_func, (minSpin, maxSpin)
    end

    # Remove duplicates (average rates at same spin)
    unique_spins = Float64[]
    unique_rates = Float64[]
    i = 1
    while i <= length(all_spins)
        current_spin = all_spins[i]
        sum_rate = all_rates[i]
        count = 1

        # Collect all points at approximately the same spin
        while i + count <= length(all_spins) && abs(all_spins[i + count] - current_spin) < 1e-6
            sum_rate += all_rates[i + count]
            count += 1
        end

        push!(unique_spins, current_spin)
        push!(unique_rates, sum_rate / count)
        i += count
    end

    if length(unique_spins) < 3
        interp_func = aspin -> 0.0
        return a_mid, interp_func, (minSpin, maxSpin)
    end

    # Transform to symlog space for interpolation
    # Set linthresh to the smallest non-zero rate value
    # This ensures smooth linear behavior through zero at the scale of real data
    nonzero_rates = abs.(unique_rates[unique_rates .!= 0])
    linthresh = isempty(nonzero_rates) ? 1e-30 : minimum(nonzero_rates)
    rates_symlog = symlog.(unique_rates, linthresh)

    # Create interpolation in symlog space
    itp = LinearInterpolation(unique_spins, rates_symlog, extrapolation_bc=Line())

    # Create the final interpolation function
    function interp_func(aspin)
        val_symlog = itp(aspin)
        return inverse_symlog(val_symlog, linthresh)
    end

    return a_mid, interp_func, (minimum(unique_spins), maximum(unique_spins))
end

"""
    compute_sr_rates_smooth(qtm_cfigs, M_BH, aBH, alph; cheby=true)

Drop-in replacement for compute_sr_rates() with smooth interpolation.
"""
function compute_sr_rates_smooth(qtm_cfigs, M_BH, aBH, alph; cheby::Bool=true)
    SR_rates = zeros(length(qtm_cfigs))
    interp_functions = []

    for (idx, (n, l, m, alph_threshold)) in enumerate(qtm_cfigs)
        if m > 5
            interp_func = aspin -> 1e-100
            push!(interp_functions, interp_func)
            SR_rates[idx] = interp_func(aBH)
            continue
        end

        if alph >= alph_threshold
            interp_func = aspin -> 1e-100
            push!(interp_functions, interp_func)
            SR_rates[idx] = interp_func(aBH)
            continue
        end

        a_mid, interp_func, a_range = pre_computed_sr_rates_unified(n, l, m, alph, M_BH; cheby=cheby)
        push!(interp_functions, interp_func)
        SR_rates[idx] = interp_func(aBH)
    end

    interp_dict = Dict{Tuple{Int,Int,Int}, Function}()
    for (i, (n, l, m, _)) in enumerate(qtm_cfigs)
        interp_dict[(n, l, m)] = interp_functions[i]
    end

    return SR_rates, interp_functions, interp_dict
end

end # module RateInterpolation
