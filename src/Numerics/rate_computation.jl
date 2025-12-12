"""
    rate_computation.jl

Rate coefficient computation module for superradiance calculations.

This module exports key rate computation functions that were previously defined
in solve_sr_rates.jl. It provides a clean interface for rate computations while
maintaining backward compatibility.

**Main Functions**:
- `compute_sr_rates()` - Orchestrates rate computation with interpolation
- `pre_computed_sr_rates()` - Loads and interpolates pre-computed rate data (in solve_sr_rates.jl)
- `precomputed_spin1()` - Computes spin-1 field rates analytically (in solve_sr_rates.jl)

**Note**: The actual implementations of pre_computed_sr_rates() and precomputed_spin1()
remain in solve_sr_rates.jl as they depend on its numerical infrastructure. This module
provides compute_sr_rates() which orchestrates their use.
"""

module RateComputation

using Interpolations
using NPZ
using DelimitedFiles

# Import constants from Core/constants module
include(joinpath(@__DIR__, "..", "Core/constants.jl"))
# Import state utilities
include(joinpath(@__DIR__, "..", "state_utils.jl"))

# Export public interface
export compute_sr_rates, pre_computed_sr_rates

# Helper function for loading pre-computed rate data
# This is a local implementation that the compute_sr_rates function depends on
function pre_computed_sr_rates(n, l, m, alph, M; n_high=20, n_low=20, delt_a=0.001, cheby=true)
    state_str = format_state_for_filename(n, l, m)
    if !cheby
        fn = "rate_sve/Imag_zero_$(state_str).dat"
    else
        fn = "rate_sve/Imag_zeroC_$(state_str).dat"
    end
    if isfile(fn)
        zerolist = readdlm(fn)

        itp = LinearInterpolation(zerolist[:, 1], zerolist[:, 2], extrapolation_bc=Line())
        a_mid = itp(alph)

        if a_mid .< 0
            a_mid = 0.01
        elseif a_mid .> maxSpin
            a_mid = maxSpin
        end
    else
        a_mid = 0.5
    end

    run_high = true
    run_low = true
    if (a_mid + delt_a) .> maxSpin
        run_high = false
    elseif (a_mid - delt_a) .< minSpin
        run_low = false
    end

    if run_high
        a_list_high = LinRange(a_mid + delt_a, maxSpin, n_high)
        if !cheby
            fn = "rate_sve/Imag_erg_pos_$(state_str).npz"
        else
            fn = "rate_sve/Imag_ergC_pos_$(state_str).npz"
        end
        if isfile(fn)
            file_in = npzread(fn)

            file_in[:, :, 3] .*= 2.0
            alpha_load = unique(file_in[:,:,1])
            a_load = unique(file_in[:,:,2])
            file_in[file_in[:, :, 3] .<= 0.0, 3] .= 1e-100
            itp = LinearInterpolation((alpha_load, a_load), log10.(file_in[:, :, 3] ./ (GNew * M)), extrapolation_bc=Line())
            out_high = 10 .^itp(alph, a_list_high)
        else
            a_list_high = []
            out_high = []
            run_high = false
        end
    else
        a_list_high = []
        out_high = []
    end

    if run_low
        a_list_low = LinRange(minSpin, a_mid - delt_a, n_low)
        if !cheby
            fn = "rate_sve/Imag_erg_neg_$(state_str).npz"
        else
            fn = "rate_sve/Imag_ergC_neg_$(state_str).npz"
        end
        if isfile(fn)
            file_in = npzread(fn)

            file_in[:, :, 3] .*= 2.0
            alpha_load = unique(file_in[:,:,1])
            a_load = unique(file_in[:,:,2])
            file_in[file_in[:, :, 3] .<= 0.0, 3] .= 1e-100
            itp = LinearInterpolation((alpha_load, a_load), log10.(file_in[:, :, 3] ./ (GNew * M)), extrapolation_bc=Line())
            out_low = 10 .^itp(alph, a_list_low)
        else
            a_list_low = []
            out_low = []
            run_low = false
        end
    else
        a_list_low = []
        out_low = []
    end

    return a_mid, run_high, a_list_high, out_high, run_low, a_list_low, out_low
end

"""
    compute_sr_rates(qtm_cfigs, M_BH, aBH, alph; delt_a=0.0001, cheby=true)

Compute superradiance rates for multiple quantum configurations.

Creates interpolation functions for rate coefficients as a function of black hole spin,
handling special cases and boundary conditions.

**Arguments**:
- `qtm_cfigs::Vector`: Quantum configurations [(n, l, m, alph_threshold), ...]
- `M_BH::Float64`: Black hole mass in solar masses
- `aBH::Float64`: Black hole spin parameter (dimensionless)
- `alph::Float64`: Fine structure constant (α ≈ 1/137)

**Keyword Arguments**:
- `delt_a::Float64`: Spin step for rate table (default: 0.0001)
- `cheby::Bool`: Use Chebyshev decomposition (default: true)

**Returns**:
- `SR_rates::Vector`: Rate values for current black hole spin
- `interp_functions::Vector`: Interpolation functions [(n,l,m) -> rate]
- `interp_dict::Dict`: Dictionary access to interpolation functions

**Notes**:
- Special handling for (3,2,2) quantum state with asymmetric buffer region
- Returns sentinel values (1e-100) for missing data files
- Linear interpolation with line extrapolation for robustness
"""
function compute_sr_rates(qtm_cfigs, M_BH, aBH, alph; delt_a=0.0001, cheby=true)
    # Initialize output array
    SR_rates = zeros(length(qtm_cfigs))
    interp_functions = []

    # Process each configuration
    for (idx, (n, l, m, alph_threshold)) in enumerate(qtm_cfigs)
        if m > 5  # file not computed!
            interp_func = LinearInterpolation([0.0, 1.0], [1e-100, 1e-100], extrapolation_bc=Interpolations.Line())
            push!(interp_functions, interp_func)
            SR_rates[idx] = interp_func(aBH)
            continue
        end

        # Get minimum guess and pre-computed rates
        amin_guess, run_high, a_list_high, out_high, run_low, a_list_low, out_low =
            pre_computed_sr_rates(n, l, m, alph, M_BH; n_high=200, n_low=200, delt_a=delt_a, cheby=cheby)

        # Handle upper interpolation (values >= threshold)
        if alph < alph_threshold
            if !run_high
                a_list_high = [amin_guess, amin_guess * 1.01]
                out_high = [1e-100, 1e-100]
            elseif length(a_list_high) == 1.0
                a_list_high = [a_list_high[1], a_list_high[1] * 1.01]
                out_high = [out_high[1], out_high[1]]
            end
        else
            amin_guess = 1.0
            a_list_high = [amin_guess, amin_guess * 1.01]
            out_high = [1e-100, 1e-100]
        end

        # Create upper interpolation
        itp_upper = LinearInterpolation(a_list_high, log10.(out_high), extrapolation_bc=Interpolations.Line())

        # Handle lower interpolation
        if !run_low
            a_list_low = [amin_guess * 0.99, amin_guess]
            out_low = [1e-100, 1e-100]
        elseif length(a_list_low) == 1.0
            a_list_low = [a_list_low[1] * 0.99, a_list_low[1]]
            out_low = [out_low[1], out_low[1]]
        end

        # Create lower interpolation
        itp_lower = LinearInterpolation(a_list_low, log10.(out_low), extrapolation_bc=Interpolations.Line())

        # Define boundaries
        abndry_upper = a_list_high[1]
        abndry_lower = a_list_low[end]

        # Special case for (3,2,2) which has a different buffer region behavior
        is_special_case = (n == 3 && l == 2 && m == 2)

        # Create an interpolation function with closure over the variables we need
        function create_interp_function(itp_u, itp_l, a_u, a_l, special)
            function interp_func(aspin)
                if aspin >= a_u
                    return 10.0 ^ itp_u(aspin)
                elseif aspin <= a_l
                    return -10.0 ^ itp_l(aspin)
                else
                    # Buffer region
                    if special
                        return (-10.0 ^ itp_l(aspin) + 10.0 ^ itp_u(aspin)) / 2.0
                    else
                        return 0.0
                    end
                end
            end
            return interp_func
        end

        # Create interpolation function with proper closures
        interp_func = create_interp_function(itp_upper, itp_lower, abndry_upper, abndry_lower, is_special_case)

        # Store this function for later use
        push!(interp_functions, interp_func)

        # Calculate and store the rate for the current aBH
        SR_rates[idx] = interp_func(aBH)
    end

    # Create a Dict to access interpolation functions by quantum numbers
    interp_dict = Dict()
    for (i, (n, l, m, _)) in enumerate(qtm_cfigs)
        interp_dict[(n, l, m)] = interp_functions[i]
    end

    return SR_rates, interp_functions, interp_dict
end



end # module RateComputation
