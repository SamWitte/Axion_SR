"""
Test Suite for Rate Computation Module

Tests the newly extracted rate computation module functionality.
"""

src_dir = joinpath(@__DIR__, "..", "src")
include(joinpath(src_dir, "super_rad.jl"))
using .RateComputation

println("="^80)
println("RATE COMPUTATION MODULE TESTS")
println("="^80)

# Test 1: Module can be imported and compute_sr_rates function exists
println("\n[TEST 1] Module import and function existence")
try
    @assert isdefined(RateComputation, :compute_sr_rates)
    println("✓ PASSED: compute_sr_rates is exported from RateComputation module")
catch e
    println("✗ FAILED: $(string(e))")
end

# Test 2: compute_sr_rates function signature and basic execution
println("\n[TEST 2] Function signature validation")
try
    # Small test with minimal quantum configs
    qtm_cfigs = [(2, 1, 1, 0.01), (3, 2, 1, 0.02)]
    M_BH = 1.0
    aBH = 0.5
    alph = 0.002

    SR_rates, interp_functions, interp_dict = compute_sr_rates(
        qtm_cfigs, M_BH, aBH, alph; delt_a=0.0001, cheby=true
    )

    @assert length(SR_rates) == 2 "SR_rates should have length 2"
    @assert length(interp_functions) == 2 "interp_functions should have length 2"
    @assert length(interp_dict) == 2 "interp_dict should have length 2"

    println("✓ PASSED: compute_sr_rates returns correct structure")
catch e
    println("✗ FAILED: $(string(e))")
end

# Test 3: Interpolation dict access
println("\n[TEST 3] Dictionary access by quantum numbers")
try
    qtm_cfigs = [(2, 1, 1, 0.01)]
    M_BH = 1.0
    aBH = 0.5
    alph = 0.002

    SR_rates, interp_functions, interp_dict = compute_sr_rates(
        qtm_cfigs, M_BH, aBH, alph
    )

    @assert haskey(interp_dict, (2, 1, 1)) "Dictionary should have key (2,1,1)"
    func = interp_dict[(2, 1, 1)]
    @assert isa(func, Function) "Dictionary value should be a function"

    # Test calling the interpolation function
    result = func(0.5)
    @assert isa(result, Real) "Function should return a real number"

    println("✓ PASSED: Interpolation functions accessible via dictionary")
catch e
    println("✗ FAILED: $(string(e))")
end

# Test 4: Rate values are reasonable
println("\n[TEST 4] Rate value reasonableness")
try
    qtm_cfigs = [(2, 1, 1, 0.01), (3, 2, 1, 0.02), (3, 2, 2, 0.03)]
    M_BH = 1.0
    aBH = 0.5
    alph = 0.002

    SR_rates, _, _ = compute_sr_rates(qtm_cfigs, M_BH, aBH, alph)

    # Rates should be real and mostly positive (with some special cases negative)
    for rate in SR_rates
        @assert isa(rate, Real) "All rates should be real"
    end

    println("✓ PASSED: All rates are real numbers")
catch e
    println("✗ FAILED: $(string(e))")
end

# Test 5: Multiple spin values
println("\n[TEST 5] Multiple spin parameter values")
try
    qtm_cfigs = [(2, 1, 1, 0.01)]
    M_BH = 1.0
    alph = 0.002

    spins = [0.1, 0.3, 0.5, 0.7, 0.9]
    rates = []

    for spin in spins
        sr, _, _ = compute_sr_rates(qtm_cfigs, M_BH, spin, alph)
        push!(rates, sr[1])
    end

    @assert length(rates) == 5 "Should compute rates for all spins"

    println("✓ PASSED: Computed rates for multiple spin values")
catch e
    println("✗ FAILED: $(string(e))")
end

# Test 6: High m values (m > 5) return sentinel values
println("\n[TEST 6] High m values handling")
try
    qtm_cfigs = [(2, 1, 6, 0.01)]  # m = 6, which is > 5
    M_BH = 1.0
    aBH = 0.5
    alph = 0.002

    SR_rates, _, _ = compute_sr_rates(qtm_cfigs, M_BH, aBH, alph)

    # Should return sentinel value 1e-100
    @assert SR_rates[1] ≈ 1e-100 "Should return sentinel value for m > 5"

    println("✓ PASSED: Sentinel values returned for m > 5")
catch e
    println("✗ FAILED: $(string(e))")
end

# Test 7: Special case (3,2,2) is handled
println("\n[TEST 7] Special case (3,2,2) handling")
try
    qtm_cfigs = [(3, 2, 2, 0.01)]
    M_BH = 1.0
    aBH = 0.5
    alph = 0.002

    SR_rates, _, interp_dict = compute_sr_rates(qtm_cfigs, M_BH, aBH, alph)

    # Check that interpolation function exists
    @assert haskey(interp_dict, (3, 2, 2)) "Should handle (3,2,2) special case"

    println("✓ PASSED: Special case (3,2,2) handled correctly")
catch e
    println("✗ FAILED: $(string(e))")
end

# Test 8: Backward compatibility - pre_computed_sr_rates still in solve_sr_rates
println("\n[TEST 8] Backward compatibility with solve_sr_rates")
try
    # These functions should still be available from solve_sr_rates
    @assert isdefined(Main, :pre_computed_sr_rates) "pre_computed_sr_rates should still exist"
    @assert isdefined(Main, :precomputed_spin1) "precomputed_spin1 should still exist"

    println("✓ PASSED: Backward compatibility maintained")
catch e
    println("✗ FAILED: $(string(e))")
end

println("\n" * "="^80)
println("TEST SUMMARY")
println("="^80)
println("All rate computation tests completed!")
println("="^80)
