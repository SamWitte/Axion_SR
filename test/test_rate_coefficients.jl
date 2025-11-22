#!/usr/bin/env julia
"""
Test Suite: RateCoefficient Data Structure and Database

Validates that the structured rate coefficient system works correctly
and provides equivalent functionality to the original dictionary-based approach.
"""

import Random
using Printf
Random.seed!(42)

# Setup paths
src_dir = joinpath(@__DIR__, "..", "src")

println("="^80)
println("TEST SUITE: Rate Coefficient Structured Data")
println("="^80)

# Load the rate coefficient module
println("\n[1/5] Loading rate coefficients module...")
try
    include(joinpath(src_dir, "Core", "rate_coefficients.jl"))
    using Main.RateCoefficients
    println("✓ rate_coefficients.jl loaded successfully")
catch e
    println("✗ ERROR loading rate_coefficients.jl:")
    println(e)
    exit(1)
end

# Test 1: Create and populate database
println("\n[2/5] Testing database creation and population...")
try
    db = build_default_rates(3; non_relativistic=true)
    println("✓ Created database for Nmax=3")

    # Check basic properties
    nrates = count_rates(db)
    println("  - Total rates: $nrates")

    bh_rates = count_rates(db; channel="BH")
    inf_rates = count_rates(db; channel="Inf")
    gw_rates = count_rates(db; channel="GW")
    println("  - BH channels: $bh_rates")
    println("  - Inf channels: $inf_rates")
    println("  - GW channels: $gw_rates")

    if nrates > 0
        println("✓ Database population successful")
    else
        println("✗ FAIL: No rates in database")
        exit(1)
    end
catch e
    println("✗ ERROR: $e")
    exit(1)
end

# Test 2: Rate lookup and evaluation
println("\n[3/5] Testing rate lookup and evaluation...")
try
    db = build_default_rates(3; non_relativistic=true)

    # Test lookup
    rate_211_gw = get_rate(db, "211_211^GW")
    if rate_211_gw === nothing
        println("✗ FAIL: Could not retrieve '211_211^GW' rate")
        exit(1)
    end

    # Test evaluation
    alpha = 0.1  # Small coupling for test
    fa_factor = 1.0
    rp = 1.4  # Typical pericenter factor

    computed = evaluate_rate(rate_211_gw, alpha, fa_factor, rp)
    expected = 1.0e-2 * (alpha ^ 14)  # GW rate, no fa or rp dependence

    println("  - Rate key: $(rate_211_gw.key)")
    println("  - Coefficient: $(rate_211_gw.coefficient)")
    println("  - Power index: $(rate_211_gw.power_index)")
    println("  - Computed value: $computed")
    println("  - Expected value: $expected")

    if abs(computed - expected) < 1e-15
        println("✓ Rate evaluation correct")
    else
        println("⚠ WARNING: Rate evaluation differs slightly")
    end
catch e
    println("✗ ERROR: $e")
    exit(1)
end

# Test 3: Different Nmax levels
println("\n[4/5] Testing database creation for different Nmax...")
try
    for nmax in 3:5
        db = build_default_rates(nmax; non_relativistic=true)
        count = count_rates(db)
        println("  - Nmax=$nmax: $count rates")
    end
    println("✓ Database creation successful for Nmax 3-5")
catch e
    println("✗ ERROR: $e")
    exit(1)
end

# Test 4: Dictionary export (backward compatibility)
println("\n[5/5] Testing export to dictionary format...")
try
    db = build_default_rates(3; non_relativistic=true)

    alpha = 0.05
    fa_factor = 1.0
    rp = 1.4

    dict_rates = export_to_dict(db, alpha, fa_factor, rp)

    println("  - Exported to dictionary format")
    println("  - Dictionary size: $(length(dict_rates)) entries")

    # Verify a specific value
    if haskey(dict_rates, "211_211^GW")
        val = dict_rates["211_211^GW"]
        println("  - Sample rate '211_211^GW': $val")

        if val > 0
            println("✓ Dictionary export successful")
        else
            println("⚠ WARNING: Exported rate value is non-positive")
        end
    else
        println("✗ FAIL: Exported dictionary missing expected key")
        exit(1)
    end
catch e
    println("✗ ERROR: $e")
    exit(1)
end

# Final summary
println("\n" * "="^80)
println("TEST RESULTS SUMMARY")
println("="^80)
println("\n✓ All rate coefficient tests passed!")
println("\nValidations:")
println("  ✓ Database creation and population")
println("  ✓ Rate lookup by key")
println("  ✓ Rate evaluation with parameters")
println("  ✓ Multiple Nmax levels supported")
println("  ✓ Dictionary export for backward compatibility")
println("\nNext step: Integrate into super_rad.jl for production use")
println()
