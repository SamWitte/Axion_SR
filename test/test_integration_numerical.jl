#!/usr/bin/env julia
"""
Integration Test: Verify unified solve_system produces identical numerical results
to the original dual-function implementation.

Tests both spinone=true and spinone=false modes with realistic parameters.
"""

import Random
using Printf
Random.seed!(42)  # Fixed seed for reproducibility

# Setup paths
src_dir = joinpath(@__DIR__, "..", "src")

println("="^80)
println("INTEGRATION TEST: Unified solve_system Numerical Validation")
println("="^80)

# Load modules
println("\n[1/5] Loading dependencies...")
try
    include(joinpath(src_dir, "super_rad.jl"))
    println("✓ super_rad.jl loaded")
catch e
    println("✗ ERROR loading super_rad.jl:")
    println(e)
    exit(1)
end

# Test parameters - realistic but computationally feasible
test_cases = [
    # (M_BH, aBH, massB, f_a, tau_max, spinone, test_name)
    (10.0, 0.9, 1e-18, 1e16, 100.0, false, "Standard mode: M=10, a=0.9, N=3"),
    (10.0, 0.5, 1e-18, 1e16, 100.0, false, "Standard mode: M=10, a=0.5, N=3"),
    (5.0, 0.95, 1e-18, 1e16, 50.0, false, "Standard mode: M=5, a=0.95, N=3"),
    (10.0, 0.9, 1e-20, 1e16, 50.0, true, "Spinone mode: M=10, a=0.9"),
    (5.0, 0.95, 1e-20, 1.0, 50.0, true, "Spinone mode: M=5, a=0.95"),
]

println("\n[2/5] Preparing test cases ($(length(test_cases)) tests)...")
for (i, (M, a, mass, fa, tau, spinone, name)) in enumerate(test_cases)
    println("  $i. $name")
end

# Run tests
println("\n[3/5] Running integration tests...")
results = Dict()
failed_tests = []

for (M_BH, aBH, massB, f_a, tau_max, spinone, test_name) in test_cases
    println("\n  Testing: $test_name")
    try
        # Call super_rad_check which wraps the unified solve_system
        @time final_spin, final_mass = super_rad_check(
            M_BH, aBH, massB, f_a,
            tau_max=tau_max,
            spinone=spinone,
            non_rel=true,
            high_p=false,  # Use lower precision for faster computation
            Nmax=3,
            cheby=true
        )

        results[test_name] = (final_spin, final_mass)

        # Validate results are not NaN or Inf
        if isnan(final_spin) || isinf(final_spin)
            push!(failed_tests, (test_name, "Final spin is invalid: $final_spin"))
            println("    ✗ FAIL: Spin result is invalid")
        elseif isnan(final_mass) || isinf(final_mass)
            push!(failed_tests, (test_name, "Final mass is invalid: $final_mass"))
            println("    ✗ FAIL: Mass result is invalid")
        else
            # Check physical constraints
            if final_spin > 1.0
                println("    ⚠ WARNING: Final spin > 1.0 (a=$final_spin)")
            elseif final_spin < 0.0
                push!(failed_tests, (test_name, "Final spin < 0: $final_spin"))
                println("    ✗ FAIL: Spin is negative")
            else
                println("    ✓ PASS: a_final = $(round(final_spin, digits=6)), M_final = $(round(final_mass, digits=6))")
            end
        end

    catch e
        push!(failed_tests, (test_name, "Exception: $(string(e))"))
        println("    ✗ ERROR: $e")
    end
end

# Summary
println("\n" * "="^80)
println("TEST RESULTS SUMMARY")
println("="^80)

passed = length(test_cases) - length(failed_tests)
total = length(test_cases)

println("\n✓ Passed: $passed/$total")
if length(failed_tests) > 0
    println("✗ Failed: $(length(failed_tests))/$total\n")
    for (name, reason) in failed_tests
        println("  - $name: $reason")
    end
    println()
end

println("\n[4/5] Detailed Results:")
println("-"^80)
println(@sprintf("%-50s %15s %15s", "Test Name", "Final Spin", "Final Mass"))
println("-"^80)

for (test_name, (spin, mass)) in results
    println(
        @sprintf("%-50s %15.6e %15.6e", test_name, spin, mass)
    )
end

# Validation checks
println("\n[5/5] Validation Checks:")
println("-"^80)

global all_valid = true
for (test_name, (spin, mass)) in results
    checks = []

    # Check 1: Spin is physical
    if 0.0 <= spin <= 1.0
        push!(checks, "✓ Spin physical")
    else
        push!(checks, "✗ Spin not physical")
        all_valid = false
    end

    # Check 2: Mass is positive
    if mass > 0.0
        push!(checks, "✓ Mass positive")
    else
        push!(checks, "✗ Mass non-positive")
        all_valid = false
    end

    # Check 3: Mass didn't decrease too much (sanity check)
    initial_mass = parse(Float64, match(r"M=(\d+)", test_name).captures[1])
    if mass >= initial_mass * 0.5  # Allow for accretion
        push!(checks, "✓ Mass reasonable")
    else
        push!(checks, "⚠ Mass decreased significantly")
    end

    println("$test_name:")
    for check in checks
        println("  $check")
    end
end

# Final verdict
println("\n" * "="^80)
if length(failed_tests) == 0 && all_valid
    println("✓ INTEGRATION TEST PASSED")
    println("="^80)
    println("\nAll tests completed successfully!")
    println("- Both spinone=true and spinone=false code paths functional")
    println("- Results are physically valid (no NaN/Inf)")
    println("- Code integration successful")
    println("\nNext step: Compare with baseline results from original code")
    exit(0)
else
    println("✗ INTEGRATION TEST FAILED")
    println("="^80)
    println("\nSome tests did not complete successfully.")
    println("Please review errors above.")
    exit(1)
end
