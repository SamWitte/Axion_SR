#!/usr/bin/env julia
"""
Quick Smoke Test: Verify both code paths (spinone and standard) can execute
without errors with minimal computational cost.
"""

src_dir = joinpath(@__DIR__, "..", "src")

println("="^80)
println("SMOKE TEST: Quick Validation of Both Code Paths")
println("="^80)

println("\n[1/3] Loading super_rad...")
include(joinpath(src_dir, "super_rad.jl"))
println("✓ Loaded successfully")

println("\n[2/3] Testing Standard Mode (spinone=false)...")
try
    # Use minimal parameters for speed
    spin_in, mass_in = 0.9, 10.0
    println("  Input: M=$mass_in, a=$spin_in")

    # Short integration time for quick test
    result_spin, result_mass = super_rad_check(
        mass_in, spin_in, 1e-18, 1e16,
        tau_max=10.0,      # Very short time
        spinone=false,
        non_rel=true,
        high_p=false,
        Nmax=3
    )

    println("  Output: a_final=$(round(result_spin, digits=4)), M_final=$(round(result_mass, digits=4))")

    if !isnan(result_spin) && !isinf(result_spin) && result_spin >= 0
        println("  ✓ PASS: Standard mode executed successfully")
    else
        println("  ✗ FAIL: Invalid result - spin=$result_spin")
        exit(1)
    end

catch e
    println("  ✗ ERROR: $e")
    exit(1)
end

println("\n[3/3] Testing Spinone Mode (spinone=true)...")
try
    spin_in, mass_in = 0.95, 5.0
    println("  Input: M=$mass_in, a=$spin_in")

    result_spin, result_mass = super_rad_check(
        mass_in, spin_in, 1e-20, 1.0,
        tau_max=10.0,      # Very short time
        spinone=true
    )

    println("  Output: a_final=$(round(result_spin, digits=4)), M_final=$(round(result_mass, digits=4))")

    if !isnan(result_spin) && !isinf(result_spin) && result_spin >= 0
        println("  ✓ PASS: Spinone mode executed successfully")
    else
        println("  ✗ FAIL: Invalid result - spin=$result_spin")
        exit(1)
    end

catch e
    println("  ✗ ERROR: $e")
    exit(1)
end

println("\n" * "="^80)
println("✓ SMOKE TEST PASSED")
println("="^80)
println("\nBoth code paths executed successfully without errors:")
println("  ✓ Standard mode (multiple quantum levels)")
println("  ✓ Spinone mode (single quantum level)")
println("\nIntegration is working correctly!")
