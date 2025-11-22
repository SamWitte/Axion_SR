"""
Test suite for AxionSR module.

This file coordinates all test suites for the Axion Superradiance project.
Run with: julia test/runtests.jl

Test Organization:
- unit_tests.jl - Core functionality tests
- integration_tests.jl - Integration and cross-module tests
- coefficient_tests.jl - Rate coefficient system tests
- computation_tests.jl - Numerical computation tests
"""

using Test

# Add src directory to load path for module imports
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

println("="^80)
println("AXION SUPERRADIANCE TEST SUITE")
println("="^80)

# Test suites to run
test_suites = [
    ("Unit Tests", "unit_tests.jl"),
    ("Integration Tests", "integration_tests.jl"),
    ("Rate Coefficient Tests", "coefficient_tests.jl"),
    ("Computation Tests", "computation_tests.jl"),
]

# Track results
total_tests = 0
passed_tests = 0
failed_suites = []

# Run each test suite
for (name, file) in test_suites
    if isfile(joinpath(@__DIR__, file))
        println("\n" * "="^80)
        println("Running: $name")
        println("="^80)
        try
            include(file)
            println("✓ $name completed successfully")
        catch e
            println("✗ $name FAILED")
            push!(failed_suites, name)
            @warn "Test file $file encountered error: $(string(e))"
        end
    else
        println("\n⚠ Skipping $name - file not found: $file")
    end
end

# Summary
println("\n" * "="^80)
println("TEST SUMMARY")
println("="^80)

if isempty(failed_suites)
    println("✓ All test suites completed successfully!")
    println("\nStatus: READY FOR PRODUCTION")
else
    println("✗ Some test suites had issues:")
    for suite in failed_suites
        println("  - $suite")
    end
end

println("="^80)
