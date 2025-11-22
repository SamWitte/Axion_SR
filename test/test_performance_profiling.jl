#!/usr/bin/env julia
"""
Phase 7: Performance Profiling Test

Profiles the AxionSR module to identify performance bottlenecks.
Provides detailed timing information for major computational functions.

Runtime: ~15-30 minutes depending on system
"""

import Random
using Printf, Statistics, BenchmarkTools
Random.seed!(42)

println("="^80)
println("PERFORMANCE PROFILING: AxionSR Module")
println("="^80)

# Setup paths
src_dir = joinpath(@__DIR__, "..", "src")

println("\n[1/5] Loading AxionSR module...")
try
    include(joinpath(src_dir, "AxionSR.jl"))
    using .AxionSR
    println("✓ AxionSR module loaded successfully")
catch e
    println("✗ ERROR loading AxionSR module:")
    println(e)
    exit(1)
end

# Profiling parameters
println("\n[2/5] Setting up profiling tests...")
profile_params = [
    (1.0, 0.9, 1e-20, 1e16, 10.0, "Minimum mass, high spin"),
    (5.0, 0.7, 1e-19, 1e16, 20.0, "Medium mass, medium spin"),
    (10.0, 0.95, 1e-18, 1e16, 50.0, "High mass, high spin"),
]

results = []

println("\n[3/5] Running performance benchmarks...")
println("-"^80)

for (M_BH, aBH, massB, f_a, tau_max, desc) in profile_params
    println("\nTest: $desc")
    println("  M=$M_BH, a=$aBH, μ=$massB, f_a=$f_a, τ=$tau_max")

    try
        # Benchmark super_rad_check
        bench = @benchmarkable super_rad_check(
            $M_BH, $aBH, $massB, $f_a,
            tau_max=$tau_max,
            spinone=false,
            non_rel=true,
            Nmax=3,
            cheby=true
        ) seconds=5

        result = run(bench)

        println("  Timing statistics:")
        println("    Mean:    $(round(mean(result.times)/1e9, digits=4)) s")
        println("    Median:  $(round(median(result.times)/1e9, digits=4)) s")
        println("    Min:     $(round(minimum(result.times)/1e9, digits=4)) s")
        println("    Max:     $(round(maximum(result.times)/1e9, digits=4)) s")
        println("    StdDev:  $(round(std(result.times)/1e9, digits=4)) s")

        push!(results, (desc, mean(result.times)/1e9))

    catch e
        println("  ✗ ERROR: $e")
    end
end

# Profile major functions
println("\n[4/5] Function-level profiling...")
println("-"^80)

# Test solve_system directly
println("\nProfile: solve_system() - Core evolution solver")
try
    M_BH, aBH, massB, f_a = 5.0, 0.8, 1e-19, 1e16
    tau_max = 10.0

    bench = @benchmarkable solve_system(
        $massB, $f_a, $aBH, $M_BH, $tau_max,
        spinone=false,
        non_rel=true,
        high_p=false,
        Nmax=3,
        cheby=true
    ) seconds=5

    result = run(bench)
    println("  Mean time: $(round(mean(result.times)/1e9, digits=4)) s")

catch e
    println("  Could not profile solve_system: $e")
end

# Summary and recommendations
println("\n[5/5] PROFILING SUMMARY")
println("="^80)

println("\nTiming Results:")
for (desc, time) in results
    println("  $desc: $(round(time, digits=4)) s")
end

if !isempty(results)
    avg_time = mean([time for (_, time) in results])
    min_time = minimum([time for (_, time) in results])
    max_time = maximum([time for (_, time) in results])

    println("\nAggregate Statistics:")
    println("  Average time: $(round(avg_time, digits=4)) s")
    println("  Min time:     $(round(min_time, digits=4)) s")
    println("  Max time:     $(round(max_time, digits=4)) s")
    println("  Ratio max/min: $(round(max_time/min_time, digits=2))x")
end

println("\n" * "="^80)
println("PROFILING COMPLETE")
println("="^80)
println("\nNext steps:")
println("1. Review timing results above")
println("2. Identify functions taking >1 second")
println("3. Profile those functions with @profile macro")
println("4. Implement optimizations based on profiles")
println("5. Re-run benchmarks to measure improvements")

exit(0)
