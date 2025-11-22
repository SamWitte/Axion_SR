#!/usr/bin/env julia
"""
Detailed Performance Analysis - Phase 7.1

Provides detailed timing analysis for AxionSR functions
using multiple iterations and warm-up runs.
"""

import Random
using Printf, Statistics
Random.seed!(42)

src_dir = joinpath(@__DIR__, "..", "src")

println("="^80)
println("DETAILED PERFORMANCE ANALYSIS: AxionSR Module")
println("="^80)

# Load module
println("\n[1/4] Loading module...")
include(joinpath(src_dir, "AxionSR.jl"))
using .AxionSR
println("✓ Module loaded")

# Test parameters with varying complexity
test_cases = [
    ("Light: M=1, a=0.9, τ=5", 1.0, 0.9, 1e-20, 1e16, 5.0),
    ("Medium: M=5, a=0.7, τ=10", 5.0, 0.7, 1e-19, 1e16, 10.0),
    ("Heavy: M=10, a=0.95, τ=20", 10.0, 0.95, 1e-18, 1e16, 20.0),
    ("Extreme: M=20, a=0.99, τ=50", 20.0, 0.99, 1e-17, 1e16, 50.0),
]

println("\n[2/4] Running warm-up passes...")
for (desc, M, a, mass, fa, tau) in test_cases[1:1]
    try
        @time super_rad_check(M, a, mass, fa, tau_max=tau, spinone=false, non_rel=true, Nmax=3, cheby=true)
    catch e
        println("Warm-up error: $e")
    end
end
println("✓ Warm-up complete")

println("\n[3/4] DETAILED TIMING ANALYSIS")
println("="^80)

results_summary = []

for (desc, M, a, mass, fa, tau) in test_cases
    println("\n$desc")
    println("-"^50)

    times = []
    for iter in 1:3
        try
            t_start = time()
            final_a, final_M = super_rad_check(
                M, a, mass, fa,
                tau_max=tau,
                spinone=false,
                non_rel=true,
                Nmax=3,
                cheby=true
            )
            t_end = time()
            elapsed = (t_end - t_start)
            push!(times, elapsed)
            println("  Iteration $iter: $(round(elapsed, digits=4)) s")
        catch e
            println("  Iteration $iter: ERROR - $e")
        end
    end

    if !isempty(times)
        mean_time = mean(times)
        min_time = minimum(times)
        max_time = maximum(times)
        std_time = std(times)

        println("  Summary:")
        println("    Mean:   $(round(mean_time, digits=4)) s")
        println("    Min:    $(round(min_time, digits=4)) s")
        println("    Max:    $(round(max_time, digits=4)) s")
        println("    StdDev: $(round(std_time, digits=4)) s")

        push!(results_summary, (desc, mean_time, M, a, tau))
    end
end

println("\n[4/4] PERFORMANCE SUMMARY")
println("="^80)

if !isempty(results_summary)
    println("\nComparative Analysis:")
    println("-"^80)
    println(@sprintf("%-30s %12s %8s %8s %8s", "Test Case", "Time (s)", "M", "a", "τ"))
    println("-"^80)

    for (desc, time, M, a, tau) in results_summary
        println(@sprintf("%-30s %12.4f %8.1f %8.2f %8.1f", desc, time, M, a, tau))
    end

    println("\nScaling Analysis:")
    println("-"^80)

    # Check if timing scales with tau_max
    if length(results_summary) >= 2
        time_ratios = []
        tau_ratios = []

        for i in 1:length(results_summary)-1
            desc1, t1, M1, a1, tau1 = results_summary[i]
            desc2, t2, M2, a2, tau2 = results_summary[i+1]

            # If similar mass and spin, check tau scaling
            if abs(a1 - a2) < 0.05
                time_ratio = t2 / t1
                tau_ratio = tau2 / tau1
                push!(time_ratios, time_ratio)
                push!(tau_ratios, tau_ratio)

                linear_scaling = time_ratio / tau_ratio
                println("  Transition $i→$(i+1):")
                println("    Time ratio:   $(round(time_ratio, digits=3))x")
                println("    τ ratio:      $(round(tau_ratio, digits=3))x")
                println("    Linear scaling factor: $(round(linear_scaling, digits=3))")

                if linear_scaling < 1.1
                    println("    → Excellent linear scaling ✓")
                elseif linear_scaling < 1.3
                    println("    → Good scaling with minor overhead")
                else
                    println("    → Sublinear scaling detected")
                end
            end
        end
    end

    # Overall statistics
    all_times = [time for (_, time, _, _, _) in results_summary]
    println("\nOverall Statistics:")
    println("-"^80)
    println("  Fastest: $(round(minimum(all_times), digits=4)) s")
    println("  Slowest: $(round(maximum(all_times), digits=4)) s")
    println("  Average: $(round(mean(all_times), digits=4)) s")
    println("  Ratio:   $(round(maximum(all_times)/minimum(all_times), digits=1))x")
end

println("\n" * "="^80)
println("ANALYSIS COMPLETE")
println("="^80)
println("\nKey Findings:")
println("1. Check if timing scales linearly with τ_max (tau parameter)")
println("2. Check if timing varies significantly with mass and spin")
println("3. Identify which parameter regime is slowest")
println("4. These results can guide optimization efforts")

exit(0)
