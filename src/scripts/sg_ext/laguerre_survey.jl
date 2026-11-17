# Accuracy of the hydrogenic radial function with FLOAT64 radii (radial_bound_NR of solve_sr_rates.jl, which
# evaluates the Laguerre polynomial as (-1)^n / n! * HypergeometricFunctions.U(-n, a + 1, x)) against a stable
# NOTE: the production rate code calls it with BigFloat radii, for which the old form is accurate; this survey shows
# what happens with Float64 arguments only.
# evaluation (three-term recurrence in BigFloat), for every level of the extended system with n > 8.
using HypergeometricFunctions, Printf
lag_old(n, α, x) = (-1) .^ n ./ factorial(big(n)) .* HypergeometricFunctions.U.(-n, α + 1, x)
function lag_new(n, α, x::AbstractVector)
    out = zeros(BigFloat, length(x))
    for (i, xi) in enumerate(x)
        xb = big(xi); L0 = big(1.0); L1 = 1 + α - xb
        if n == 0; out[i] = L0; continue; end
        for k in 1:n-1
            L0, L1 = L1, ((2k + 1 + α - xb) * L1 - (k + α) * L0) / (k + 1)
        end
        out[i] = L1
    end
    return out
end
function radial(n, l, r_a0, lag)           # r in units of a0; common normalisation
    x = 2 .* r_a0 ./ n
    return Float64.(sqrt((2 / n)^3 * factorial(big(n - l - 1)) / (2n * factorial(big(n + l)))) .* exp.(-big.(r_a0) ./ n) .* big.(x) .^ l .* lag(n - l - 1, 2l + 1, x))
end
levels = [parse.(Int, split(l, "-")) for l in readlines("/data/dust/user/samwitte/Axion_SR_gw_scan/sg_ext/levels_sgx.txt") if !isempty(strip(l))]
r = collect(10 .^ range(log10(0.03), log10(3000.0), 1200))          # in a0: the production grids start at the horizon (~0.03 a0 at alpha = 0.15)
println("# n l m  nodes  max|err|/max|R|   r below which |err| > 1e-3 max|R| [a0]   2n^2 [a0]")
bad = 0
open("survey.tsv", "w") do io
    for (n, l, m) in levels
        n > 8 || continue
        Rn = radial(n, l, r, lag_new)
        Ro = try radial(n, l, r, lag_old) catch e; fill(NaN, length(r)) end
        sc = maximum(abs.(Rn))
        err = abs.(Ro .- Rn) ./ sc
        err[isnan.(err)] .= Inf
        worst = maximum(err)
        k = findlast(e -> e > 1e-3, err)
        rb = k === nothing ? 0.0 : r[k]
        nanfrac = count(isinf, err) / length(err)
        println(io, join((n, l, m, n - l - 1, worst, rb, nanfrac), "\t"))
        if worst > 1e-3
            global bad += 1
            @printf("%2d %2d %2d  %2d   %.2e   %8.2f   %5d   %s\n", n, l, m, n - l - 1, worst, rb, 2n^2, nanfrac > 0 ? "NaN/Inf in $(round(100nanfrac, digits=1))% of points" : "")
        end
    end
end
println("levels with n > 8: ", count(x -> x[1] > 8, levels), "   with error > 1e-3 of the maximum: ", bad)
