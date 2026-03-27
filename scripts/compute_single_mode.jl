#!/usr/bin/env julia
"""
Compute a single eigenvalue at a specific (n, l, m) mode and (alpha, a) point
using the Chebyshev method.

Usage: edit the parameters below and run with `julia compute_single_mode.jl`
"""

src_dir = joinpath(@__DIR__, "..", "src")
include(joinpath(src_dir, "Core/constants.jl"))
include(joinpath(src_dir, "heunc.jl"))
include(joinpath(src_dir, "solve_sr_rates.jl"))

# ── Parameters ──────────────────────────────────────────────────────────────
n, l, m = 8, 7, 7          # quantum numbers (877 mode)
alpha    = 0.3              # dimensionless gravitational coupling α = μ M G_N
a        = 0.01              # dimensionless BH spin
M_BH     = 1.0              # BH mass in solar masses

# Chebyshev solver settings
Npoints  = 70               # number of Chebyshev nodes
Iter     = 30               # max Newton iterations
L        = 8                # spherical harmonic truncation
prec     = 200              # BigFloat precision (bits)
cvg_acc  = 1e-10
der_acc  = 1e-20
debug    = true
# ─────────────────────────────────────────────────────────────────────────────

mu = alpha / (M_BH * GNew)

println("="^60)
println("Single-mode Chebyshev eigenvalue computation")
println("  Mode   : (n,l,m) = ($n,$l,$m)")
println("  α      = $alpha")
println("  a      = $a")
println("  μ      = $mu eV")
println("="^60)

@time wR, wI = eigensys_Cheby(
    M_BH, a, mu, n, l, m;
    prec     = prec,
    L        = L,
    Npoints  = Npoints,
    Iter     = Iter,
    debug    = debug,
    der_acc  = der_acc,
    cvg_acc  = cvg_acc,
    return_wf = false,
    sfty_run  = false,
)

println("\n── Result ──────────────────────────────────────────────")
println("  Re(ω)  = $wR")
println("  Im(ω)  = $wI")
println("  Im(ω)/α = $(wI / alpha)   [superradiance rate in units of 1/M]")
println("="^60)
