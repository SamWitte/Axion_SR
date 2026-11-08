# Source term of the s = -2 Teukolsky equation: the kernel form used for extended sources
# (teuk_source_kernels / teuk_Z_extended) against the direct operator form of ST03 (eqs. 9-12, 18)
# and Teukolsky's original NP source T_4, on a smooth test source off the equator with a ≠ 0.
#   julia test/test_teukolsky_source.jl
using Test, QuadGK
include(joinpath(@__DIR__, "..", "src", "Numerics", "teukolsky.jl"))

a, ω, m, l = 0.9, 0.37, 2, 3
r0, w = 8.0, 1.0
bump(r) = exp(-(r - r0)^2 / w^2)
# arbitrary smooth components, regular at the poles for azimuthal number m = 2: a quantity of spin
# weight s behaves as sin(θ/2)^|m+s| cos(θ/2)^|m-s|, so T_nn (s = 0) must vanish like sin^2 θ while
# T_m̄m̄ (s = -2) may be finite at θ = 0. Without this the integrations by parts leave boundary terms.
Tnn(r, θ) = bump(r) * sin(θ)^2 * (1 + 0.3cos(θ)) * (1 + 0.5im)
Tnm(r, θ) = bump(r) * sin(θ) * (1 + 0.5cos(θ)) * (0.7 - 0.4im)
Tmm(r, θ) = bump(r) * sin(θ)^2 * (1 - 0.3cos(θ)) * (0.3 + 0.6im)

@testset "Teukolsky source" begin
    @testset "ST03 operator form = Teukolsky's NP source T_4" begin
        B2 = teuk_source_B2(Tnn, Tnm, Tmm, a, ω, m)
        T4 = teuk_source_T4(Tnn, Tnm, Tmm, a, ω, m)
        worst = 0.0; scale = 0.0
        for r in (6.0, 8.0, 9.5), θ in (0.4, 1.1, 1.9, 2.8)
            b = B2(r, θ); t = T4(r, θ)
            worst = max(worst, abs(b - t)); scale = max(scale, abs(t))
        end
        println("  |B'_2 + B'*_2 - T_4| / max|T_4| = ", worst / scale)
        @test worst < 1e-9 * scale
    end

    @testset "kernel form = direct form of Z^∞" begin
        tr = teuk_radial(-2, l, m, ω, a)
        xr, wr = gauss(80); rn = r0 .+ 5w .* xr; rw = 5w .* wr
        xθ, wθ = gauss(40); θn = π / 2 .* (xθ .+ 1); θw = π / 2 .* wθ
        Zk = teuk_Z_extended(tr, Tnn, Tnm, Tmm, rn, rw, θn, θw)
        Zd = teuk_Z_direct(tr, Tnn, Tnm, Tmm, rn, rw, θn, θw)
        println("  Z_kernels = ", Zk, "   Z_direct = ", Zd, "   ratio - 1 = ", Zk / Zd - 1)
        @test abs(Zk / Zd - 1) < 1e-6
    end
end
