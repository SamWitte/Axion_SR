# Radial Teukolsky solver and point-particle source (src/Numerics/teukolsky.jl).
#   julia test/test_teukolsky.jl
using Test, QuadGK
include(joinpath(@__DIR__, "..", "src", "Numerics", "teukolsky.jl"))

@testset "Teukolsky" begin

    @testset "angular eigenvalue convention and normalisation" begin
        for (s, l, m, c) in ((-2, 2, 2, 0.0), (-2, 3, 2, 0.37), (-2, 4, -3, -0.2), (-2, 2, 1, 0.8))
            λ, S = teuk_angular(s, l, m, c)
            # angular Teukolsky equation for s = -2, ST03 eq. 15, with λ in Teukolsky's convention
            res = maximum(range(0.2, π - 0.2, length=30)) do θ
                S0 = S(θ); S1 = S(θ; theta_derivative=1); S2 = S(θ; theta_derivative=2)
                abs(S2 + cot(θ) * S1 + (-c^2 * sin(θ)^2 - (m - 2cos(θ))^2 / sin(θ)^2 + 4c * cos(θ) - 2 + 2m * c + λ) * S0)
            end
            @test res < 1e-9 * max(1, abs(λ))
            @test 2π * quadgk(θ -> S(θ)^2 * sin(θ), 0, π)[1] ≈ 1 atol = 1e-10
            c == 0 && @test λ ≈ (l - s) * (l + s + 1)
        end
    end

    @testset "boundary series solve the radial equation" begin
        for (a, l, m, ω) in ((0.0, 2, 2, 0.3), (0.9, 2, 2, 0.3), (0.9, 3, -1, 0.05), (0.99, 2, 2, 1e-3), (0.5, 5, 4, 2.0))
            tr = teuk_radial(-2, l, m, ω, a)
            rp, rm = kerr_rpm(a)
            # horizon series beyond its matching point (inside its radius of convergence) vs the integrated ODE
            r1 = rp + 0.3 * (rp - rm)
            Rs, dRs = teuk_eval_horizon(tr.fh, -2, m, ω, a, r1)
            Ro, dRo, _ = teuk_in(tr, r1)
            @test abs(Rs / Ro - 1) < 1e-10
            @test abs(dRs / dRo - 1) < 1e-10
            # asymptotic series inside r_inf vs the integrated ODE
            r2 = 0.7 * tr.r_inf
            Us, dUs = teuk_eval_infinity(tr.finf, -2, m, ω, a, r2)
            Uo, dUo, _ = teuk_up(tr, r2)
            @test abs(Us / Uo - 1) < 1e-10
            @test abs(dUs / dUo - 1) < 1e-10
            # Wronskian constant over [r_m, r_inf]
            @test tr.W_spread < 1e-7
            # R'' from the ODE agrees with a finite difference of R'
            r3 = 3 * rp
            h = 1e-4
            _, _, d2 = teuk_in(tr, r3)
            @test abs((teuk_in(tr, r3 + h)[2] - teuk_in(tr, r3 - h)[2]) / (2h) / d2 - 1) < 1e-6
        end
    end

    @testset "BigFloat path agrees with Float64" begin
        setprecision(BigFloat, 80) do
            tr64 = teuk_radial(-2, 2, 2, 0.3, 0.9)
            trB = teuk_radial(-2, 2, 2, big"0.3", big"0.9"; T=BigFloat, rtol=big"1e-16")
            @test abs(trB.W / tr64.W - 1) < 1e-7       # the Float64 W carries the ~1e-8 seen in W_spread
            @test abs(teuk_in(trB, big(10.0))[1] / teuk_in(tr64, 10.0)[1] - 1) < 1e-10
            @test abs(teuk_up(trB, trB.r_inf / 2)[1] / teuk_up(tr64, tr64.r_inf / 2)[1] - 1) < 1e-10
        end
    end

    @testset "circular equatorial orbit vs post-Newtonian luminosity (ST03 eqs. 178, 181)" begin
        for (a, r0, lmax) in ((0.0, 1000.0, 4), (0.0, 100.0, 6), (0.9, 1000.0, 4), (-0.7, 1000.0, 4))
            P, modes = teuk_circular_flux(a, r0; lmax=lmax)
            ratio = P / pn_flux_newtonian(a, r0) / pn_flux_circular(a, r0)
            println("  circular orbit a=$a r0=$r0: P/PN_series - 1 = $(ratio - 1)")
            @test abs(ratio - 1) < 3e-8
            @test modes[(2, -2)] ≈ modes[(2, 2)] rtol = 1e-12       # (m, ω) -> (-m, -ω) symmetry
            @test modes[(3, -3)] ≈ modes[(3, 3)] rtol = 1e-12
            # Newtonian limit: the quadrupole flux sits in (2, ±2), l = 3 enters at O(v^2)
            @test abs(2 * modes[(2, 2)] / pn_flux_newtonian(a, r0) - 1) < 10 / r0
            @test modes[(3, 3)] / modes[(2, 2)] < 5 / r0
        end
    end
end
