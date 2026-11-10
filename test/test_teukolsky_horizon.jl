# Horizon (absorbed) GW flux from the Teukolsky solver (src/Numerics/teukolsky.jl: R^up over the source,
# Z^H and the Teukolsky-Press factor α_lmω) against the post-Newtonian absorption rate of a particle on a
# circular equatorial orbit (Tagoshi, Mano & Takasugi 1997, PTP 98, 829, eq. 52, to relative order v^5),
# and the horizon part of the cloud's annihilation rate.
#   julia test/test_teukolsky_horizon.jl
using Test, Suppressor
using SpecialFunctions: digamma
include(joinpath(@__DIR__, "..", "src", "Core", "constants.jl"))
include(joinpath(@__DIR__, "..", "src", "state_utils.jl"))
@suppress include(joinpath(@__DIR__, "..", "src", "Numerics", "gw_rates.jl"))
include(joinpath(@__DIR__, "..", "src", "Numerics", "gw_kerr_build.jl"))

"""
(dE/dt)_H / (32/5 (μ/M)² v^10) for the prograde circular equatorial orbit at r0, v = (M/r0)^{1/2},
TMT97 eq. 52 through relative order v^5 (the v^6 and higher terms are dropped).
"""
function pn_flux_horizon_circular(q, r0)
    v = 1 / sqrt(r0)
    κ = sqrt(1 - q^2)
    B(n) = q == 0 ? 0.0 : imag(digamma(3 + im * n * q / κ))
    B1, B2 = B(1), B(2)
    br = -3 // 4 * q^3 - q / 4 +
         (-q - 33 // 16 * q^3) * v^2 +
         (7 // 2 * q^4 + 2q * B2 + 1 // 2 + 6q^3 * B2 + 85 // 12 * q^2 + 3q^4 * κ + κ / 2 + 13 // 2 * κ * q^2) * v^3 +
         (-4651 // 336 * q^3 - 43 // 7 * q - 17 // 56 * q^5) * v^4 +
         (569 // 24 * q^2 + 371 // 48 * q^4 + 18q^3 * B2 - 3 // 4 * q^3 * B1 + 2κ + 2 + 33 // 4 * q^4 * κ + 6q * B2 +
          163 // 8 * κ * q^2 + q * B1) * v^5
    return v^5 * br
end

@testset "Teukolsky horizon flux" begin
    @testset "α_lmω: Schwarzschild Starobinsky constant and sign" begin
        tr = teuk_radial(-2, 2, 2, 0.1, 0.0)
        ω = 0.1; l = 2
        C2 = ((l - 1) * l * (l + 1) * (l + 2))^2 + 144ω^2
        @test teuk_horizon_alpha(tr) ≈ 256 * (2 * 2)^5 * ω * (ω^2 + 4 / 64) * (ω^2 + 16 / 64) * ω^3 / C2 rtol = 1e-12   # r_+ = 2, ε = 1/8
        trs = teuk_radial(-2, 2, 2, 0.05, 0.9)          # ω < m Ω_H: superradiant, α < 0
        @test trs.k < 0 && teuk_horizon_alpha(trs) < 0
    end

    @testset "R^up carried to the source and the Wronskian" begin
        tr = teuk_radial(-2, 2, 2, 0.02, 0.9; r_up_min=2.0)
        rp = 1 + sqrt(1 - 0.81)
        for r in (2.0, 5.0, 50.0)
            Rin, dRin = teuk_in(tr, r); Rup, dRup = teuk_up(tr, r)
            W = kerr_Delta(r, 0.9)^(-1) * (Rin * dRup - Rup * dRin)
            @test abs(W / tr.W - 1) < 1e-6
        end
        @test tr.r_m <= 2.0
    end

    @testset "circular orbit: absorbed flux vs TMT97 post-Newtonian expansion" begin
        # the PN series is truncated at relative v^5 (TMT97 eq. 52 to v^8 total); at a = 0 the next term is 172/7 v^7
        for (a, r0, tol) in ((0.0, 400.0, 1e-3), (0.9, 100.0, 1e-3), (0.5, 60.0, 2e-3), (-0.7, 100.0, 2e-3))
            PH, modes = teuk_circular_flux_H(a, r0; lmax=5)
            ratio = PH / (pn_flux_newtonian(0.0, r0) * pn_flux_horizon_circular(a, r0))
            println("  a=$a r0=$r0: P_H/(32/5 v^10) = $(PH / (32 / 5 * r0^-5)),  Teukolsky/PN = $ratio,  (2,2)/(2,1) = $(modes[(2,2)] / modes[(2,1)])")
            @test abs(ratio - 1) < tol
            a > 0 && @test PH < 0              # prograde orbit around a spinning hole at v ≪ 1: superradiant absorption
            a == 0 && @test PH > 0
            @test modes[(2, -2)] ≈ modes[(2, 2)] rtol = 1e-10
        end
    end

    @testset "cloud: direct source form vs kernel form (windowed source, both forms valid)" begin
        # With a window that vanishes quadratically at r_+, the by-parts kernel form is valid for the horizon
        # amplitude too, and must equal the direct (phase-factored jet) form; without a window only the direct
        # form is (the kernel form drops boundary terms at r_+ that R^up does not suppress).
        for (sa, sb, kind, α, a) in (((2, 1, 1), (2, 1, 1), :ann, 0.3, 0.9), ((3, 2, 2), (2, 1, 1), :tr, 0.2, 0.9),
                                     ((2, 1, 1), (2, 1, 1), :ann, 0.3, 0.0))
            ma, _ = kerr_mode(sa..., α, a); mb = sa == sb ? ma : kerr_mode(sb..., α, a)[1]
            ω̃0 = kind == :tr ? real(ma.ω) - real(mb.ω) : real(ma.ω) + real(mb.ω)
            rn, rw, θn, θw = gw_kerr_grid(ma, mb, ω̃0)
            rp = 1 + sqrt(1 - a^2); w = 0.5
            win(r) = (x = (r - rp) / w; e = exp(-x^2); (1 - e, 2x / w * e, (2 / w^2 - 4x^2 / w^2) * e))
            Tnn, Tnm, Tmm, ω̃, m̃ = gw_kerr_tetrad(ma, mb, kind, rn, θn)
            Wv = [win(r)[1] for r in rn]
            Tnn .*= Wv; Tnm .*= Wv; Tmm .*= Wv
            Shat, Ph = gw_kerr_source_hat(ma, mb, kind, rn, θn; window=win)
            l = max(2, abs(m̃))
            tr = teuk_radial(-2, l, m̃, ω̃, a; r_in_max=1.05rn[end], tstops_in=rn, r_up_min=rn[1], tstops_up=rn)
            ZHk = teuk_Z_grid(tr, Tnn, Tnm, Tmm, rn, rw, θn, θw; horizon=true)
            ZHd = teuk_Z_hat(tr, Shat, Ph, rn, rw, θn, θw; horizon=true)
            ZIk = teuk_Z_grid(tr, Tnn, Tnm, Tmm, rn, rw, θn, θw)
            ZId = teuk_Z_hat(tr, Shat, Ph, rn, rw, θn, θw)
            println("  windowed $(kind) $(join(sa))x$(join(sb)) α=$α a=$a ℓ=$l: Z^H kernel/direct - 1 = $(abs(ZHk / ZHd - 1)),  Z^∞: $(abs(ZIk / ZId - 1))")
            @test abs(ZHk / ZHd - 1) < 1e-4
            @test abs(ZIk / ZId - 1) < 1e-6
        end
    end

    @testset "cloud: horizon part of the rate" begin
        # 2p x 2p at a = 0.9, α = 0.3: the level and the graviton are superradiant (k < 0), the hole feeds
        # the radiation (Γ_H < 0) at the per-cent level; the direct-form P_∞ reproduces the kernel form
        ma, _ = kerr_mode(2, 1, 1, 0.3, 0.9)
        Γ, P, modes, ω̃, m̃, ΓH, modesH, chk = gw_kerr_rate(ma, ma, :ann; horizon=true)
        println("  211x211 α=0.3 a=0.9: Γ_∞ = $Γ, Γ_H/Γ_∞ = $(ΓH / Γ), chk = $chk")
        @test -0.05 < ΓH / Γ < -0.01
        @test abs(chk) < 1e-5
        @test abs(Γ / gw_kerr_rate(ma, ma, :ann)[1] - 1) < 1e-6        # the infinity flux does not depend on horizon=true
        # Schwarzschild, α = 0.1: absorption (Γ_H > 0), small
        ma0, _ = kerr_mode(2, 1, 1, 0.1, 0.0)
        Γ0, _, _, _, _, ΓH0, _, chk0 = gw_kerr_rate(ma0, ma0, :ann; horizon=true)
        println("  211x211 α=0.1 a=0: Γ_H/Γ_∞ = $(ΓH0 / Γ0), chk = $chk0")
        @test 0 < ΓH0 / Γ0 < 1e-2 && abs(chk0) < 1e-5
        # transition 322 -> 211 at α = 0.2, a = 0.9: m̃ = 1 graviton far below m̃ Ω_H, horizon part ~ -4e-4
        mb, _ = kerr_mode(3, 2, 2, 0.2, 0.9); mc, _ = kerr_mode(2, 1, 1, 0.2, 0.9)
        Γt, _, _, _, _, ΓHt, _, chkt = gw_kerr_rate(mb, mc, :tr; horizon=true)
        println("  322->211 α=0.2 a=0.9: Γ_H/Γ_∞ = $(ΓHt / Γt), chk = $chkt")
        @test -1e-3 < ΓHt / Γt < 0 && abs(chkt) < 1e-10
    end
end
