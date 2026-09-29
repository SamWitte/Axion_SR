# Validation of src/Numerics/gw_rates.jl
#   julia test/test_gw_rates.jl
using Test
using QuadGK
using SpecialFunctions: sphericalbesselj
using ForwardDiff

include(joinpath(@__DIR__, "..", "src", "Core", "constants.jl"))
include(joinpath(@__DIR__, "..", "src", "state_utils.jl"))
include(joinpath(@__DIR__, "..", "src", "Numerics", "gw_rates.jl"))

evalR(R, r) = sum(c * r^p for (p, c) in zip(R.pows, R.coefs)) * exp(-R.beta * r)

@testset "GW rates" begin

    @testset "hydrogen radial normalisation + Schrodinger" begin
        for (n, l) in [(2, 1), (3, 1), (3, 2), (6, 4), (12, 5), (18, 17)]
            R = gw_hydrogen_radial(n, l)
            @test quadgk(r -> evalR(R, r)^2 * r^2, 0, Inf)[1] ≈ 1 atol = 1e-10
            # -1/2 (R'' + 2R'/r - l(l+1)R/r^2) - R/r = -1/(2n^2) R
            r0 = 0.7 * n
            f(r) = evalR(R, r)
            d1 = ForwardDiff.derivative(f, r0)
            d2 = ForwardDiff.derivative(r -> ForwardDiff.derivative(f, r), r0)
            lhs = -0.5 * (d2 + 2d1 / r0 - l * (l + 1) * f(r0) / r0^2) - f(r0) / r0
            @test lhs ≈ -f(r0) / (2n^2) rtol = 1e-8
        end
    end

    @testset "spherical harmonics orthonormal" begin
        xs, ws = gauss(40)
        for (l1, l2, m) in [(3, 3, 2), (5, 3, 1), (7, 7, -4), (17, 17, 17)]
            v = 2π * sum(w * gw_theta_lm(l1, m, acos(x)) * gw_theta_lm(l2, m, acos(x)) for (x, w) in zip(xs, ws))
            @test v ≈ (l1 == l2 ? 1.0 : 0.0) atol = 1e-12
        end
    end

    @testset "gradient formula vs finite differences" begin
        x0 = [0.9, -1.3, 0.6] .* 2
        for (n, l, m) in [(2, 1, 1), (3, 2, 1), (5, 4, 3), (6, 3, 2)]
            g = gw_grad_psi(n, l, m, x0)
            h = 1e-5
            fd = [(gw_psi(n, l, m, x0 .+ h .* e) - gw_psi(n, l, m, x0 .- h .* e)) / (2h)
                  for e in ([1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0])]
            @test maximum(abs.(g .- fd)) < 1e-8 * max(1, maximum(abs.(fd)))
        end
    end

    @testset "closed-form radial integral" begin
        for (p, J, β, k) in [(3, 1, 1.0, 2.0), (4, 4, 1.5, 3.0), (6, 2, 0.8, 5.0), (5, 5, 1.0, 7.0), (9, 3, 0.7, 4.0)]
            num = quadgk(r -> r^p * exp(-β * r) * sphericalbesselj(J, k * r), 0, Inf, rtol=1e-12)[1]
            @test gw_Itilde(p, J, β / k) / k^(p + 1) ≈ num rtol = 1e-8
        end
    end

    @testset "2p annihilation vs Arvanitaki+ 2015 Table VI (flat space, all orders)" begin
        # dP/dΩ N^-2 = α^18 G/r_g^4 [6α^3+40α-3(α^2+4)^2 atan(2/α)]^2 (28cos2θ+cos4θ+35)/(2^24 π (α^2+4)^4)
        # solver rate = X/(2α^2),  X = ∫dΩ (dP/dΩ N^-2) r_g^4/G
        ang = 2π * (16 / 5 + 32 + 16)
        arv(α) = α^18 * (6α^3 + 40α - 3(α^2 + 4)^2 * atan(2 / α))^2 * ang / (2^24 * π * (α^2 + 4)^4) / (2α^2)
        for α in (0.02, 0.1, 0.3)
            ours = gw_annihilation_rate_flat((2, 1, 1), (2, 1, 1), α)
            println("  2p ann α=$α  ours=$(ours)  Arvanitaki=$(arv(α))  ratio=$(ours / arv(α))")
            @test ours ≈ arv(α) rtol = 1e-6
        end
        C, p = gw_annihilation_leading((2, 1, 1), (2, 1, 1))
        println("  2p leading: C=$C p=$p")
        @test p == 16
    end

    @testset "3d, 4f annihilation leading terms vs Arvanitaki+ 2015 Table VI" begin
        # 3d: α^20 sin^4(28cos2θ+cos4θ+35)/(2^4 3^16 π);  4f: α^24 sin^8(...) 5^2/(2^44 π)
        I4 = 2π * 2 * (8 / 9 + 32 / 7 - 16 + 32 / 3 + 8)
        f8(x) = (1 - x^2)^4 * (8x^4 + 48x^2 + 8)
        I8 = 2π * quadgk(f8, -1, 1)[1]
        C3, p3 = gw_annihilation_leading((3, 2, 2), (3, 2, 2))
        C4, p4 = gw_annihilation_leading((4, 3, 3), (4, 3, 3))
        ref3 = I4 / (2^4 * 3.0^16 * π) / 2
        # (the PDF text of Table VI garbles the 4f prefactor; it is 5^2/(4^22 π), cf. 3d: 1/(2^4 3^16 π))
        ref4 = I8 * 5.0^2 / (2.0^44 * π) / 2
        println("  3d: C=$C3 p=$p3 (ref $(ref3), p=18);  4f: C=$C4 p=$p4 (ref $(ref4), p=22)")
        @test p3 == 18
        @test C3 ≈ ref3 rtol = 1e-3
        @test p4 == 22
        @test C4 ≈ ref4 rtol = 1e-3
    end

    @testset "current multipole vs brute-force 3D quadrature" begin
        # J^(m) = (1/L) ∫ (x × p)·∇f,  f = r^L Y*_Lm,  p = (i/2)(ψa ∇ψb* - ψb* ∇ψa)
        sa, sb, L = (3, 2, 2), (2, 1, 1), 2
        m = sa[3] - sb[3]
        function integrand(r, θ, φ)
            x = r .* [sin(θ)cos(φ), sin(θ)sin(φ), cos(θ)]
            ψa = gw_psi(sa..., x); ψb = gw_psi(sb..., x)
            ga = gw_grad_psi(sa..., x); gb = gw_grad_psi(sb..., x)
            pv = (im / 2) .* (ψa .* conj.(gb) .- conj(ψb) .* ga)
            xp = [x[2] * pv[3] - x[3] * pv[2], x[3] * pv[1] - x[1] * pv[3], x[1] * pv[2] - x[2] * pv[1]]
            fgrad = ForwardDiff.gradient(y -> begin
                rr = sqrt(sum(abs2, y)); tt = acos(y[3] / rr); pp = atan(y[2], y[1])
                rr^L * gw_theta_lm(L, m, tt) * cos(m * pp)
            end, x) .- im .* ForwardDiff.gradient(y -> begin
                rr = sqrt(sum(abs2, y)); tt = acos(y[3] / rr); pp = atan(y[2], y[1])
                rr^L * gw_theta_lm(L, m, tt) * sin(m * pp)
            end, x)
            return sum(xp .* fgrad) * r^2 * sin(θ) / L
        end
        rs, wr = gauss(60, 0, 60.0); ts, wt = gauss(30, 0, π); ps, wp = gauss(30, 0, 2π)
        J = sum(wr[i] * wt[j] * wp[k] * integrand(rs[i], ts[j], ps[k]) for i in 1:60, j in 1:30, k in 1:30)
        ch = gw_transition_channels(sa, sb)
        Kc = [K for (LL, kind, K) in ch if LL == L && kind == :current][1]
        println("  current quadrupole |J|^2: brute=$(abs2(J))  table=$(Kc / (2 * gw_ccurr(L)))")
        @test abs2(J) ≈ Kc / (2 * gw_ccurr(L)) rtol = 1e-4
    end

    @testset "transitions vs exact quadrupole formula (Arvanitaki & Dubovsky 2011, eq. 39)" begin
        # dN1/dt = N1 N0 (2 G Δω^5 / 5) Q_ij Q_ij*,  Q_ij = μ ∫ ψa ψb* (x_i x_j - δ_ij r^2/3)
        # -> solver rate = (2/5) δ^5 α^{-2} |Q|^2 (Q in Bohr units). Brute-force Q on a grid.
        α = 0.01
        for (sa, sb) in [((6, 4, 4), (5, 4, 4)), ((4, 3, 3), (3, 1, 1)), ((5, 2, 2), (3, 2, 2))]
            rmax = 6.0 * max(sa[1], sb[1])^2
            rs, wr = gauss(90, 0, rmax); ts, wt = gauss(40, 0, π); ps, wp = gauss(40, 0, 2π)
            Q = zeros(ComplexF64, 3, 3)
            for i in 1:90, j in 1:40, k in 1:40
                x = rs[i] .* [sin(ts[j])cos(ps[k]), sin(ts[j])sin(ps[k]), cos(ts[j])]
                w = wr[i] * wt[j] * wp[k] * rs[i]^2 * sin(ts[j]) * gw_psi(sa..., x) * conj(gw_psi(sb..., x))
                Q .+= w .* (x * x' .- (rs[i]^2 / 3) .* [1.0 0 0; 0 1 0; 0 0 1])
            end
            δ = α^2 / 2 * (1 / sb[1]^2 - 1 / sa[1]^2)
            @test gw_transition_rate(gw_transition_channels(sa, sb), α, δ) ≈ 2 / 5 * δ^5 * α^(-2) * sum(abs2, Q) rtol = 1e-8
        end
    end

    @testset "2p annihilation vs Arvanitaki & Dubovsky 2011 eq. (44)" begin
        # dP/dΩ = N^2 9π G α^18/(2^26 r_g^4) (35 + 28cos2θ + cos4θ)  ->  C α^16 with C = X/2
        C, p = gw_annihilation_leading((2, 1, 1), (2, 1, 1))
        @test p == 16
        @test C ≈ 9π / 2.0^26 * 2π * (16 / 5 + 32 + 16) / 2 rtol = 1e-4
    end

    @testset "literature l=1 override" begin
        C, p = gw_literature_annihilation((2, 1, 1), (2, 1, 1))
        @test p == 14 && C ≈ (484 + 9π^2) / 46080
        # same n-scaling as the flat-space coefficients
        f(a, b) = gw_annihilation_leading(a, b)[1]
        lit(a, b) = gw_literature_annihilation(a, b)[1]
        @test lit((3, 1, 1), (5, 1, 1)) / lit((2, 1, 1), (2, 1, 1)) ≈ f((3, 1, 1), (5, 1, 1)) / f((2, 1, 1), (2, 1, 1)) rtol = 1e-4
        @test gw_literature_annihilation((2, 1, 1), (3, 2, 2)) === nothing
    end

    @testset "transitions vs literature" begin
        # 322 -> 211: Baryakhtar+ 2021 Table IV: 5e-6 α^10
        ch = gw_transition_channels((3, 2, 2), (2, 1, 1))
        α = 0.01
        δ = α^2 / 2 * (1 / 4 - 1 / 9)
        r = gw_transition_rate(ch, α, δ)
        println("  322->211 channels $(ch);  rate/α^10 = $(r / α^10)  (Baryakhtar+: 5e-6)")
        @test r / α^10 ≈ 5e-6 rtol = 0.1
        # 6g -> 5g vs Arvanitaki+ 2015 Table VII (flat-space kinetic prescription):
        #   dP/dΩ/(N N') = 2^28 3^4 5^5 α^12 sin^4θ/(11^22 π) G/r_g^4  ->  rate = X/(α^2 δ)
        ch6b = gw_transition_channels((6, 4, 4), (5, 4, 4))
        δ6 = α^2 / 2 * (1 / 25 - 1 / 36)
        arv = 2.0^28 * 3^4 * 5^5 / (11.0^22 * π) * (32π / 15) * α^12 / (α^2 * δ6)
        ours6 = gw_transition_rate(ch6b, α, δ6)
        println("  6g->5g (644->544): ours rate/α^8 = $(ours6 / α^8),  flat-kinetic (Arvanitaki+) = $(arv / α^8)")
        # the flat-space kinetic-stress prescription misses the Newtonian stress; for
        # same-l transitions it is exactly 4x the quadrupole formula
        @test arv / ours6 ≈ 4 rtol = 1e-6
    end
end
