# Quasi-bound Kerr modes (src/Numerics/kerr_modes.jl): eigenfrequencies, Klein-Gordon residuals, norm.
#   julia test/test_kerr_modes.jl
using Test, ForwardDiff
include(joinpath(@__DIR__, "..", "src", "Numerics", "teukolsky.jl"))
include(joinpath(@__DIR__, "..", "src", "Numerics", "kerr_modes.jl"))

@testset "Kerr modes" begin
    cases = ((2, 1, 1, 0.1, 0.9), (2, 1, 1, 0.3, 0.9), (3, 2, 2, 0.2, 0.5), (4, 3, 3, 0.4, 0.95), (3, 1, 1, 0.1, 0.7), (2, 1, 1, 0.05, 0.0))
    for (n, l, m, alph, a) in cases
        md, res = kerr_mode(n, l, m, alph, a)
        @testset "($n,$l,$m) α=$alph a=$a" begin
            @test abs(res) < 1e-10
            @test real(md.ω) < alph
            # hydrogenic expansion to O(α^5): the difference is O(α^6)
            @test abs(real(md.ω) - real(km_omega_guess(n, l, m, alph, a))) < 2alph^6 * 10
            # growing (superradiant) iff m Ω_H > Re ω
            @test sign(imag(md.ω)) == sign(m * kerr_OmegaH(a) - real(md.ω)) || a == 0
            # radial Klein-Gordon residual in Boyer-Lindquist form at several radii
            rp = 1 + sqrt(1 - a^2)
            for r in (rp + 1e-2, 3.0, n^2 / alph^2, 4n^2 / alph^2)
                f, f1, f2 = kerr_radial_f(md, r)
                ω = md.ω; Δ = r^2 - 2r + a^2; K = (r^2 + a^2) * ω - a * m; dK = 2r * ω; dΔ = 2(r - 1)
                ph = exp(-im * ω * kerr_rstar(r, a) + im * m * kerr_rsharp(r, a))
                R = ph * f; R1 = ph * (f1 - im * K / Δ * f)
                R2 = ph * ((-im * K / Δ) * (f1 - im * K / Δ * f) + f2 - im * (dK / Δ - K * dΔ / Δ^2) * f - im * K / Δ * f1)
                lhs = Δ * (dΔ * R1 + Δ * R2) + (K^2 - Δ * (alph^2 * r^2 + a^2 * ω^2 - 2a * m * ω + md.Λ)) * R
                @test abs(lhs) < 1e-11 * (abs(md.Λ * Δ * R) + abs(K^2 * R))
                # kerr_radial_BL agrees with the phase-converted f
                Rb, dRb = kerr_radial_BL(md, r)
                @test Rb ≈ R rtol = 1e-12
                @test dRb ≈ R1 rtol = 1e-12
            end
            # angular equation with the spheroidal eigenvalue
            c2 = a^2 * (md.ω^2 - alph^2)
            for θ in (0.4, 1.0, 1.6, 2.5)
                S, dS = kerr_angular(md, θ)
                d2S = sum(b * ForwardDiff.derivative(t -> ForwardDiff.derivative(s -> km_theta_lm(lp, m, s), t), θ) for (lp, b) in md.sph)
                lhs = d2S + cot(θ) * dS + (c2 * cos(θ)^2 - m^2 / sin(θ)^2 + md.Λ) * S
                @test abs(lhs) < 1e-11 * abs(md.Λ * S)
            end
            # unit Klein-Gordon norm, and its non-relativistic limit ∫|ψ|² d³x
            @test kerr_mode_norm(md) ≈ 1 rtol = 1e-8
            rn, rw = km_log_panels(rp, rp + 80n^2 / alph^2, 400)
            xθ, wθ = gauss(48); θs = π / 2 .* (xθ .+ 1); wθs = π / 2 .* wθ
            ang = sum(wθs .* sin.(θs) .* abs2.(first.(kerr_angular.(Ref(md), θs))))
            Nnr = 2π * ang * sum(rw .* [abs2(kerr_radial_BL(md, r)[1]) * r^2 for r in rn])
            @test abs(Nnr - 1) < 3alph^2 + 1e-3
            @test 2π * ang ≈ 1 atol = 1e-8
        end
    end

    @testset "mirror symmetry (a, m) -> (-a, -m)" begin
        m1, _ = kerr_mode(3, 2, 2, 0.2, 0.7)
        m2, _ = kerr_mode(3, 2, -2, 0.2, -0.7)
        @test m1.ω ≈ m2.ω rtol = 1e-12
        @test abs(kerr_radial_f(m1, 30.0)[1]) ≈ abs(kerr_radial_f(m2, 30.0)[1]) rtol = 1e-10
    end

    # superradiance tables of the evolution code (gitignored; skipped when absent)
    tab = joinpath(@__DIR__, "..", "src", "rate_sve", "Imag_ergC_pos_211.npz")
    if isfile(tab)
        using NPZ
        A = npzread(tab)
        alphas = sort(unique(A[:, :, 1])); spins = sort(unique(A[:, :, 2]))
        @testset "Im ω vs the production table (211)" begin
            for (α0, a0) in ((0.1, 0.9), (0.2, 0.7), (0.3, 0.95))
                α = alphas[findmin(abs.(alphas .- α0))[2]]; a = spins[findmin(abs.(spins .- a0))[2]]
                t = reshape(A[:, :, 3], length(alphas), length(spins))[findfirst(==(α), alphas), findfirst(==(a), spins)]
                md, _ = kerr_mode(2, 1, 1, α, a)
                @test imag(md.ω) ≈ t rtol = 1e-5
            end
        end
    else
        @info "Imag_ergC_pos_211.npz not present; skipping the table comparison"
    end
end
