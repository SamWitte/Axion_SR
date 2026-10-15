# Checks of the Kerr modes: eigenvalue vs find_im_part, |R|^2 vs solve_radial, Klein-Gordon residual in ingoing coords.
include("solver.jl"); include("modes.jl")
const SRC = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(SRC, "Core/constants.jl")); include(joinpath(SRC, "solve_sr_rates.jl"))
using Printf

for (n, l, m, alph) in ((2, 1, 1, 0.1), (3, 2, 2, 0.3), (2, 1, 1, 0.3), (5, 4, 4, 0.5))
    a = 0.95
    mu = alph / GNew
    wR, wI = find_im_part(mu, 1.0, a, n, l, m; Ntot_force=5000, for_s_rates=true, return_both=true)
    w0 = ComplexF64(Float64(wR) + im * Float64(wI))
    t0 = time()
    md, res = kerr_mode(n, l, m, alph, a, w0)
    @printf("(%d%d%d) alpha=%.2f  w(find_im_part) = %.15f %+.3ei   w(refined) = %.15f %+.3ei  |res|=%.1e  K=%d  [%.1fs]\n",
            n, l, m, alph, real(w0), imag(w0), real(md.w), imag(md.w), abs(res), length(md.b) - 1, time() - t0)
    # compare |R|^2 shape with solve_radial
    rl, rf, _ = solve_radial(mu, 1.0, a, n, l, m; rpts=400, return_erg=true, Ntot_safe=5000, use_heunc=true, pre_compute_erg=w0)
    rl = Float64.(real.(rl)); rf = ComplexF64.(rf)
    mine = [abs2(radial_f(md, r)[1]) for r in rl]
    theirs = abs2.(rf)
    sel = theirs .> 1e-6 * maximum(theirs)
    @printf("   max | |R|^2_mine/|R|^2_solve_radial - 1 | over the cloud = %.2e\n", maximum(abs.(mine[sel] ./ theirs[sel] .- 1)))
    # Klein-Gordon residual (□ - mu^2) phi / (mu^2 phi) at a few points, with Re(w)
    Tm = SWSHTable(m, l + 16, 30; smax=0)
    S, dS, d2S = spheroidal_on_nodes(md.sph, Tm)
    wr = md.w
    worst = 0.0
    for r in (rplus(a) + 1e-3, 3.0, n^2 / alph^2, 3n^2 / alph^2), j in (3, 15, 27)
        th = Tm.th[j]
        f, f1, f2 = radial_f(md, r)
        d1 = [-im * wr * f * S[j], f1 * S[j], f * dS[j], im * m * f * S[j]]
        k = [-im * wr, 0, 0, im * m]
        dd = zeros(ComplexF64, 4, 4)
        for i in 1:4, jj in 1:4
            if i in (1, 4) && jj in (1, 4); dd[i, jj] = k[i] * k[jj] * f * S[j]
            elseif i in (1, 4); dd[i, jj] = k[i] * d1[jj]
            elseif jj in (1, 4); dd[i, jj] = k[jj] * d1[i]
            else dd[i, jj] = (i, jj) == (2, 2) ? f2 * S[j] : (i, jj) == (3, 3) ? f * d2S[j] : f1 * dS[j]
            end
        end
        bg = bg_point(a, r, th)
        box = sum(bg.gi[μ, ν] * (dd[μ, ν] - sum(bg.Gam[λ, μ, ν] * d1[λ] for λ in 1:4)) for μ in 1:4, ν in 1:4)
        rel = abs(box - alph^2 * f * S[j]) / abs(alph^2 * f * S[j])
        worst = max(worst, rel)
    end
    @printf("   KG residual max |(□-mu^2)phi| / |mu^2 phi| = %.2e\n", worst)
end
