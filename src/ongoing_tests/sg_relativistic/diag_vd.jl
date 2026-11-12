# Compare the discrete near-null mode v_d (dominant part of A^{-1} e for near-horizon forcing) with the analytic
# displacement gauge field L_{∇Φ} g (Φ: regular l = 1 massless solution at the graviton frequency).
include("exchange.jl"); include("hybrid.jl")
const SRC = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(SRC, "Core/constants.jl")); include(joinpath(SRC, "solve_sr_rates.jl"))
using Printf, Random
LinearAlgebra.BLAS.set_num_threads(2)
a = parse(Float64, ARGS[1]); alph = 0.2; mu = alph / GNew
g1 = find_im_part(mu, 1.0, 0.95, 3, 2, 2; Ntot_force=5000, for_s_rates=true, return_both=true)
g3 = find_im_part(mu, 1.0, 0.95, 2, 1, 1; Ntot_force=5000, for_s_rates=true, return_both=true)
md1, _ = kerr_mode(3, 2, 2, alph, a, ComplexF64(Float64(g1[1]) + im * Float64(g1[2])))
md3, _ = kerr_mode(2, 1, 1, alph, a, ComplexF64(Float64(g3[1]) + im * Float64(g3[2])))
es = setup(md1, md3; Lmax=5, drho=0.04)
A = assemble(es.a, es.wg, es.mg, es.ab, es.rg)
F = lu(A)
n = NC(es.ab); N = length(es.rg.r)
Random.seed!(1)
e = zeros(ComplexF64, n, N); for i in 1:6; e[:, i] = randn(ComplexF64, n); end
x = reshape(F \ vec(e), n, N)
x = reshape(F \ vec(x), n, N)      # one more inverse iteration
# analytic L_{∇Φ} g (no window: R1 huge) projected on the basis at real radii
rc = real(es.rg.r[es.ireal])
hg = HybridGauge(es.a, es.wg, es.mg; R1=10 * rc, R2=20 * rc)
# need f on [r+, ...]: rebuild Chebyshev on a near range instead
sol, _ = massless_f(es.a, es.wg, es.mg, 1; r1=1.0, r2=rc)
rp = rplus(es.a)
hgn = HybridGauge(es.a, es.wg, es.mg, 1e9, 2e9, Cheb(r -> sol(max(r, rp + 1.1e-6))[1], rp, rc, 120), hg.Ylm)
G = zeros(ComplexF64, n, es.ireal)
for i in 1:es.ireal
    r = real(es.rg.r[i]); i == 1 && (r += 1e-5)
    Tc = zeros(ComplexF64, 10, length(es.ab.T.th))
    for (j, th) in enumerate(es.ab.T.th), (c, (p, q)) in enumerate(PAIRS)
        Tc[c, j] = hzeta(hgn, p, q, r, th)
    end
    G[:, i] = project_tensor(es.ab, Tc)
end
for (lo, hi) in ((1, 20), (1, 60), (1, es.ireal ÷ 2), (es.ireal ÷ 2, es.ireal))
    xs = vec(x[:, lo:hi]); gs = vec(G[:, lo:hi])
    c = dot(gs, xs) / dot(gs, gs)
    @printf("r in [%.2f, %.1f]:  |x - c g| / |x| = %.3e   (c = %.3e)\n", real(es.rg.r[lo]), real(es.rg.r[hi]), norm(xs - c * gs) / norm(xs), abs(c))
end
