# Which near-horizon l = 1 structure does A^{-1} amplify?  Apply A^{-1} to localised random forcing and split the
# l = 1 part of the response near the hole into even (gradient) and odd (curl) vector parts.
include("exchange.jl")
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
idx = es.ab.idx
u(k, l) = findfirst(==((k, l)), idx)
for trial in 1:3
    e = zeros(ComplexF64, n, N)
    for i in 1:6
        e[:, i] = randn(ComplexF64, n)
    end
    x = reshape(F \ vec(e), n, N)
    # l = 1 parts at radial points near the hole
    for i in (1, 10, 40)
        c = x[:, i]
        ev = abs(c[u(5, 1)] - c[u(6, 1)]) + abs(c[u(7, 1)] - c[u(8, 1)])
        od = abs(c[u(5, 1)] + c[u(6, 1)]) + abs(c[u(7, 1)] + c[u(8, 1)])
        sc = abs(c[u(1, 1)]) + abs(c[u(2, 1)]) + abs(c[u(3, 1)]) + abs(c[u(4, 1)]) / real(es.rg.r[i])^2
        rest = norm(c[[k for k in 1:n if idx[k][2] != 1]])
        @printf("trial %d  r=%7.2f  |x|=%.2e   l=1 scalar %.2e  vector even %.2e  odd %.2e   l!=1 %.2e\n", trial, real(es.rg.r[i]),
                norm(c), sc, ev, od, rest)
    end
end
