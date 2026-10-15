# Diagnose near-null modes of the Lorenz-gauge operator for a dipole channel.
include("exchange.jl")
const SRC = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(SRC, "Core/constants.jl")); include(joinpath(SRC, "solve_sr_rates.jl"))
using Arpack, Printf
LinearAlgebra.BLAS.set_num_threads(2)
alph = parse(Float64, ARGS[1]); kap0 = parse(Float64, ARGS[2]); kp = parse(Int, ARGS[3])
a = 0.95; mu = alph / GNew
g1 = find_im_part(mu, 1.0, a, 3, 2, 2; Ntot_force=5000, for_s_rates=true, return_both=true)
g3 = find_im_part(mu, 1.0, a, 2, 1, 1; Ntot_force=5000, for_s_rates=true, return_both=true)
md1, _ = kerr_mode(3, 2, 2, alph, a, ComplexF64(Float64(g1[1]) + im * Float64(g1[2])))
md3, _ = kerr_mode(2, 1, 1, alph, a, ComplexF64(Float64(g3[1]) + im * Float64(g3[2])))
es = setup(md1, md3; Lmax=9, drho=0.03)
kap = r -> kap0 / (1 + abs(r))^kp
X, A, b = solve_h(es; kap=kap)
n = NC(es.ab); N = length(es.rg.r)
println("graviton w = ", es.wg, "  m = ", es.mg)
# smallest eigenvalues (shift-invert around 0)
vals, vecs = eigs(A, nev=6, sigma=0.0, which=:LM)
for k in 1:length(vals)
    V = reshape(vecs[:, k], n, N)
    # angular composition: norm per (k,l)
    comp = [sqrt(sum(abs2, V[u, :])) for u in 1:n]
    top = sortperm(comp, rev=true)[1:4]
    rad = [sqrt(sum(abs2, V[:, i])) for i in 1:N]
    ipk = argmax(rad)
    @printf("eig %d: lambda = %.3e %+.3ei   peak at r = %s;  top comps: %s\n", k, real(vals[k]), imag(vals[k]),
            string(round(es.rg.r[ipk], sigdigits=3)), join([string(es.ab.idx[u], ":", round(comp[u] / comp[top[1]], sigdigits=2)) for u in top], "  "))
end
# Z residual and S4 profile of the solution
S4, Zrel = mode4_source(es, X)
for i in [1, 2, 3, 5, 10, 20, 40, 60, 80, 100, 120, 140, 160, 180, es.ireal]
    i > es.ireal && continue
    @printf("  r=%9.3f  Zrel=%.2e  |S4|max=%.3e  |h|=%.3e\n", real(es.rg.r[i]), Zrel[i], maximum(abs.(S4[i, :])), maximum(abs.(X[:, i])))
end
