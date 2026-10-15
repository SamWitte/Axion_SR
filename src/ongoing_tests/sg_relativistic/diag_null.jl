# Near-null modes of the row-scaled operator; characterise each by its Lorenz violation and Einstein residual.
include("exchange.jl")
const SRC = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(SRC, "Core/constants.jl")); include(joinpath(SRC, "solve_sr_rates.jl"))
using Arpack, Printf
LinearAlgebra.BLAS.set_num_threads(2)
alph = 0.2; a = parse(Float64, ARGS[1]); mu = alph / GNew
g1 = find_im_part(mu, 1.0, max(a, 0.5), 3, 2, 2; Ntot_force=5000, for_s_rates=true, return_both=true)
g3 = find_im_part(mu, 1.0, max(a, 0.5), 2, 1, 1; Ntot_force=5000, for_s_rates=true, return_both=true)
md1, _ = kerr_mode(3, 2, 2, alph, a, ComplexF64(Float64(g1[1]) + im * Float64(g1[2])))
md3, _ = kerr_mode(2, 1, 1, alph, a, ComplexF64(Float64(g3[1]) + im * Float64(g3[2])))
es = setup(md1, md3; Lmax=7, drho=0.04)
A = assemble(es.a, es.wg, es.mg, es.ab, es.rg)
n = NC(es.ab); N = length(es.rg.r)
# row scaling: each block row by 1/(max abs entry in that block row)
d = zeros(size(A, 1))
rows = rowvals(A); vals = nonzeros(A)
for j in 1:size(A, 2), k in nzrange(A, j)
    d[rows[k]] = max(d[rows[k]], abs(vals[k]))
end
for i in 1:N
    blk = (i - 1) * n .+ (1:n); d[blk] .= maximum(d[blk])
end
As = spdiagm(0 => 1 ./ d) * A
println("graviton w = ", es.wg, "  N = ", N, "  NC = ", n)
vals_, vecs = eigs(As, nev=8, sigma=0.0, which=:LM)
Hj = jet_arrays(es.wg, es.mg)
for k in sortperm(abs.(vals_))
    V = reshape(vecs[:, k], n, N)
    rad = [sqrt(sum(abs2, V[:, i]) ) for i in 1:N]
    ipk = argmax(rad)
    comp = [sqrt(sum(abs2, V[u, :])) for u in 1:n]; top = sortperm(comp, rev=true)[1:3]
    # Lorenz violation and Einstein-residual relative at peak (real part of grid only)
    i = min(ipk, es.ireal); zr = 0.0; gr = 0.0; hr = 0.0
    for (j, th) in enumerate(es.ab.T.th)
        bg = bg_point(es.a, real(es.rg.r[i]), th)
        D1, D2 = cov_derivs(bg, Hj...); _, dG, Z = einstein_jets(bg, D1, D2)
        J = h_jets(es.ab, es.rg, V, i, j)
        zr = max(zr, maximum(abs(sum(Z[ν, :] .* J)) for ν in 1:4))
        gr = max(gr, maximum(abs(sum(dG[μ, ν, :] .* J)) for μ in 1:4, ν in 1:4))
        hr = max(hr, maximum(abs(sum(D1[ρ, μ, ν, :] .* J)) for ρ in 1:4, μ in 1:4, ν in 1:4))
    end
    @printf("lambda=%.3e  peak r=%-18s  |Z|/|∇h|=%.2e |dG|/|∇h|=%.2e   top %s  (horizon/peak amplitude %.2e)\n", abs(vals_[k]),
            string(round(es.rg.r[ipk], sigdigits=3)), zr / hr, gr / hr, join(string.(es.ab.idx[top]), " "), rad[1] / rad[ipk])
end
