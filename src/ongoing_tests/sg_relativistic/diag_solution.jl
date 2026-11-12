include("exchange.jl")
const SRC = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(SRC, "Core/constants.jl")); include(joinpath(SRC, "solve_sr_rates.jl"))
using Printf
LinearAlgebra.BLAS.set_num_threads(2)
st(s) = (parse(Int, s[1]), parse(Int, s[2]), parse(Int, s[3]))
s1 = st(ARGS[1]); s3 = st(ARGS[2]); alph = parse(Float64, ARGS[3]); drho = parse(Float64, ARGS[4]); Lx = parse(Int, ARGS[5])
a = 0.95; mu = alph / GNew
g1 = find_im_part(mu, 1.0, a, s1...; Ntot_force=5000, for_s_rates=true, return_both=true)
g3 = find_im_part(mu, 1.0, a, s3...; Ntot_force=5000, for_s_rates=true, return_both=true)
md1, _ = kerr_mode(s1..., alph, a, ComplexF64(Float64(g1[1]) + im * Float64(g1[2])))
md3, _ = kerr_mode(s3..., alph, a, ComplexF64(Float64(g3[1]) + im * Float64(g3[2])))
es = setup(md1, md3; Lmax=s1[2] + s3[2] + Lx, drho=drho)
X, A, b = solve_h(es)
n = NC(es.ab); N = length(es.rg.r)
println("relative residual of linear solve: ", norm(A * vec(X) - vec(b)) / norm(vec(b)))
# dominant components at a few radii, in "dimensionless" form: divide angular-ish comps by r^p
pw = [0, 0, 0, 2, 1, 1, 1, 1, 2, 2]   # h_vv,h_vr,h_rr ~ 1; trace ~ r^2; vA ~ r; rA ~ r; AB ~ r^2
for i in [1, 5, 20, 60, 100, 140, es.ireal, es.ireal + 40, N - 50]
    r = es.rg.r[i]
    v = [abs(X[u, i]) / abs(r)^pw[es.ab.idx[u][1]] for u in 1:n]
    top = sortperm(v, rev=true)[1:5]
    @printf("r=%-22s  %s\n", string(round(r, sigdigits=4)), join([@sprintf("%s:%.2e", string(es.ab.idx[u]), v[u]) for u in top], "  "))
end
# Einstein residual dG[h] - 8πT, relative, at a few radii (pointwise at the angular nodes)
Hj = jet_arrays(es.wg, es.mg)
w1 = real(md1.w); w3 = real(md3.w)
for i in [1, 3, 10, 30, 60, 100, 140]
    r = real(es.rg.r[i]); f1 = radial_f(md1, r); f3 = radial_f(md3, r)
    worst = 0.0; scale = 0.0; zw = 0.0
    for (j, th) in enumerate(es.ab.T.th)
        bg = bg_point(a, r, th)
        D1, D2 = cov_derivs(bg, Hj...)
        E, dG, Z = einstein_jets(bg, D1, D2)
        J = h_jets(es.ab, es.rg, X, i, j)
        d1, _ = scalar_jets(f1..., es.S1[1][j], es.S1[2][j], es.S1[3][j], w1, md1.m)
        d3, _ = scalar_jets(f3..., es.S3[1][j], es.S3[2][j], es.S3[3][j], w3, md3.m)
        d3c = conj.(d3); dot = sum(bg.gi[μ, ν] * d1[μ] * d3c[ν] for μ in 1:4, ν in 1:4)
        phi13 = f1[1] * es.S1[1][j] * conj(f3[1] * es.S3[1][j])
        for μ in 1:4, ν in 1:4
            T = 8pi * (d1[μ] * d3c[ν] + d1[ν] * d3c[μ] - bg.g[μ, ν] * (dot + alph^2 * phi13))
            worst = max(worst, abs(sum(dG[μ, ν, :] .* J) - T)); scale = max(scale, abs(T), abs(sum(dG[μ, ν, :] .* J)))
        end
    end
    @printf("  r=%8.3f  max|dG - 8πT| = %.3e   scale = %.3e\n", r, worst, scale)
end
