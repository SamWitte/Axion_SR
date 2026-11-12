# Manufactured-solution test of the 2D solver: pick smooth coefficient profiles c_kl(r) (decaying along the complex
# contour), evaluate E[h] pointwise with exact derivatives, project, solve, compare.
include("solver.jl")
using ForwardDiff, Printf
LinearAlgebra.BLAS.set_num_threads(2)
const FD = ForwardDiff.derivative

function run_test(a, w, m, Lmax; drho=0.03, ord=6, rend=400.0, rcs=60.0, Ltest=4)
    ab = AngBasis(m, Lmax)
    rg = RadGrid(a, w; rcs=rcs, rend=rend, drho=drho, ord=ord)
    n = NC(ab); N = length(rg.r)
    # test profiles: only (k,l) with l <= Ltest
    prof(u, r) = (k, l) = ab.idx[u]
    fprof(u, r) = begin
        k, l = ab.idx[u]
        l > Ltest && return zero(r)
        (1 + 0.1im * k + 0.05l) * r^2 / (r^2 + 7) * exp(-(r - 2) / (6 + k)) * (1 + 0.3k * 5 / (r + 5))
    end
    Xt = zeros(ComplexF64, n, N)
    for i in 1:N, u in 1:n
        Xt[u, i] = fprof(u, rg.r[i])
    end
    # rhs: pointwise E[h_test] from exact jets
    Hj = jet_arrays(w, m)
    b = zeros(ComplexF64, n, N)
    for i in 1:N-1
        r = rg.r[i]
        c0 = [fprof(u, r) for u in 1:n]
        c1 = [FD(t -> fprof(u, r + t), 0.0) for u in 1:n]
        c2 = [FD(t -> FD(t2 -> fprof(u, r + t + t2), 0.0), 0.0) for u in 1:n]
        Tc = zeros(ComplexF64, 10, length(ab.T.th))
        for (j, th) in enumerate(ab.T.th)
            bg = bg_point(a, r, th)
            D1, D2 = cov_derivs(bg, Hj...)
            E, _, _ = einstein_jets(bg, D1, D2)
            v = ab.B[j] * c0; vt = ab.Bt[j] * c0; vtt = ab.Btt[j] * c0
            vr = ab.B[j] * c1; vrt = ab.Bt[j] * c1; vrr = ab.B[j] * c2
            J = zeros(ComplexF64, NJ)
            for cc in 1:10
                J[jet(cc, 1)] = v[cc]; J[jet(cc, 2)] = vr[cc]; J[jet(cc, 3)] = vt[cc]
                J[jet(cc, 4)] = vrr[cc]; J[jet(cc, 5)] = vrt[cc]; J[jet(cc, 6)] = vtt[cc]
            end
            for (cp, (μ, ν)) in enumerate(PAIRS)
                Tc[cp, j] = sum(E[μ, ν, :] .* J)
            end
        end
        b[:, i] = project_tensor(ab, Tc)
    end
    b[:, N] = Xt[:, N]
    t0 = time()
    A = assemble(a, w, m, ab, rg)
    t1 = time()
    X = reshape(A \ vec(b), n, N)
    t2 = time()
    err = maximum(abs.(X .- Xt)) / maximum(abs.(Xt))
    R = reshape(A * vec(Xt) .- vec(b), n, N); rr = [maximum(abs.(R[:, i])) for i in 1:N]; bb = [maximum(abs.(b[:, i])) for i in 1:N]
    for i in [1, 2, 3, 5, 10, 20, N ÷ 2, N - 20, N - 3, N - 2, N - 1]; @printf("   i=%d r=%s  |res|=%.2e  |b|=%.2e  |X-Xt|=%.2e |Xt|=%.2e\n", i, string(round(rg.r[i], sigdigits=4)), rr[i], bb[i], maximum(abs.(X[:, i] .- Xt[:, i])), maximum(abs.(Xt[:, i]))); end
    resid = maximum(abs.(A * vec(Xt) .- vec(b))) / maximum(abs.(b))
    @printf("a=%.2f w=%.3f m=%d Lmax=%d N=%d n=%d nnz=%.2e  assemble %.1fs solve %.1fs   max rel err %.3e  (residual of exact: %.2e)\n",
            a, w, m, Lmax, N, n, nnz(A), t1 - t0, t2 - t1, err, resid)
    return err
end

println("Q * Qinv check: ", maximum(abs.(Qmat(0.7) * Qinv3(0.7)[1] - I)))
run_test(0.0, 0.05, -1, 6; drho=0.06)



run_test(0.0, 0.05, -1, 6; drho=0.03)
run_test(0.7, 0.05, -1, 6; drho=0.03)
run_test(0.7, 0.05, -1, 10; drho=0.03)
