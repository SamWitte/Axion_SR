# 2D (r, theta) frequency-domain solver for the Lorenz-gauge linearised Einstein equations on Kerr,
#     E_μν[h] = 8π T_μν,       h, T ∝ e^{-i w v + i m phi}       (ingoing Kerr coordinates, G = M = 1)
# Angular: the 10 components are recombined into spin-weighted quantities (s = 0,0,0,0,±1,±1,±2) and expanded in
# sY_lm, l <= Lmax; equations are Galerkin-projected on the same functions (Gauss-Legendre quadrature).
# Radial: finite differences of order `ord` on a uniform grid in rho, r(rho) = r+ + delta (e^rho - 1)(1 + i z sig(rho))
# with exterior complex scaling (sig: smooth step at rho_c) so that outgoing waves decay; h = 0 at the last point.
# At the horizon the equation itself is imposed (regularity in ingoing coordinates selects the ingoing solution).

using SparseArrays, LinearAlgebra
include("swsh.jl"); include("kerr.jl"); include("jets.jl")

const SPIN = [0, 0, 0, 0, 1, -1, 1, -1, 2, -2]

# spin-weighted combinations k = Q[k,c] X_c  and inverse X_c = Qi[c,k] η_k  (c: vv vr vth vph rr rth rph thth thph phph)
function Qmat(th)
    s = sin(th); Q = zeros(ComplexF64, 10, 10)
    Q[1, 1] = 1; Q[2, 2] = 1; Q[3, 5] = 1
    Q[4, 8] = 1; Q[4, 10] = 1 / s^2
    Q[5, 3] = 1; Q[5, 4] = im / s;  Q[6, 3] = 1; Q[6, 4] = -im / s
    Q[7, 6] = 1; Q[7, 7] = im / s;  Q[8, 6] = 1; Q[8, 7] = -im / s
    Q[9, 8] = 1; Q[9, 10] = -1 / s^2; Q[9, 9] = 2im / s
    Q[10, 8] = 1; Q[10, 10] = -1 / s^2; Q[10, 9] = -2im / s
    return Q
end
# Qi and its first two theta derivatives (exact)
function Qinv3(th)
    s, c = sin(th), cos(th)
    S1 = (s, c, -s)                 # sin and derivatives
    S2 = (s^2, 2s * c, 2(c^2 - s^2))  # sin^2 and derivatives
    out = [zeros(ComplexF64, 10, 10) for _ in 1:3]
    for d in 1:3
        Qi = out[d]
        d == 1 && (Qi[1, 1] = 1; Qi[2, 2] = 1; Qi[5, 3] = 1)
        if d == 1
            Qi[3, 5] = 0.5; Qi[3, 6] = 0.5
            Qi[6, 7] = 0.5; Qi[6, 8] = 0.5
            Qi[8, 4] = 0.5; Qi[8, 9] = 0.25; Qi[8, 10] = 0.25
        end
        Qi[4, 5] = S1[d] / (2im); Qi[4, 6] = -S1[d] / (2im)
        Qi[7, 7] = S1[d] / (2im); Qi[7, 8] = -S1[d] / (2im)
        Qi[10, 4] = S2[d] / 2; Qi[10, 9] = -S2[d] / 4; Qi[10, 10] = -S2[d] / 4
        Qi[9, 9] = S1[d] / (4im); Qi[9, 10] = -S1[d] / (4im)
    end
    return out
end

struct AngBasis
    m::Int
    Lmax::Int
    T::SWSHTable
    idx::Vector{Tuple{Int,Int}}           # (k, l) for each angular unknown
    B::Vector{Matrix{ComplexF64}}         # per node: 10 x NC values of coordinate comps
    Bt::Vector{Matrix{ComplexF64}}        # d_theta
    Btt::Vector{Matrix{ComplexF64}}       # d_theta^2
    P::Vector{Matrix{ComplexF64}}         # per node: NC x 10 projection (incl. 2π w_j)
end
NC(ab::AngBasis) = length(ab.idx)

function AngBasis(m, Lmax; nth=Lmax + 16)
    T = SWSHTable(m, Lmax, nth)
    idx = [(k, l) for k in 1:10 for l in lmin_of(SPIN[k], m):Lmax]
    n = length(idx)
    B = Matrix{ComplexF64}[]; Bt = similar(B); Btt = similar(B); P = similar(B)
    for (j, th) in enumerate(T.th)
        Qi, dQi, d2Qi = Qinv3(th); Q = Qmat(th)
        b = zeros(ComplexF64, 10, n); bt = similar(b); btt = similar(b); p = zeros(ComplexF64, n, 10)
        for (u, (k, l)) in enumerate(idx)
            s = SPIN[k]
            Y = T.Y[(s, l)][j]; dY = T.dY[(s, l)][j]; d2Y = T.d2Y[(s, l)][j]
            b[:, u] = Qi[:, k] .* Y
            bt[:, u] = dQi[:, k] .* Y .+ Qi[:, k] .* dY
            btt[:, u] = d2Qi[:, k] .* Y .+ 2 .* dQi[:, k] .* dY .+ Qi[:, k] .* d2Y
            p[u, :] = (2pi * T.w[j] * Y) .* Q[k, :]
        end
        push!(B, b); push!(Bt, bt); push!(Btt, btt); push!(P, p)
    end
    return AngBasis(m, Lmax, T, idx, B, Bt, Btt, P)
end

# ------------------------------------------------------------------ radial grid
"""Fornberg finite-difference weights for derivatives 0..mder at x0 on nodes xs."""
function fdweights(x0, xs, mder)
    n = length(xs); c = zeros(n, mder + 1)
    c1 = 1.0; c4 = xs[1] - x0; c[1, 1] = 1.0
    for i in 2:n
        mn = min(i, mder + 1); c2 = 1.0; c5 = c4; c4 = xs[i] - x0
        for j in 1:i-1
            c3 = xs[i] - xs[j]; c2 *= c3
            if j == i - 1
                for k in mn:-1:2
                    c[i, k] = c1 * ((k - 1) * c[i-1, k-1] - c5 * c[i-1, k]) / c2
                end
                c[i, 1] = -c1 * c5 * c[i-1, 1] / c2
            end
            for k in mn:-1:2
                c[j, k] = (c4 * c[j, k] - (k - 1) * c[j, k-1]) / c3
            end
            c[j, 1] = c4 * c[j, 1] / c3
        end
        c1 = c2
    end
    return c
end

struct RadGrid
    rho::Vector{Float64}
    r::Vector{ComplexF64}
    dr::Vector{ComplexF64}      # dr/drho
    d2r::Vector{ComplexF64}
    D1::SparseMatrixCSC{Float64,Int}   # d/drho
    D2::SparseMatrixCSC{Float64,Int}   # d^2/drho^2
end

"""
    RadGrid(a, w; delta=1.0, rcs, rend, drho=0.03, z=1.0, wid=0.5, ord=6)
Uniform rho grid from r+ to rend; complex scaling switched on smoothly (C^inf, exactly zero below) between
r = rcs and rho_cs + wid, Im r / Re r -> zeta sign(w).
"""
function RadGrid(a, w; delta=1.0, rcs, rend, drho=0.03, zeta=1.0, wid=1.5, ord=6)
    rp = rplus(a)
    rho1 = log(rend / delta + 1)
    N = Int(ceil(rho1 / drho)) + 1
    rho = collect(range(0, rho1, length=N))
    rhoc = log(rcs / delta + 1)
    z = zeta * sign(w)
    r = zeros(ComplexF64, N); dr = similar(r); d2r = similar(r)
    # C^infinity step, exactly 0 for rho <= rhoc and 1 for rho >= rhoc + wid
    bump(x) = x <= 0 ? 0.0 : exp(-1 / x)
    function step3(x)
        h = 1e-4
        st(y) = (b0 = bump(y); b1 = bump(1 - y); b0 / (b0 + b1))
        y = (x - rhoc) / wid
        y <= 0 && return (0.0, 0.0, 0.0)
        y >= 1 && return (1.0, 0.0, 0.0)
        d1 = ForwardDiff.derivative(st, y) / wid
        d2 = ForwardDiff.derivative(t -> ForwardDiff.derivative(st, t), y) / wid^2
        return (st(y), d1, d2)
    end
    for (i, x) in enumerate(rho)
        sg, sg1, sg2 = step3(x)
        u = delta * (exp(x) - 1); u1 = delta * exp(x); u2 = u1
        q = 1 + im * z * sg; q1 = im * z * sg1; q2 = im * z * sg2
        r[i] = rp + u * q; dr[i] = u1 * q + u * q1; d2r[i] = u2 * q + 2u1 * q1 + u * q2
    end
    # FD matrices on the uniform rho grid (one-sided near the ends)
    h = rho[2] - rho[1]
    I1 = Int[]; J1 = Int[]; V1 = Float64[]; V2 = Float64[]
    half = ord ÷ 2
    for i in 1:N
        lo = clamp(i - half, 1, N - ord); hi = lo + ord
        nodes = collect(lo:hi)
        W = fdweights(rho[i], rho[nodes], 2)
        for (q, jn) in enumerate(nodes)
            push!(I1, i); push!(J1, jn); push!(V1, W[q, 2]); push!(V2, W[q, 3])
        end
    end
    D1 = sparse(I1, J1, V1, N, N); D2 = sparse(I1, J1, V2, N, N)
    return RadGrid(rho, r, dr, d2r, D1, D2)
end

# ------------------------------------------------------------------ assembly
"""Angular blocks (A0, A1, A2) at radius r: equation = A0 c + A1 ∂_r c + A2 ∂_r² c (projected)."""
function angular_blocks(a, r, w, m, ab::AngBasis, Hj; kap=0.0)
    H, dH, ddH = Hj
    n = NC(ab)
    A0 = zeros(ComplexF64, n, n); A1 = zeros(ComplexF64, n, n); A2 = zeros(ComplexF64, n, n)
    C = zeros(ComplexF64, 10, NJ)
    for (j, th) in enumerate(ab.T.th)
        bg = bg_point(a, r, th)
        D1, D2 = cov_derivs(bg, H, dH, ddH)
        E, _, _ = einstein_jets(bg, D1, D2; kap=kap(r), nlow=ks_normal(a, r, th))
        for (cp, (μ, ν)) in enumerate(PAIRS)
            C[cp, :] = E[μ, ν, :]
        end
        cv(d) = C[:, [jet(c, d) for c in 1:10]]
        M0 = cv(1) * ab.B[j] .+ cv(3) * ab.Bt[j] .+ cv(6) * ab.Btt[j]
        M1 = cv(2) * ab.B[j] .+ cv(5) * ab.Bt[j]
        M2 = cv(4) * ab.B[j]
        A0 .+= ab.P[j] * M0; A1 .+= ab.P[j] * M1; A2 .+= ab.P[j] * M2
    end
    return A0, A1, A2
end

"""Assemble the global sparse operator (Dirichlet at the outer end)."""
function assemble(a, w, m, ab::AngBasis, rg::RadGrid; drop=0.0, kap=r -> 0.0)
    n = NC(ab); N = length(rg.r)
    Hj = jet_arrays(w, m)
    I = Int[]; J = Int[]; V = ComplexF64[]
    for i in 1:N
        if i == N
            for u in 1:n
                push!(I, (i - 1) * n + u); push!(J, (i - 1) * n + u); push!(V, 1.0)
            end
            continue
        end
        A0, A1, A2 = angular_blocks(a, rg.r[i], w, m, ab, Hj; kap=kap)
        rp1 = rg.dr[i]; rp2 = rg.d2r[i]
        Bρ = A1 ./ rp1 .- A2 .* (rp2 / rp1^3)     # coefficient of ∂_rho c
        Bρρ = A2 ./ rp1^2                         # coefficient of ∂_rho^2 c
        tol = drop * max(maximum(abs, A0), maximum(abs, Bρ), maximum(abs, Bρρ))
        cols1, vals1 = findnz(rg.D1[i, :]); cols2, vals2 = findnz(rg.D2[i, :])
        cols = union(cols1, cols2, [i])
        for jn in cols
            blk = zeros(ComplexF64, n, n)
            jn == i && (blk .+= A0)
            k1 = findfirst(==(jn), cols1); k1 === nothing || (blk .+= vals1[k1] .* Bρ)
            k2 = findfirst(==(jn), cols2); k2 === nothing || (blk .+= vals2[k2] .* Bρρ)
            for v in 1:n, u in 1:n
                x = blk[u, v]
                abs(x) > tol || continue
                push!(I, (i - 1) * n + u); push!(J, (jn - 1) * n + v); push!(V, x)
            end
        end
    end
    return sparse(I, J, V, N * n, N * n)
end

"""Project a tensor field given at the angular nodes (Tc[c, j] coordinate comps) onto the basis -> NC vector."""
function project_tensor(ab::AngBasis, Tc)
    out = zeros(ComplexF64, NC(ab))
    for j in eachindex(ab.T.th)
        out .+= ab.P[j] * Tc[:, j]
    end
    return out
end

"""Reconstruct jets of h (60 values) at radial index i and node j from the solution coefficients X (n x N)."""
function h_jets(ab::AngBasis, rg::RadGrid, X, i, j)
    c = X[:, i]
    dρ = X * rg.D1[i, :]; dρρ = X * rg.D2[i, :]
    rp1 = rg.dr[i]; rp2 = rg.d2r[i]
    cr = dρ ./ rp1
    crr = dρρ ./ rp1^2 .- dρ .* (rp2 / rp1^3)
    J = zeros(ComplexF64, NJ)
    v = ab.B[j] * c; vt = ab.Bt[j] * c; vtt = ab.Btt[j] * c
    vr = ab.B[j] * cr; vrt = ab.Bt[j] * cr; vrr = ab.B[j] * crr
    for cc in 1:10
        J[jet(cc, 1)] = v[cc]; J[jet(cc, 2)] = vr[cc]; J[jet(cc, 3)] = vt[cc]
        J[jet(cc, 4)] = vrr[cc]; J[jet(cc, 5)] = vrt[cc]; J[jet(cc, 6)] = vtt[cc]
    end
    return J
end
