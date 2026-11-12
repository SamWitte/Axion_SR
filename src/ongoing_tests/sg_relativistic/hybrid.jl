# Black-hole-frame ("pinned") formulation for dipole transitions (|m_g| = 1).
#
# In Lorenz gauge with decaying boundary conditions the hole is free (centre-of-mass frame): its displacement b is a
# near-zero mode of the operator (eigenvalue ~ w^2), and near the horizon the corresponding gauge field ~ b M/r^2 is
# ~alpha^-4 larger than the cloud's own field, so discretisation errors there are amplified by ~1/w^2.  We therefore
# solve for
#     h_hyb = h_CM - beta L_zeta g,      zeta = W(r) ∇Φ,
# where Φ e^{-i w v + i m phi} is the l = 1 massless scalar solution regular at the horizon (so ∇Φ is an exact
# Lorenz-compatible "translation"), W = 1 near the hole and 0 beyond the cloud.  h_hyb is pinned (no l = 1 horizon
# displacement) and satisfies  E[h_hyb] + beta E[L_zeta g] = 8π T,  with E[L_zeta g] supported where W' != 0.
# Gauge invariance: zeta decays, so Inf amplitudes are unchanged; near the hole h_hyb is the hole's rest frame.

using OrdinaryDiffEq

# massless ingoing-smooth radial function: Δ f'' + (Δ' - 2iK) f' - (2i w r + lam) f = 0, K = (r^2+a^2) w - a m
function massless_f(a, w, m, l; r1, r2)
    rp = rplus(a)
    lam = l * (l + 1) + a^2 * w^2 - 2a * m * w          # spheroidal eigenvalue ~ l(l+1) for c = a w -> 0
    Del(r) = r^2 - 2r + a^2; dDel(r) = 2r - 2; K(r) = (r^2 + a^2) * w - a * m
    rhs(u, p, r) = [u[2], ((2im * w * r + lam) * u[1] - (dDel(r) - 2im * K(r)) * u[2]) / Del(r)]
    e = 1e-6
    f0 = 1.0 + 0im; f1 = (2im * w * rp + lam) / (dDel(rp) - 2im * K(rp)) * f0
    u0 = [f0 + e * f1, f1]
    sol = solve(ODEProblem(rhs, u0, (rp + e, r2)), Vern9(); reltol=1e-13, abstol=1e-30)
    return sol, lam
end

# Chebyshev fit on [A, B] of a sampled function (for AD-differentiable evaluation)
struct Cheb
    A::Float64; B::Float64; c::Vector{ComplexF64}
end
function Cheb(f, A, B, n)
    k = 0:n-1
    x = cos.(pi .* (k .+ 0.5) ./ n)
    fv = [f((B - A) / 2 * xi + (B + A) / 2) for xi in x]
    c = [2 / n * sum(fv .* cos.(j .* pi .* (k .+ 0.5) ./ n)) for j in 0:n-1]
    c[1] /= 2
    return Cheb(A, B, c)
end
function (ch::Cheb)(r)
    x = (2r - ch.A - ch.B) / (ch.B - ch.A)
    b1 = zero(x) * 0im; b2 = zero(x) * 0im
    for j in length(ch.c):-1:2
        b1, b2 = 2x * b1 - b2 + ch.c[j], b1
    end
    return x * b1 - b2 + ch.c[1]
end

# C^inf window: 1 for r <= R1, 0 for r >= R2
function window(r, R1, R2)
    y = (r - R1) / (R2 - R1)
    bump(t) = t <= 0 ? zero(t) : exp(-1 / t)
    y <= 0 && return one(y)
    y >= 1 && return zero(y)
    b0 = bump(1 - y); b1 = bump(y)
    return b0 / (b0 + b1)
end

struct HybridGauge
    a::Float64; w::Float64; m::Int
    R1::Float64; R2::Float64
    fch::Cheb
    Ylm::Function
end

function HybridGauge(a, w, m; R1, R2)
    sol, _ = massless_f(a, w, m, 1; r1=R1, r2=R2 * 1.05)
    fch = Cheb(r -> sol(r)[1], R1 * 0.95, R2 * 1.05, 80)
    # Y_{1m}(θ) (spheroidal corrections O((a w)^2) neglected)
    Y(th) = m == 0 ? sqrt(3 / (4pi)) * cos(th) : -sign(m) * sqrt(3 / (8pi)) * sin(th)
    return HybridGauge(a, w, m, R1, R2, fch, Y)
end

# zeta^μ(r, θ) = W(r) g^{μν} ∂_ν Φ   (common factor e^{-i w v + i m φ} dropped)
function zeta_up(hg::HybridGauge, r, th)
    a = hg.a
    F(rr, tt) = hg.fch(rr) * hg.Ylm(tt)
    dPhi = [-im * hg.w * F(r, th), ForwardDiff.derivative(t -> F(r + t, th), 0.0), ForwardDiff.derivative(u -> F(r, th + u), 0.0), im * hg.m * F(r, th)]
    G = [gcomp(i, j, a, r, th) for i in 1:4, j in 1:4]
    return window(r, hg.R1, hg.R2) .* (inv(G) * dPhi)
end

# h = L_zeta g, components (i, j) as a function of (r, θ)
function hzeta(hg::HybridGauge, i, j, r, th)
    a = hg.a; D = ForwardDiff.derivative
    z = zeta_up(hg, r, th)
    dz(k) = k == 1 ? -im * hg.w .* z : k == 4 ? im * hg.m .* z : k == 2 ? D(t -> zeta_up(hg, r + t, th), 0.0) : D(u -> zeta_up(hg, r, th + u), 0.0)
    dg(k) = k == 2 ? D(t -> gcomp(i, j, a, r + t, th), 0.0) : k == 3 ? D(u -> gcomp(i, j, a, r, th + u), 0.0) : zero(r)
    dzi = dz(i); dzj = dz(j)
    return sum(z[l] * dg(l) + gcomp(l, j, a, r, th) * dzi[l] + gcomp(i, l, a, r, th) * dzj[l] for l in 1:4)
end

"""Projected E[L_zeta g] (NC x N), nonzero only in the window region."""
function zeta_source(es, hg::HybridGauge)
    ab = es.ab; rg = es.rg; n = NC(ab); N = length(rg.r)
    Hj = jet_arrays(es.wg, es.mg)
    s = zeros(ComplexF64, n, N)
    D = ForwardDiff.derivative
    for i in 1:N
        r = real(rg.r[i])
        (imag(rg.r[i]) != 0 || r <= hg.R1 || r >= hg.R2) && continue
        Tc = zeros(ComplexF64, 10, length(ab.T.th))
        for (j, th) in enumerate(ab.T.th)
            bg = bg_point(es.a, r, th)
            D1, D2 = cov_derivs(bg, Hj...)
            E, _, _ = einstein_jets(bg, D1, D2)
            J = zeros(ComplexF64, NJ)
            for (c, (p, q)) in enumerate(PAIRS)
                f(t, u) = hzeta(hg, p, q, r + t, th + u)
                J[jet(c, 1)] = f(0.0, 0.0); J[jet(c, 2)] = D(t -> f(t, 0.0), 0.0); J[jet(c, 3)] = D(u -> f(0.0, u), 0.0)
                J[jet(c, 4)] = D(t -> D(t2 -> f(t + t2, 0.0), 0.0), 0.0); J[jet(c, 5)] = D(t -> D(u -> f(t, u), 0.0), 0.0)
                J[jet(c, 6)] = D(u -> D(u2 -> f(0.0, u + u2), 0.0), 0.0)
            end
            for (cp, (μ, ν)) in enumerate(PAIRS)
                Tc[cp, j] = sum(E[μ, ν, :] .* J)
            end
        end
        s[:, i] = project_tensor(ab, Tc)
    end
    return s
end

"""Row functional: l = 1 (Y_{1m}) projection of h^{rr} at the horizon (i = 1) -> coefficient vector over X[:, 1]."""
function pin_row(es)
    ab = es.ab; n = NC(ab)
    T = ab.T
    row = zeros(ComplexF64, n)
    for (j, th) in enumerate(T.th)
        bg = bg_point(es.a, rplus(es.a), th)
        # h^{rr} = g^{rα} g^{rβ} h_{αβ};  h_{αβ} = B[j] * c
        wts = [bg.gi[2, α] * bg.gi[2, β] for α in 1:4, β in 1:4]
        coeff = zeros(ComplexF64, 10)
        for α in 1:4, β in 1:4
            coeff[CIDX[α, β]] += wts[α, β]
        end
        row .+= 2pi * T.w[j] * T.Y[(0, 1)][j] .* vec(transpose(coeff) * ab.B[j])
    end
    return row
end

"""Solve the bordered system for (h_hyb, beta)."""
function solve_h_hybrid(es, hg::HybridGauge; verbose=true)
    t0 = time()
    b = rhs_T31(es)
    A = assemble(es.a, es.wg, es.mg, es.ab, es.rg)
    t1 = time()
    s = zeta_source(es, hg)
    t2 = time()
    n = NC(es.ab); N = length(es.rg.r)
    pr = zeros(ComplexF64, n * N); pr[1:n] = pin_row(es)
    scale = maximum(abs, A) / maximum(abs, pr)
    Ab = [A sparse(vec(s)); sparse(transpose(pr .* scale)) 0]
    sol = Ab \ [vec(b); 0]
    X = reshape(sol[1:end-1], n, N); beta = sol[end]
    t3 = time()
    verbose && println("  assemble $(round(t1-t0,digits=1))s  zeta-source $(round(t2-t1,digits=1))s  solve $(round(t3-t2,digits=1))s  beta = $beta")
    return X, beta
end
