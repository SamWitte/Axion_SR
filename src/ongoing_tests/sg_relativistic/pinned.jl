# Black-hole-frame solution for dipole transitions (|m_g| = 1) without near-zero modes.
#
# With decaying boundary conditions the Lorenz-gauge operator A has two near-zero l = 1 modes for |m| = 1: the hole's
# displacement (eigenvalue ~ w^2) and its tilt / transverse angular momentum (~ w).  Any unknown that is fixed only
# through these modes is determined with ~1/w^2 error amplification.  We instead
#   (1) pin the hole: p_d·h = 0 (no l = 1 horizon displacement) and p_t·h = 0 (no odd l = 1 horizon rotation), and let
#       auxiliary "struts" z_d, z_t (a force and a torque on the hole) take up the force balance:
#           A h + tau_d z_d + tau_t z_t = rhs,   p·h = given   (bordered system, well conditioned);
#   (2) solve with rhs = 8πT  -> (h0, tau0);
#   (3) solve for the frame fields q_k = W U_k + q'_k, U_d = uniform acceleration (h_tt = 2 r sinθ e^{imφ}),
#       U_t = rotating frame (h_ti = (Ω × x)_i), with rhs = -δG[W U_k] and pin data -p·(W U_k)  -> (q'_k, tau_k);
#   (4) h = h0 + sum_k a_k q_k with the struts cancelled: tau0 + sum_k a_k tau_k = 0.
# h then satisfies δG[h] = 8πT exactly (no strut), is pinned at the horizon, and decays: it is the physical solution in
# the hole's frame.  The result does not depend on the choice of U_k (any field with nonzero strut works).

function frame_field(kind::Symbol, m, i, j, r, th)
    # components in ingoing Kerr coords (v, r, θ, φ); common factor e^{-i w v + i m φ} dropped; dt = dv - dr
    s = sin(th); c = cos(th)
    i > j && return frame_field(kind, m, j, i, r, th)
    if kind == :d || kind == :d2
        F = kind == :d ? 2 * r * s : 2 * r * s * (1 + 3 / (r + 1)) * (1 + 0.2 * c^2)
        (i, j) == (1, 1) && return F + 0im
        (i, j) == (1, 2) && return -F + 0im
        (i, j) == (2, 2) && return F + 0im
        return zero(r) * 0im
    else   # :t, (Ω×x)·dx for Ω = x̂ ± iŷ:  r^2 (±i dθ - s c dφ) e^{±iφ}
        Ath = sign(m) * im * r^2; Aph = -r^2 * s * c + 0im
        (i, j) == (1, 3) && return Ath
        (i, j) == (1, 4) && return Aph
        (i, j) == (2, 3) && return -Ath
        (i, j) == (2, 4) && return -Aph
        return zero(r) * 0im
    end
end

windowed(kind, m, R1, R2) = (i, j, r, th) -> window(r, R1, R2) * frame_field(kind, m, i, j, r, th)

"""Jets (60) of an analytic tensor field hfun(i, j, r, θ) at a point."""
function analytic_jets(hfun, r, th)
    D = ForwardDiff.derivative
    J = zeros(ComplexF64, NJ)
    for (c, (p, q)) in enumerate(PAIRS)
        f(t, u) = hfun(p, q, r + t, th + u)
        J[jet(c, 1)] = f(0.0, 0.0); J[jet(c, 2)] = D(t -> f(t, 0.0), 0.0); J[jet(c, 3)] = D(u -> f(0.0, u), 0.0)
        J[jet(c, 4)] = D(t -> D(t2 -> f(t + t2, 0.0), 0.0), 0.0); J[jet(c, 5)] = D(t -> D(u -> f(t, u), 0.0), 0.0)
        J[jet(c, 6)] = D(u -> D(u2 -> f(0.0, u + u2), 0.0), 0.0)
    end
    return J
end

"""Projected -δG[hfun] (NC x N) over the real radial points where hfun is nonzero (r < R2)."""
function minus_dG_source(es, hfun, R2)
    ab = es.ab; rg = es.rg; n = NC(ab); N = length(rg.r)
    Hj = jet_arrays(es.wg, es.mg)
    s = zeros(ComplexF64, n, N)
    for i in 1:N
        r = rg.r[i]
        (imag(r) != 0 || real(r) >= R2) && continue
        r = real(r)
        Tc = zeros(ComplexF64, 10, length(ab.T.th))
        for (j, th) in enumerate(ab.T.th)
            bg = bg_point(es.a, r, th)
            D1, D2 = cov_derivs(bg, Hj...)
            _, dG, _ = einstein_jets(bg, D1, D2)
            J = analytic_jets(hfun, r, th)
            for (cp, (μ, ν)) in enumerate(PAIRS)
                Tc[cp, j] = -sum(dG[μ, ν, :] .* J)
            end
        end
        s[:, i] = project_tensor(ab, Tc)
    end
    return s
end

"""Basis coefficients of an analytic tensor field at radial index i."""
function analytic_coeffs(es, hfun, i)
    r = real(es.rg.r[i])
    Tc = zeros(ComplexF64, 10, length(es.ab.T.th))
    for (j, th) in enumerate(es.ab.T.th), (c, (p, q)) in enumerate(PAIRS)
        Tc[c, j] = hfun(p, q, r, th)
    end
    return project_tensor(es.ab, Tc)
end

"""Pin rows (on X[:, 1]) and strut columns (localised at the first few radial points)."""
function pins_and_struts(es; nstrut=3)
    ab = es.ab; n = NC(ab); N = length(es.rg.r)
    u(k, l) = findfirst(==((k, l)), ab.idx)
    pd = zeros(ComplexF64, n * N); pd[1:n] = pin_row(es)
    pt = zeros(ComplexF64, n * N); pt[u(5, 1)] = 1; pt[u(6, 1)] = 1          # odd l = 1 part of h_vA at r+
    zd = zeros(ComplexF64, n * N); zt = zeros(ComplexF64, n * N)
    for i in 1:nstrut
        zd[(i - 1) * n + u(1, 1)] = 1.0          # vv equation, l = 1
        zt[(i - 1) * n + u(5, 1)] = 1.0; zt[(i - 1) * n + u(6, 1)] = 1.0   # odd vA equation, l = 1
    end
    return pd, pt, zd, zt
end

"""
Solve for the hole-frame metric perturbation.  Returns the coefficient matrix X (numerical part), the analytic part as a
function hA(i, j, r, θ), and diagnostics.
"""
function solve_h_pinned(es; R1, R2, verbose=true, odd=:auto, frame_scale=1.0, nstrut=3, dkind=:d)
    # odd = :frame  pin the odd l = 1 mode and cancel its torque with the rotating-frame field (Kerr: tilt is gauge)
    #       :pin    pin it and leave the (small) torque uncompensated (Schwarzschild: transverse dJ is physical, ~0 here)
    #       :none   do not pin it
    odd == :auto && (odd = es.a != 0 ? :frame : :pin)
    t0 = time()
    n = NC(es.ab); N = length(es.rg.r)
    b = rhs_T31(es)
    A = assemble(es.a, es.wg, es.mg, es.ab, es.rg)
    pd, pt, zd, zt = pins_and_struts(es; nstrut=nstrut)
    sc = maximum(abs, A)
    pdn = pd .* (sc / maximum(abs, pd)); ptn = pt .* sc
    P = odd == :none ? reshape(pdn, :, 1) : [pdn ptn]
    Z = odd == :none ? reshape(zd, :, 1) : [zd zt]
    np_ = size(P, 2)
    kinds = odd == :frame ? [dkind, :t] : [dkind]
    Ab = [A sparse(Z); sparse(transpose(P)) spzeros(np_, np_)]
    F = lu(Ab)
    t1 = time()
    fields = [windowed(k, es.mg, R1, R2) for k in kinds]
    xs = []
    x0 = F \ [vec(b); zeros(np_)]
    for f in fields
        sk = minus_dG_source(es, f, R2)
        ck = zeros(ComplexF64, n * N); ck[1:n] = analytic_coeffs(es, f, 1)
        push!(xs, F \ [vec(sk); -(transpose(P) * ck)])
    end
    nk = length(kinds)
    tau(x) = x[n*N+1:n*N+nk]                         # struts to be cancelled (first nk)
    M = reduce(hcat, [tau(x) for x in xs])
    acoef = -(M \ tau(x0)) .* frame_scale
    Xv = x0[1:n*N]; tres = x0[n*N+1:end]
    for k in 1:nk
        Xv = Xv .+ acoef[k] .* xs[k][1:n*N]
        tres = tres .+ acoef[k] .* xs[k][n*N+1:end]
    end
    X = reshape(Xv, n, N)
    hA = (i, j, r, th) -> sum(acoef[k] * fields[k](i, j, r, th) for k in 1:nk)
    t2 = time()
    verbose && println("  pinned (odd=$odd, frame fields $(kinds)): assemble+factor $(round(t1-t0,digits=1))s, frame+solves $(round(t2-t1,digits=1))s;  ",
                       "tau0 = ", round.(x0[n*N+1:end], sigdigits=3), "  residual struts = ", round.(tres, sigdigits=3),
                       "  a = ", round.(acoef, sigdigits=4))
    return X, hA, acoef
end

"""Mode-4 source with the analytic frame part added to the numerical jets."""
function mode4_source_pinned(es::ExchangeSetup, X, hA)
    ab = es.ab; rg = es.rg
    w1 = real(es.md1.w); m1 = es.md1.m
    S1, dS1, d2S1 = es.S1
    Hj = jet_arrays(es.wg, es.mg)
    nth = length(ab.T.th)
    S4 = zeros(ComplexF64, es.ireal, nth)
    for i in 1:es.ireal
        r = real(rg.r[i])
        f1 = radial_f(es.md1, r)
        for (j, th) in enumerate(ab.T.th)
            bg = bg_point(es.a, r, th)
            D1, D2 = cov_derivs(bg, Hj...)
            _, _, Z = einstein_jets(bg, D1, D2)
            J = h_jets(ab, rg, X, i, j) .+ analytic_jets(hA, r, th)
            d1, dd = scalar_jets(f1..., S1[j], dS1[j], d2S1[j], w1, m1)
            K = box_coupling(bg, Hj[1], D1, Z, d1, dd)
            S4[i, j] = sum(K .* J)
        end
    end
    return S4
end
