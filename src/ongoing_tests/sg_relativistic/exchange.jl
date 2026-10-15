# Exact one-graviton exchange for N1 x N1 -> N3 x {BH, Inf} on Kerr.
#   1. T31_μν = ∂_μφ1 ∂_νφ3* + ∂_νφ1 ∂_μφ3* - g_μν (∂φ1·∂φ3* + mu^2 φ1 φ3*)      (frequency w1 - w3, m1 - m3)
#   2. E_μν[h] = 8π T31_μν  (Lorenz gauge, ingoing Kerr coords, complex scaling; solver.jl)
#   3. source for mode 4:  (□ - mu^2) φ4 = S4 = -δ□[h] φ1 = h^μν ∇_μ∇_ν φ1 + Z^ν ∂_ν φ1
# In the Newtonian limit S4 -> -4 mu^4 φ1 ∫ φ1 φ3*/|x-y|, so J_GR ≡ S4/(-4 alpha^4) is normalised like the Newtonian
# self-gravity source J_N = φ1 Φ[φ1 φ3*] of sg_vs_si.jl (lengths in GM), and A_GR/A_N = J_GR/J_N after projection
# with the Kerr Green's function of mode 4.

using Interpolations
include("solver.jl"); include("modes.jl")

# scalar-field derivatives (∂_μ, ∂_μ∂_ν) of F(r) S(θ) e^{-i w v + i m φ} from radial and angular jets
function scalar_jets(f, f1, f2, S, dS, d2S, w, m)
    k = (-im * w, 0.0im, 0.0im, im * m)
    d1 = [k[1] * f * S, f1 * S, f * dS, k[4] * f * S]
    dd = zeros(ComplexF64, 4, 4)
    for i in 1:4, j in 1:4
        if i in (1, 4) && j in (1, 4)
            dd[i, j] = k[i] * k[j] * f * S
        elseif i in (1, 4)
            dd[i, j] = k[i] * d1[j]
        elseif j in (1, 4)
            dd[i, j] = k[j] * d1[i]
        else
            dd[i, j] = (i, j) == (2, 2) ? f2 * S : (i, j) == (3, 3) ? f * d2S : f1 * dS
        end
    end
    return d1, dd
end

struct ExchangeSetup
    md1::KerrMode; md3::KerrMode
    a::Float64; alph::Float64
    wg::Float64; mg::Int                 # graviton frequency and azimuthal number
    ab::AngBasis; rg::RadGrid
    S1::NTuple{3,Vector{ComplexF64}}; S3::NTuple{3,Vector{ComplexF64}}
    ireal::Int                           # last radial index with real r and nonzero source region
    delta::Float64
end

function setup(md1, md3; Lmax=nothing, drho=0.03, ord=6, zeta=1.0, delta=1.0, src_tol=1e-15, csfac=1.5, wfac=30.0, wid=1.5, tail=50.0)
    a = md1.a; alph = md1.alph
    wg = real(md1.w) - real(md3.w); mg = md1.m - md3.m
    Lmax = isnothing(Lmax) ? md1.l + md3.l + 6 : Lmax
    ab = AngBasis(mg, Lmax)
    # radial extent: source negligible beyond rcut (|f1 f3| ~ e^{-(k1+k3) r}); complex scaling beyond
    k1 = real(sqrt(alph^2 - real(md1.w)^2 + 0im)); k3 = real(sqrt(alph^2 - real(md3.w)^2 + 0im))
    # the mode-4 source h ∇∇φ1 must be kept until φ1 itself is negligible: in the hole's frame the recoil part grows
    # like r φ1 and is projected with weight ~ r^(l4 + 3), so cut at k1 r = tail rather than where |φ1 φ3| is small
    rcut = max(-log(src_tol) / (k1 + k3), tail / k1)
    rcs = csfac * rcut
    rend = rcs * 3 + wfac / abs(wg)
    rg = RadGrid(a, wg; delta=delta, rcs=rcs, rend=rend, drho=drho, zeta=zeta, ord=ord, wid=wid)
    ireal = findlast(x -> imag(x) == 0 && real(x) <= rcut, rg.r)
    nth = length(ab.T.th)
    T1 = SWSHTable(md1.m, maximum(keys(md1.sph)), nth; smax=0)
    T3 = SWSHTable(md3.m, maximum(keys(md3.sph)), nth; smax=0)
    S1 = spheroidal_on_nodes(md1.sph, T1); S3 = spheroidal_on_nodes(md3.sph, T3)
    return ExchangeSetup(md1, md3, a, alph, wg, mg, ab, rg, S1, S3, ireal, delta)
end

"""Projected right-hand side 8π T31 (NC x N)."""
function rhs_T31(es::ExchangeSetup)
    ab = es.ab; rg = es.rg; n = NC(ab); N = length(rg.r)
    w1 = real(es.md1.w); w3 = real(es.md3.w); m1 = es.md1.m; m3 = es.md3.m; mu2 = es.alph^2
    b = zeros(ComplexF64, n, N)
    S1, dS1, d2S1 = es.S1; S3, dS3, d2S3 = es.S3
    for i in 1:es.ireal
        r = real(rg.r[i])
        f1 = radial_f(es.md1, r); f3 = radial_f(es.md3, r)
        Tc = zeros(ComplexF64, 10, length(ab.T.th))
        for (j, th) in enumerate(ab.T.th)
            bg = bg_point(es.a, r, th)
            d1, _ = scalar_jets(f1..., S1[j], dS1[j], d2S1[j], w1, m1)
            d3, _ = scalar_jets(f3..., S3[j], dS3[j], d2S3[j], w3, m3)
            d3c = conj.(d3)                                   # ∂_μ φ3*  (real frequencies)
            dot = sum(bg.gi[μ, ν] * d1[μ] * d3c[ν] for μ in 1:4, ν in 1:4)
            phi13 = f1[1] * S1[j] * conj(f3[1] * S3[j])
            for (c, (μ, ν)) in enumerate(PAIRS)
                Tc[c, j] = 8pi * (d1[μ] * d3c[ν] + d1[ν] * d3c[μ] - bg.g[μ, ν] * (dot + mu2 * phi13))
            end
        end
        b[:, i] = project_tensor(ab, Tc)
    end
    return b
end

"""Solve for h; returns coefficient matrix X (NC x N) and the residual of the linear solve."""
function solve_h(es::ExchangeSetup; verbose=true, kap=r -> 0.0)
    t0 = time()
    b = rhs_T31(es)
    t1 = time()
    A = assemble(es.a, es.wg, es.mg, es.ab, es.rg; kap=kap)
    t2 = time()
    F = lu(A)
    X = reshape(F \ vec(b), NC(es.ab), length(es.rg.r))
    t3 = time()
    verbose && println("  rhs $(round(t1-t0,digits=1))s  assemble $(round(t2-t1,digits=1))s  solve $(round(t3-t2,digits=1))s  N=$(length(es.rg.r)) NC=$(NC(es.ab)) nnz=$(nnz(A))")
    return X, A, b
end

"""
Source for mode 4 on the real radial points 1..ireal and the angular nodes:
returns S4[i, j] (ingoing coords) and the gauge-condition residual |Z| relative to |∇h|.
"""
function mode4_source(es::ExchangeSetup, X)
    ab = es.ab; rg = es.rg
    w1 = real(es.md1.w); m1 = es.md1.m
    S1, dS1, d2S1 = es.S1
    Hj = jet_arrays(es.wg, es.mg)
    nth = length(ab.T.th)
    S4 = zeros(ComplexF64, es.ireal, nth)
    Zrel = zeros(es.ireal)
    for i in 1:es.ireal
        r = real(rg.r[i])
        f1 = radial_f(es.md1, r)
        zmax = 0.0; dmax = 0.0
        for (j, th) in enumerate(ab.T.th)
            bg = bg_point(es.a, r, th)
            D1, D2 = cov_derivs(bg, Hj...)
            _, _, Z = einstein_jets(bg, D1, D2)
            J = h_jets(ab, rg, X, i, j)
            d1, dd = scalar_jets(f1..., S1[j], dS1[j], d2S1[j], w1, m1)
            K = box_coupling(bg, Hj[1], D1, Z, d1, dd)
            S4[i, j] = sum(K .* J)
            zmax = max(zmax, maximum(abs(sum(Z[ν, :] .* J)) for ν in 1:4))
            dmax = max(dmax, maximum(abs(sum(D1[ρ, μ, ν, :] .* J)) for ρ in 1:4, μ in 1:4, ν in 1:4))
        end
        Zrel[i] = zmax / max(dmax, 1e-300)
    end
    return S4, Zrel
end

# ---------------------------------------------------------------- Newtonian and contact sources (BL, same modes)
R_BL(md::KerrMode, r) = radial_f(md, r)[1] * exp(-im * real(md.w) * rstar(md.a, r) + im * md.m * rsharp(md.a, r))

cumtrap(f, r) = (c = zeros(ComplexF64, length(r)); for i in 2:length(r); c[i] = c[i-1] + 0.5 * (f[i] + f[i-1]) * (r[i] - r[i-1]); end; c)
revtrap(f, r) = (c = zeros(ComplexF64, length(r)); for i in length(r)-1:-1:1; c[i] = c[i+1] + 0.5 * (f[i] + f[i+1]) * (r[i+1] - r[i]); end; c)

"""
Radial sources (projected on conj(Z4) with the Sigma weight) on the fine BL grid r:
  SI : R1^2 R3* <S1 S1 S3* Z4* Σ>                          (contact; same normalisation as sg_vs_si.jl up to preFac/unitMatch)
  N  : R1 Φ[R1 R3* S1 S3*]  projected                       (Newtonian, fixed BH)
  NR : N with the BH recoil (BH-frame indirect term): Q_1 -> Q_1 - r ∫_0^∞ f_1 dr'
"""
function newtonian_sources(md1, md3, l4, m4, w4, r; nth=80, Lextra=6)
    a = md1.a; alph = md1.alph
    M_ = md1.m - md3.m
    Lmax = md1.l + md3.l + Lextra
    T1 = SWSHTable(md1.m, maximum(keys(md1.sph)), nth; smax=0)
    T3 = SWSHTable(md3.m, maximum(keys(md3.sph)), nth; smax=0)
    _, sph4 = spheroidal_coeffs(l4, m4, a^2 * (w4^2 - alph^2))
    T4 = SWSHTable(m4, maximum(keys(sph4)), nth; smax=0)
    TL = SWSHTable(M_, Lmax, nth; smax=0)
    S1 = spheroidal_on_nodes(md1.sph, T1)[1]; S3 = spheroidal_on_nodes(md3.sph, T3)[1]; S4 = spheroidal_on_nodes(sph4, T4)[1]
    x = T1.x; w = T1.w
    R1 = [R_BL(md1, ri) for ri in r]; R3 = [R_BL(md3, ri) for ri in r]
    ang(v) = 2pi * sum(w .* v)
    # contact
    A0 = ang(S1 .* S1 .* conj.(S3) .* conj.(S4)); A2 = ang(S1 .* S1 .* conj.(S3) .* conj.(S4) .* x .^ 2)
    T_SI = R1 .* R1 .* conj.(R3) .* (A0 .* r .^ 2 .+ A2 * a^2)
    # Newtonian: Φ = sum_L 4π/(2L+1) c_L Q_L(r) Y_LM,  c_L = <S1 S3* Y_LM*>
    rho = R1 .* conj.(R3)
    T_N = zeros(ComplexF64, length(r)); T_NR = zeros(ComplexF64, length(r))
    for L in abs(M_):Lmax
        cL = ang(S1 .* conj.(S3) .* TL.Y[(0, L)])
        abs(cL) < 1e-14 && continue
        ci = cumtrap(rho .* r .^ (L + 2), r); co = revtrap(rho .* r .^ (1 - L), r)
        Q = r .^ (-(L + 1)) .* ci .+ r .^ L .* co
        QR = L == 1 ? Q .- r .* (co[1] + 0 * ci[1]) : Q          # recoil: subtract r ∫_0^∞ f_1 dr' (co[1] = ∫ over all r)
        BL0 = ang(S1 .* TL.Y[(0, L)] .* conj.(S4)); BL2 = ang(S1 .* TL.Y[(0, L)] .* conj.(S4) .* x .^ 2)
        T_N .+= (4pi / (2L + 1)) * cL .* Q .* (BL0 .* r .^ 2 .+ BL2 * a^2)
        T_NR .+= (4pi / (2L + 1)) * cL .* QR .* (BL0 .* r .^ 2 .+ BL2 * a^2)
    end
    T_N .*= R1; T_NR .*= R1
    return T_SI, T_N, T_NR, sph4
end

"""GR source projected the same way, on the real points of the graviton grid (BL phase convention)."""
function gr_radial_source(es::ExchangeSetup, S4, l4, m4, w4, sph4)
    a = es.a; alph = es.alph
    ab = es.ab; nth = length(ab.T.th)
    T4 = SWSHTable(m4, maximum(keys(sph4)), nth; smax=0)
    Z4 = spheroidal_on_nodes(sph4, T4)[1]
    r = real.(es.rg.r[1:es.ireal])
    out = zeros(ComplexF64, length(r))
    for i in eachindex(r)
        ph = exp(-im * w4 * rstar(a, r[i]) + im * m4 * rsharp(a, r[i]))       # ingoing -> BL
        i == 1 && (ph = 0.0)                                                  # r = r+ exactly: phase singular (unused)
        out[i] = 2pi * sum(ab.T.w .* S4[i, :] .* conj.(Z4) .* (r[i]^2 .+ a^2 .* ab.T.x .^ 2)) * ph / (-4 * alph^4)
    end
    return r, out
end

function interp_rho(es::ExchangeSetup, rsrc, vals, r)
    # cubic interpolation in rho = log((r - r+)/delta + 1) (uniform grid); zero outside
    rho = es.rg.rho[1:length(rsrc)]
    itR = CubicSplineInterpolation(range(rho[1], rho[end], length=length(rho)), real.(vals))
    itI = CubicSplineInterpolation(range(rho[1], rho[end], length=length(rho)), imag.(vals))
    rp = rplus(es.a)
    out = zeros(ComplexF64, length(r))
    for (k, x) in enumerate(r)
        y = log((x - rp) / es.delta + 1)
        (y < rho[2] || y > rho[end]) && continue
        out[k] = itR(y) + im * itI(y)
    end
    return out
end

