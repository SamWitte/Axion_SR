"""
    teukolsky.jl

Homogeneous solutions of the radial Teukolsky equation on Kerr and the source
integrals that give the gravitational-wave amplitude at infinity. This is the
radial machinery behind the relativistic GW rates (gw_rates_relativistic).

Conventions: Sasaki & Tagoshi, Living Rev. Relativ. 6 (2003) 6 [ST03], Sec. 2,
with G = c = M = 1, a = J/M^2, Boyer-Lindquist coordinates and the Kinnersley
tetrad. For spin weight s the radial equation (ST03 eq. 14 for s = -2) is

    Δ^{-s} d/dr(Δ^{s+1} dR/dr) + U(r) R = 𝒯,
    U = [K^2 - 2is(r-1)K]/Δ + 4isωr - λ,   K = (r^2+a^2)ω - am,

λ the eigenvalue of the spin-weighted spheroidal harmonic in Teukolsky's
convention (λ -> (l-s)(l+s+1) as aω -> 0; SpinWeightedSpheroidalHarmonics.jl
returns exactly this). s = -2 is ψ_4 = ρ^4 Σ_{lm} e^{-iωt+imφ} S_lm(θ) R_lm(r),
ρ = (r - ia cosθ)^{-1}, with S normalised to ∫|S|^2 dΩ = 1.

Homogeneous solutions (ST03 eqs. 19-20, normalised here with B^trans = C^trans = 1):
    R^in  -> Δ^{-s} e^{-ikr*}                         (r -> r_+),   k = ω - mΩ_H
    R^in  -> B^ref r^{-2s-1} e^{iωr*} + B^inc r^{-1} e^{-iωr*}   (r -> ∞)
    R^up  -> r^{-2s-1} e^{iωr*}                       (r -> ∞)
W = Δ^{s+1}(R^in R^up' - R^up R^in') is constant and equals 2iω B^inc.
The amplitude of the outgoing wave at infinity for a source 𝒯 is (ST03 eq. 25)
    Z^∞ = (1/W) ∫ R^in 𝒯 Δ^{-2} dr,   ψ_4 -> (1/r) Σ Z^∞ S e^{iω(r*-t)+imφ},
and the energy flux is dE/dt = Σ |Z^∞|^2/(4πω^2) (ST03 eq. 48).

Numerics: both boundary series are generated from the ODE itself with a small
truncated-Laurent-series algebra (no hand-derived recurrences), and the ODE is
integrated with OrdinaryDiffEq in any floating type (Float64 or BigFloat).
R^in is integrated outward from the horizon, which is stable; R^up inward from
r_inf, which picks up the growing solution as (r_inf/r)^4, so it is only
carried to r_m = r_inf/10 where the Wronskian is taken.

Sources: `teuk_source_kernels` gives the kernels multiplying the tetrad
components T_nn, T_n m̄, T_m̄ m̄ after all derivatives have been moved onto
R^in and S (ST03 eqs. 36-44 with the particle's C_ab replaced by 1), so one
routine serves a point particle (test: circular equatorial orbits against the
post-Newtonian luminosity, ST03 eqs. 178 and 181) and the extended scalar cloud.
"""

using OrdinaryDiffEq
using SpinWeightedSpheroidalHarmonics
using LinearAlgebra

# ----------------------------------------------------------------------------
# Kerr geometry (M = 1)
# ----------------------------------------------------------------------------

kerr_rpm(a) = (1 + sqrt(1 - a^2), 1 - sqrt(1 - a^2))
kerr_Delta(r, a) = r^2 - 2r + a^2
kerr_OmegaH(a) = a / (2 * kerr_rpm(a)[1])
kerr_K(r, a, ω, m) = (r^2 + a^2) * ω - a * m

"""Tortoise coordinate, ST03 eq. 21 (integration constant as there)."""
function kerr_rstar(r, a)
    rp, rm = kerr_rpm(a)
    if rm == 0
        return r + 2rp / (rp - rm) * log((r - rp) / 2)
    end
    return r + 2rp / (rp - rm) * log((r - rp) / 2) - 2rm / (rp - rm) * log((r - rm) / 2)
end

# ----------------------------------------------------------------------------
# Truncated Laurent series in one variable z: Σ_{n=lo}^{hi} c_n z^n
# ----------------------------------------------------------------------------

struct LSeries{T}
    lo::Int
    c::Vector{T}   # coefficients of z^lo, z^(lo+1), ...
    hi::Int        # exponents above hi are unknown (dropped)
end

ls_coef(f::LSeries, n::Int) = (f.lo <= n <= f.lo + length(f.c) - 1) ? f.c[n - f.lo + 1] : zero(eltype(f.c))

function lseries(lo::Int, c::AbstractVector{T}, hi::Int) where T
    n = min(length(c), hi - lo + 1)
    n <= 0 && return LSeries{T}(hi, T[zero(T)], hi)
    return LSeries{T}(lo, collect(c[1:n]), hi)
end
ls_const(x::T, hi::Int) where T = lseries(0, T[x], hi)

function Base.:+(f::LSeries{T}, g::LSeries{T}) where T
    lo = min(f.lo, g.lo); hi = min(f.hi, g.hi)
    return lseries(lo, T[ls_coef(f, n) + ls_coef(g, n) for n in lo:hi], hi)
end
Base.:-(f::LSeries{T}) where T = LSeries{T}(f.lo, -f.c, f.hi)
Base.:-(f::LSeries, g::LSeries) = f + (-g)
Base.:*(x::Number, f::LSeries{T}) where T = LSeries{T}(f.lo, T(x) .* f.c, f.hi)
Base.:*(f::LSeries, x::Number) = x * f
Base.:+(f::LSeries{T}, x::Number) where T = f + ls_const(T(x), f.hi)
Base.:+(x::Number, f::LSeries) = f + x
Base.:-(f::LSeries, x::Number) = f + (-x)
Base.:-(x::Number, f::LSeries) = (-f) + x

function Base.:*(f::LSeries{T}, g::LSeries{T}) where T
    lo = f.lo + g.lo
    hi = min(f.hi + g.lo, g.hi + f.lo)       # beyond this the truncated tails contribute
    c = zeros(T, max(hi - lo + 1, 1))
    for (i, fi) in enumerate(f.c), (j, gj) in enumerate(g.c)
        n = f.lo + i - 1 + g.lo + j - 1
        n > hi && continue
        c[n - lo + 1] += fi * gj
    end
    return lseries(lo, c, hi)
end

"""d/dz of the series."""
function ls_deriv(f::LSeries{T}) where T
    c = T[(f.lo + i - 1) * f.c[i] for i in 1:length(f.c)]
    return lseries(f.lo - 1, c, f.hi - 1)
end

"""1/f; the lowest coefficient must be nonzero (f = z^lo (c_1 + c_2 z + ...))."""
function ls_inv(f::LSeries{T}) where T
    f.c[1] == 0 && error("ls_inv: lowest coefficient vanishes; construct the series with the exact lowest exponent")
    N = f.hi - f.lo
    d = zeros(T, N + 1)
    d[1] = 1 / f.c[1]
    for n in 1:N
        acc = zero(T)
        for k in 1:n
            ck = k + 1 <= length(f.c) ? f.c[k + 1] : zero(T)
            acc += ck * d[n - k + 1]
        end
        d[n + 1] = -acc / f.c[1]
    end
    return lseries(-f.lo, d, N - f.lo)
end

function ls_eval(f::LSeries, z)
    acc = zero(promote_type(eltype(f.c), typeof(z)))
    zp = z^f.lo
    for ci in f.c
        acc += ci * zp
        zp *= z
    end
    return acc
end

"""
    ls_frobenius(A2, A1, A0, N, D; tol) -> f_0..f_N

Power-series solution f = Σ f_n z^n (f_0 = 1) of A2 f'' + A1 f' + A0 f = 0 with
Laurent-series coefficients A_k. The term f_n first appears at the power z^(D+n),
with D fixed by the ansatz (the exponent of z already extracted into g, so the
indicial equation has the root 0). Coefficients that would act below z^D must
vanish; they are checked against the first indicial coefficient.
"""
function ls_frobenius(A2::LSeries{T}, A1::LSeries{T}, A0::LSeries{T}, N::Int, D::Int; tol=1e-9) where T
    indicial(n) = n * (n - 1) * ls_coef(A2, D + 2) + n * ls_coef(A1, D + 1) + ls_coef(A0, D)
    ref = abs(indicial(1))
    ref > 0 || error("ls_frobenius: vanishing indicial coefficient")
    low = max(abs(indicial(0)),
              maximum(max(abs(ls_coef(A2, M + 2)), abs(ls_coef(A1, M + 1)), abs(ls_coef(A0, M)))
                      for M in (min(A2.lo - 2, A1.lo - 1, A0.lo, D - 1)):(D - 1)))
    low <= tol * ref || error("ls_frobenius: the ansatz is not a solution at the lowest orders (|low|/|I(1)| = $(low / ref))")
    min(A2.hi - 2, A1.hi - 1, A0.hi) >= D + N ||
        error("ls_frobenius: coefficient series too short for $(N) terms")
    f = zeros(T, N + 1)
    f[1] = one(T)
    for n in 1:N
        M = D + n
        acc = zero(T)
        for j in 0:(n - 1)
            acc += f[j + 1] * (j * (j - 1) * ls_coef(A2, M - j + 2) + j * ls_coef(A1, M - j + 1) + ls_coef(A0, M - j))
        end
        f[n + 1] = -acc / indicial(n)
    end
    return f
end

# ----------------------------------------------------------------------------
# Boundary series of the radial Teukolsky equation
# ----------------------------------------------------------------------------

# Operator on f for R = g f: A2 f'' + A1 f' + A0 f with p = g'/g,
#   A2 = Δ,  A1 = 2Δp + (s+1)Δ',  A0 = Δ(p' + p^2) + (s+1)Δ'p + U.
function _teuk_operator_series(s, m, ω, a, λ, Δ, dΔ, r, r2a2, p, dp, ::Type{T}) where T
    invΔ = ls_inv(Δ)
    K = ω * r2a2 - a * m
    U = (K * K - (2im * s) * (r - 1) * K) * invΔ + (4im * s * ω) * r - λ
    A2 = Δ
    A1 = 2 * (Δ * p) + (s + 1) * dΔ
    A0 = Δ * (dp + p * p) + (s + 1) * (dΔ * p) + U
    return A2, A1, A0
end

"""
    teuk_series_horizon(s, m, ω, a, λ, N, T) -> coefficients f_n

R^in = Δ^{-s} e^{-ikr*} Σ_n f_n (r - r_+)^n, f_0 = 1 (converges for r - r_+ < r_+ - r_-).
"""
function teuk_series_horizon(s::Int, m::Int, ω, a, λ, N::Int, ::Type{T}) where T
    C = Complex{T}
    rp, rm = kerr_rpm(T(a))
    k = T(ω) - m * T(a) / (2rp)
    hi = N + 8
    Δ = lseries(1, C[rp - rm, 1], hi)                    # z (z + r_+ - r_-)
    dΔ = lseries(0, C[2(rp - 1), 2], hi)                 # 2(r - 1)
    r = lseries(0, C[rp, 1], hi)
    r2a2 = lseries(0, C[rp^2 + a^2, 2rp, 1], hi)
    invΔ = ls_inv(Δ)
    p = (-s * dΔ - (im * k) * r2a2) * invΔ               # g'/g, g = Δ^{-s} e^{-ikr*}
    dp = ls_deriv(p)
    A2, A1, A0 = _teuk_operator_series(s, m, C(ω), C(a), C(λ), Δ, dΔ, r, r2a2, p, dp, C)
    return ls_frobenius(A2, A1, A0, N, -1; tol=sqrt(eps(T)) * 1e3)   # A2 ~ z, A1 ~ 1, A0 ~ 1/z
end

"""
    teuk_series_infinity(s, m, ω, a, λ, N, T) -> coefficients f_n

R^up = r^{-2s-1} e^{iωr*} Σ_n f_n r^{-n}, f_0 = 1 (asymptotic series; use at ω r ≫ 1).
"""
function teuk_series_infinity(s::Int, m::Int, ω, a, λ, N::Int, ::Type{T}) where T
    C = Complex{T}
    hi = N + 8
    u = lseries(1, C[1], hi)                               # u = 1/r
    r = lseries(-1, C[1], hi)
    D = lseries(0, C[1, -2, a^2], hi)                      # Δ = u^{-2} D
    Δ = lseries(-2, C[1, -2, a^2], hi)
    dΔ = lseries(-1, C[2, -2], hi)                         # 2(r - 1)
    r2a2 = lseries(-2, C[1, 0, a^2], hi)
    invD = ls_inv(D)
    p = (-2s - 1) * u + (im * ω) * (lseries(0, C[1, 0, a^2], hi) * invD)   # g'/g, g = r^{-2s-1} e^{iωr*}
    dp = -(u * u * ls_deriv(p))                                             # d/dr = -u^2 d/du
    A2, A1, A0 = _teuk_operator_series(s, m, C(ω), C(a), C(λ), Δ, dΔ, r, r2a2, p, dp, C)
    # f' = -u^2 f_u, f'' = u^4 f_uu + 2u^3 f_u
    B2 = (u * u * u * u) * A2
    B1 = 2 * ((u * u * u) * A2) - (u * u) * A1
    return ls_frobenius(B2, B1, A0, N, -1; tol=sqrt(eps(T)) * 1e3)   # B2 ~ u^2, B1 ~ 1, A0 ~ 1
end

"""R^in and dR^in/dr from the horizon series at r."""
function teuk_eval_horizon(f::Vector, s::Int, m::Int, ω, a, r)
    rp, rm = kerr_rpm(a)
    k = ω - m * a / (2rp)
    x = r - rp
    F = zero(eltype(f)); dF = zero(eltype(f))
    xp = one(x)
    for (n, fn) in enumerate(f)
        F += fn * xp
        n > 1 && (dF += (n - 1) * fn * xp / x)
        xp *= x
    end
    Δ = kerr_Delta(r, a)
    g = Δ^(-s) * exp(-im * k * kerr_rstar(r, a))
    p = (-s * 2(r - 1) - im * k * (r^2 + a^2)) / Δ
    return g * F, g * (dF + p * F)
end

"""R^up and dR^up/dr from the asymptotic series at r."""
function teuk_eval_infinity(f::Vector, s::Int, m::Int, ω, a, r)
    u = 1 / r
    F = zero(eltype(f)); dFu = zero(eltype(f))
    up = one(u)
    for (n, fn) in enumerate(f)
        F += fn * up
        n > 1 && (dFu += (n - 1) * fn * up / u)
        up *= u
    end
    Δ = kerr_Delta(r, a)
    g = r^(-2s - 1) * exp(im * ω * kerr_rstar(r, a))
    p = (-2s - 1) / r + im * ω * (r^2 + a^2) / Δ
    return g * F, g * (-u^2 * dFu + p * F)
end

# ----------------------------------------------------------------------------
# ODE integration and the TeukRadial container
# ----------------------------------------------------------------------------

function teuk_rhs!(du, u, p, r)
    s, m, ω, a, λ = p
    Δ = r^2 - 2r + a^2
    K = (r^2 + a^2) * ω - a * m
    U = (K^2 - 2im * s * (r - 1) * K) / Δ + 4im * s * ω * r - λ
    du[1] = u[2]
    du[2] = -((s + 1) * 2(r - 1) * u[2] + U * u[1]) / Δ
    return nothing
end

"""Second derivative from the ODE, given (R, R')."""
function teuk_R2(s, m, ω, a, λ, r, R, dR)
    Δ = r^2 - 2r + a^2
    K = (r^2 + a^2) * ω - a * m
    U = (K^2 - 2im * s * (r - 1) * K) / Δ + 4im * s * ω * r - λ
    return -((s + 1) * 2(r - 1) * dR + U * R) / Δ
end

"""
Spheroidal eigenvalue (Teukolsky's λ) and the unit-normalised harmonic
S(θ; theta_derivative) for spin weight s and spheroidicity c = aω.
"""
function teuk_angular(s::Int, l::Int, m::Int, c)
    cF = Float64(c)
    lmin = max(abs(m), abs(s))
    N = max(10, l - lmin + 8)
    λ = spin_weighted_spheroidal_eigenvalue(s, l, m, cF; N=N)
    Sf = spin_weighted_spheroidal_harmonic(s, l, m, cF; N=N)
    S(θ; theta_derivative::Int=0) = real(Sf(Float64(θ), 0.0; theta_derivative=theta_derivative))
    return λ, S
end

struct TeukRadial{T, S1, S2}
    s::Int; l::Int; m::Int
    ω::T; a::T; λ::T; k::T
    fh::Vector{Complex{T}}; r_h0::T          # horizon series, used for r <= r_h0
    finf::Vector{Complex{T}}; r_inf::T       # asymptotic series, used for r >= r_inf
    r_m::T                                   # R^up is available on [r_m, r_inf]
    sol_in::S1                               # ODE solution of R^in on [r_h0, r_inf]
    sol_up::S2                               # ODE solution of R^up on [r_m, r_inf]
    W::Complex{T}                            # Δ^{s+1}(R^in R^up' - R^up R^in') = 2iω B^inc
    W_spread::T                              # relative variation of W over [r_m, r_inf] (diagnostic)
    Sfun                                     # S(θ; theta_derivative)
end

"""
    teuk_radial(s, l, m, ω, a; T=Float64, N=40, rtol, atol, ωr_inf=25, r_inf_min=200,
                hfrac=0.1, ωr_m=0.5, λ=nothing) -> TeukRadial

Homogeneous solutions R^in (horizon to r_inf) and R^up (r_inf down to
r_m = min(ωr_m/|ω|, r_inf/10)) for real ω ≠ 0. The horizon series is used out
to r_+ + hfrac (r_+ - r_-), the asymptotic series from r_inf = max(ωr_inf/|ω|,
r_inf_min). `λ` overrides the spheroidal eigenvalue (e.g. a BigFloat value).
"""
function teuk_radial(s::Int, l::Int, m::Int, ω, a; T=Float64, N::Int=40, rtol=nothing, atol=nothing,
                     ωr_inf=25.0, r_inf_min=200.0, hfrac=0.1, ωr_m=0.5, λ=nothing, nW=8)
    ω == 0 && error("teuk_radial: ω = 0 (static) is not handled")
    ωT, aT = T(ω), T(a)
    rp, rm = kerr_rpm(aT)
    λT = λ === nothing ? T(teuk_angular(s, l, m, aT * ωT)[1]) : T(λ)
    _, Sfun = teuk_angular(s, l, m, aT * ωT)
    k = ωT - m * aT / (2rp)
    # Vern9 steps scale as rtol^(-1/9); for BigFloat the default stops well short of eps(T)
    rtol === nothing && (rtol = T == Float64 ? 1e-13 : T(10)^(-min(25, precision(T) ÷ 4)))
    atol === nothing && (atol = rtol * T(1e-12))
    C = Complex{T}
    prm = (s, m, ωT, aT, λT)

    # R^in: series at r_h0, then outward to r_inf
    fh = teuk_series_horizon(s, m, ωT, aT, λT, N, T)
    r_h0 = rp + T(hfrac) * (rp - rm)
    r_inf = max(T(ωr_inf) / abs(ωT), T(r_inf_min))
    R0, dR0 = teuk_eval_horizon(fh, s, m, ωT, aT, r_h0)
    prob_in = ODEProblem(teuk_rhs!, C[R0, dR0], (r_h0, r_inf), prm)
    sol_in = solve(prob_in, Vern9(); reltol=rtol, abstol=atol, maxiters=10^8)
    OrdinaryDiffEq.SciMLBase.successful_retcode(sol_in) || error("teuk_radial: R^in integration failed ($(sol_in.retcode))")

    # R^up: asymptotic series at r_inf, then inward to r_m
    finf = teuk_series_infinity(s, m, ωT, aT, λT, N, T)
    # R^up is carried into the near zone, where the Wronskian has no cancellation between the
    # r^{-2s-1} parts of the two solutions (the contamination of R^up by the solution that
    # grows inward is ~ (ωr_inf)^4 × the starting error)
    r_m = max(min(T(ωr_m) / abs(ωT), r_inf / 10), 2 * r_h0)
    R1, dR1 = teuk_eval_infinity(finf, s, m, ωT, aT, r_inf)
    prob_up = ODEProblem(teuk_rhs!, C[R1, dR1], (r_inf, r_m), prm)
    sol_up = solve(prob_up, Vern9(); reltol=rtol, abstol=atol, maxiters=10^8)
    OrdinaryDiffEq.SciMLBase.successful_retcode(sol_up) || error("teuk_radial: R^up integration failed ($(sol_up.retcode))")

    # Wronskian on the overlap
    rs = exp.(range(log(r_m), log(r_inf), length=nW))
    Ws = map(rs) do r
        ui = sol_in(r); uu = sol_up(r)
        kerr_Delta(r, aT)^(s + 1) * (ui[1] * uu[2] - uu[1] * ui[2])
    end
    W = Ws[1]
    W_spread = maximum(abs.(Ws .- W)) / abs(W)
    return TeukRadial(s, l, m, ωT, aT, λT, k, fh, r_h0, finf, r_inf, r_m, sol_in, sol_up, W, W_spread, Sfun)
end

"""(R^in, R^in', R^in'') at r."""
function teuk_in(tr::TeukRadial, r)
    if r <= tr.r_h0
        R, dR = teuk_eval_horizon(tr.fh, tr.s, tr.m, tr.ω, tr.a, r)
    else
        r <= tr.r_inf || error("teuk_in: r = $r beyond r_inf = $(tr.r_inf)")
        u = tr.sol_in(r)
        R, dR = u[1], u[2]
    end
    return R, dR, teuk_R2(tr.s, tr.m, tr.ω, tr.a, tr.λ, r, R, dR)
end

"""(R^up, R^up', R^up'') at r (only for r >= r_m)."""
function teuk_up(tr::TeukRadial, r)
    if r >= tr.r_inf
        R, dR = teuk_eval_infinity(tr.finf, tr.s, tr.m, tr.ω, tr.a, r)
    else
        r >= tr.r_m || error("teuk_up: r = $r below r_m = $(tr.r_m)")
        u = tr.sol_up(r)
        R, dR = u[1], u[2]
    end
    return R, dR, teuk_R2(tr.s, tr.m, tr.ω, tr.a, tr.λ, r, R, dR)
end

"""B^inc of R^in (with R^in -> Δ^{-s} e^{-ikr*} at the horizon)."""
teuk_Binc(tr::TeukRadial) = tr.W / (2im * tr.ω)

# ----------------------------------------------------------------------------
# Source kernels (s = -2)
# ----------------------------------------------------------------------------

"""
    teuk_source_kernels(tr, r, θ) -> (K_nn, K_nm̄, K_m̄m̄)

Kernels of the amplitude at infinity for a source with tetrad components
T_ab(t, r, θ, φ) = T̂_ab(r, θ) e^{-iωt + imφ} (Kinnersley tetrad, ST03 eqs. 36-44
with the particle's C_ab -> 1 and S unit-normalised over the sphere):

    Z^∞ = (2π)^2/W ∫ dr ∫ sinθ dθ [K_nn T̂_nn + K_nm̄ T̂_nm̄ + K_m̄m̄ T̂_m̄m̄],

all r- and θ-derivatives having been moved onto R^in and S. For a point
particle with T_ab = μ C_ab δ(r-r_0)δ(θ-θ_0)δ(φ-φ(t))/sinθ on a circular orbit
(φ = Ωt, ω = mΩ) this collapses to Z^∞ = 2πμ/W [K_nn C_nn + ...] at (r_0, θ_0).
"""
function teuk_source_kernels(tr::TeukRadial, r, θ)
    tr.s == -2 || error("teuk_source_kernels: only s = -2 (ψ_4) sources are implemented")
    a, ω, m = tr.a, tr.ω, tr.m
    S = tr.Sfun(θ); dS = tr.Sfun(θ; theta_derivative=1); d2S = tr.Sfun(θ; theta_derivative=2)
    Rin, dRin, d2Rin = teuk_in(tr, r)
    sθ, cθ = sin(θ), cos(θ)
    ρ = 1 / (r - im * a * cθ); ρb = 1 / (r + im * a * cθ)
    Δ = kerr_Delta(r, a)
    K = kerr_K(r, a, ω, m)
    KΔ = K / Δ
    dKΔ = 2r * ω / Δ - K * 2(r - 1) / Δ^2
    # L_s^† = ∂_θ - m/sinθ + aω sinθ + s cotθ
    q2 = -m / sθ + a * ω * sθ + 2cθ / sθ
    dq2 = m * cθ / sθ^2 + a * ω * cθ - 2 / sθ^2
    L2S = dS + q2 * S
    dL2S = d2S + q2 * dS + dq2 * S
    # L_1^†[ρ^{-4} L_2^†(ρ^3 S)],  ∂_θ ρ = -ia sinθ ρ^2,  ∂_θ ρ^{-1} = ia sinθ
    G = L2S / ρ - 3im * a * sθ * S
    dG = im * a * sθ * L2S + dL2S / ρ - 3im * a * cθ * S - 3im * a * sθ * dS
    L1G = dG + (-m / sθ + a * ω * sθ + cθ / sθ) * G
    Ann0 = -(2 / Δ^2) * ρ^-2 * ρb^-1 * L1G
    Amn0 = (2sqrt(2) / Δ) * ρ^-3 * (L2S * (im * KΔ + ρ + ρb) - a * sθ * S * KΔ * (ρb - ρ))
    Amm0 = -ρ^-3 * ρb * S * (-im * dKΔ - KΔ^2 + 2im * ρ * KΔ)
    Amn1 = (2sqrt(2) / Δ) * ρ^-3 * (L2S + im * a * sθ * (ρb - ρ) * S)
    Amm1 = -2 * ρ^-3 * ρb * S * (im * KΔ + ρ)
    Amm2 = -ρ^-3 * ρb * S
    Knn = Rin * Ann0
    Kmn = Rin * Amn0 - dRin * Amn1
    Kmm = Rin * Amm0 - dRin * Amm1 + d2Rin * Amm2
    return Knn, Kmn, Kmm
end

# ----------------------------------------------------------------------------
# Point particle on a circular equatorial orbit (validation)
# ----------------------------------------------------------------------------

"""
Prograde circular equatorial geodesic at Boyer-Lindquist radius r0 (M = 1):
(E, Lz, Ω, dt/dτ) per unit particle mass.
"""
function teuk_circular_orbit(a, r0)
    v = 1 / sqrt(r0)
    den = sqrt(1 - 3v^2 + 2a * v^3)
    E = (1 - 2v^2 + a * v^3) / den
    Lz = r0 * v * (1 - 2a * v^3 + a^2 * v^4) / den
    Ω = 1 / (r0^(3 // 2) + a)
    Δ = kerr_Delta(r0, a)
    tdot = (-a * (a * E - Lz) + (r0^2 + a^2) * (E * (r0^2 + a^2) - a * Lz) / Δ) / r0^2
    return E, Lz, Ω, tdot
end

"""
Z^∞ for the mode (l, m) of `tr` from a unit-mass particle on the prograde
circular equatorial orbit at r0; requires tr.ω = m Ω(r0).
"""
function teuk_circular_Z(tr::TeukRadial, r0)
    E, Lz, Ω, tdot = teuk_circular_orbit(tr.a, r0)
    abs(tr.ω - tr.m * Ω) <= 1e-12 * abs(tr.ω) || error("teuk_circular_Z: tr.ω must equal m Ω(r0)")
    a = tr.a
    X = E * (r0^2 + a^2) - a * Lz          # -2Σ n·u
    Y = a * E - Lz                           # (m̄·u) = i ρ Y /√2 at the equator
    Cnn = X^2 / (4 * r0^6 * tdot)
    Cmn = -im * X * Y / (2sqrt(2) * r0^5 * tdot)
    Cmm = -Y^2 / (2 * r0^4 * tdot)
    Knn, Kmn, Kmm = teuk_source_kernels(tr, r0, oftype(r0, π) / 2)
    return 2π * (Knn * Cnn + Kmn * Cmn + Kmm * Cmm) / tr.W
end

"""
    teuk_circular_flux(a, r0; lmax=6, T=Float64, kwargs...) -> (P, modes)

GW energy flux to infinity, in units of (μ/M)^2, from a unit-mass particle on
the prograde circular equatorial orbit at r0, summed over 2 <= l <= lmax and
m = -l..l (m ≠ 0). `modes` maps (l, m) to its flux. kwargs go to teuk_radial.
"""
function teuk_circular_flux(a, r0; lmax::Int=6, T=Float64, kwargs...)
    _, _, Ω, _ = teuk_circular_orbit(T(a), T(r0))
    modes = Dict{Tuple{Int,Int}, T}()
    P = zero(T)
    for l in 2:lmax, m in -l:l
        m == 0 && continue
        ω = m * Ω
        tr = teuk_radial(-2, l, m, ω, T(a); T=T, kwargs...)
        Z = teuk_circular_Z(tr, T(r0))
        modes[(l, m)] = abs2(Z) / (4π * ω^2)
        P += modes[(l, m)]
    end
    return P, modes
end

"""
    pn_flux_circular(a, r0) -> dE/dt / (32/5 x^10)

Post-Newtonian luminosity of a particle on a circular equatorial orbit, ST03
eq. 178 (Schwarzschild, to x^11) plus the spin-dependent terms of eq. 181
(to x^8), x = Ω^{1/3}. For a ≠ 0 the spin-independent terms beyond x^8 are not
spin-complete, so the a ≠ 0 result is accurate to O(x^9).
"""
function pn_flux_circular(a, r0)
    Ω = 1 / (r0^(3 // 2) + a)
    x = Ω^(1 // 3)
    q = a
    γ = Base.MathConstants.eulergamma
    f = 1 - 1247 // 336 * x^2 + 4π * x^3 - 44711 // 9072 * x^4 - 8191 // 672 * π * x^5 +
        (6643739519 // 69854400 - 1712 // 105 * γ + 16 // 3 * π^2 - 3424 // 105 * log(2) - 1712 // 105 * log(x)) * x^6 -
        16285 // 504 * π * x^7 +
        (-323105549467 // 3178375200 + 232597 // 4410 * γ - 1369 // 126 * π^2 + 39931 // 294 * log(2) -
         47385 // 1568 * log(3) + 232597 // 4410 * log(x)) * x^8
    if a == 0
        f += (265978667519 // 745113600 * π - 6848 // 105 * π * γ - 13696 // 105 * π * log(2) - 6848 // 105 * π * log(x)) * x^9 +
             (-2500861660823683 // 2831932303200 + 916628467 // 7858620 * γ - 424223 // 6804 * π^2 -
              83217611 // 1122660 * log(2) + 47385 // 196 * log(3) + 916628467 // 7858620 * log(x)) * x^10 +
             (8399309750401 // 101708006400 * π + 177293 // 1176 * π * γ + 8521283 // 17640 * π * log(2) -
              142155 // 784 * π * log(3) + 177293 // 1176 * π * log(x)) * x^11
    else
        f += -11 // 4 * q * x^3 + 33 // 16 * q^2 * x^4 - 59 // 16 * q * x^5 + (-65 // 6 * π * q + 611 // 504 * q^2) * x^6 +
             (162035 // 3888 * q + 65 // 8 * π * q^2 - 71 // 24 * q^3) * x^7 +
             (-359 // 14 * π * q + 22667 // 4536 * q^2 + 17 // 16 * q^4) * x^8
    end
    return f
end

pn_flux_newtonian(a, r0) = 32 // 5 * (1 / (r0^(3 // 2) + a))^(10 // 3)
