"""
    kerr_modes.jl

Quasi-bound states of a massive scalar on Kerr for the relativistic GW rates
(G = c = M = 1, μ = α). The radial function is Leaver's series (Dolan 2007), as
in solve_sr_rates.jl / find_im_part, in the horizon-regular ingoing form

    ψ = f(r) S(θ) e^{-iωv + imφ̃},   f = R_BL e^{iωr* - imr#},
    R_BL = (r-r_+)^{-iσ} (r-r_-)^{iσ+χ-1} e^{qr} Σ_n b_n u^n,   u = (r-r_+)/(r-r_-),

with v = t + r*, φ̃ = φ + r#, dr#/dr = a/Δ. The minimal solution b_n comes from
the continued-fraction recurrence run backward and the complex eigenvalue ω from a
secant iteration on its eigenvalue condition, with the exact spheroidal
eigenvalue Λ of (L² - c² cos²θ) S = Λ S, c² = a²(ω² - μ²), from diagonalising
in the Y_lm basis (S = Σ b_l' Y_l'm, so θ-derivatives are exact).

Precision of the eigenvalue (2026-11-21). The first ~n terms of the recurrence of level
(n, l, m) form a binomial hump (|b_k/b_0| up to C(n-1, (n-1)/2)), and the slope of the
eigenvalue condition with respect to ω falls as Π_{j<n} j/sqrt(j² + 4σ²),
σ = (ω - mΩ_H)/(4π T_H). In Float64 the rounding noise of the condition therefore
exceeds its whole variation across a level spacing for weakly bound levels at high spin
(e.g. (20,6,6) at a = 0.998, (24,10,10) at a ≥ 0.92), and a root search returns noise that
still passes a residual test. No Float64 reformulation helps (it is the conditioning of the
recurrence). The eigenvalue condition is therefore evaluated in mixed precision
(km_cf_eval): the head n ≤ J (default 300) in BigFloat, the tail in Float64 (it
contributes below 1e-25), and every root is validated a posteriori (km_eigenfrequency):
the condition is re-evaluated at the root with twice the precision, a three times longer
head and twice as many terms, each difference divided by the slope is an error of the
root in units of the level spacing, and precision, head or tail are raised until the
errors are below `tol`. The secant runs on the pole-free form D = α_0 b_1 + β_0 b_0
(β_0/α_0 - t_1 has a pole next to the root for nodeless levels).

Normalisation: N_j quanta correspond to the field √(N_j/(2ω_j)) (ψ_j + ψ_j*)
with ψ_j of unit Klein-Gordon norm on a v = const slice (finite at the
horizon; the Boyer-Lindquist slice diverges there for real ω). In the
non-relativistic limit this is ∫|ψ|² d³x = 1.

Radial overtones (N = n - l - 1 ≥ 2, kerr_mode radial=:auto): the series Σ b_n uⁿ
behaves as (1-u)^N g(u) with large alternating b_n, so at large r (u → 1) it cancels
to ~r^{-N} of its terms: in Float64 the sum is wrong beyond a few overtones, and since
the Klein-Gordon norm is set there, the normalised profile came out too small by many
orders of magnitude in the inner region that overlaps compact levels. For N ≥ 2 the
series is used only up to u = 0.2 (with the BigFloat head coefficients; the handover data
are summed in BigFloat), and f is continued outward by integrating its ODE,

    f'' + (Δ' - 2iK)/Δ f' - [iK' + μ²r² + a²ω² - 2amω + Λ]/Δ f = 0,   K = (r²+a²)ω - am,

in the variable h = f e^{-iωr}, which is regular at the horizon like f but does not carry
the far-field phase e^{iωr*} of f (2n²/α radians out to the turning point, 1e5 for n = 45
at α = 0.03):

    h'' + [2(r-1) - 2i(2ωr - am)]/Δ h' + [-q²r² + 2ω²r - 2iω - Λ]/Δ h = 0,   q² = μ² - ω².

The integration is stable up to the outer turning point (the regular solution dominates or
oscillates there) and is stopped 6/|q| beyond it; further out f follows the leading
asymptotic form (r-r_-)^{iσ+χ-1} e^{(q+iω)r}. The Klein-Gordon norm of such a mode is
integrated on the steps of the ODE solution (a fixed number of nodes cannot resolve the
n-l-1 radial nodes of a high overtone: 400 nodes gave errors of 1-20 % at n ~ 39).
"""

using LinearAlgebra
using ForwardDiff
using QuadGK: gauss
using OrdinaryDiffEq: ODEProblem, solve, Vern9

# ----------------------------------------------------------------------------
# Spherical harmonics (θ part, Condon-Shortley; works for ForwardDiff duals)
# ----------------------------------------------------------------------------

function km_theta_lm(l::Int, m::Int, θ)
    abs(m) > l && return zero(cos(θ))
    m < 0 && return (isodd(m) ? -1 : 1) * km_theta_lm(l, -m, θ)
    x = cos(θ); s = sin(θ)
    pmm = one(x) / sqrt(4π)
    for i in 1:m
        pmm *= -sqrt((2i + 1) / (2i)) * s
    end
    l == m && return pmm
    pm1 = x * sqrt(2m + 3.0) * pmm
    l == m + 1 && return pm1
    p_prev, p_cur = pmm, pm1
    for ll in (m + 2):l
        a = sqrt((4ll^2 - 1) / (ll^2 - m^2))
        b = sqrt(((ll - 1)^2 - m^2) / (4 * (ll - 1)^2 - 1))
        p_prev, p_cur = p_cur, a * (x * p_cur - b * p_prev)
    end
    return p_cur
end

# ----------------------------------------------------------------------------
# Spheroidal harmonics (s = 0, massive: c² = a²(ω² - μ²))
# ----------------------------------------------------------------------------

"""
    km_spheroidal(l, m, c2; nextra=8, nq=64) -> (Λ, Dict l' => b_l')

Eigenfunction of (L² - c² cos²θ) S = Λ S continuously connected to Y_lm,
S = Σ b_l' Y_l'm with Σ b² = 1 (so ∫ S² dΩ = 1).
"""
function km_spheroidal(l::Int, m::Int, c2; nextra::Int=8, nq::Int=64)
    ls = [lp for lp in abs(m):(l + 2nextra) if iseven(lp - l)]
    x, w = gauss(nq)                       # x = cosθ
    Y = [km_theta_lm(lp, m, acos(xi)) for xi in x, lp in ls]
    n = length(ls)
    A = zeros(ComplexF64, n, n)
    for i in 1:n, j in 1:n
        cc = 2π * sum(w .* Y[:, i] .* Y[:, j] .* x .^ 2)
        A[i, j] = (i == j ? ls[i] * (ls[i] + 1) : 0.0) - c2 * cc
    end
    E = eigen(A)
    i0 = findfirst(==(l), ls)
    k = argmax(abs.(E.vectors[i0, :]))
    v = E.vectors[:, k]
    v ./= sqrt(sum(v .^ 2))
    v .*= sign(real(v[i0]))
    return E.values[k], Dict(ls[i] => v[i] for i in 1:n)
end

# ----------------------------------------------------------------------------
# Leaver series
# ----------------------------------------------------------------------------

"""Parameters and recurrence coefficients of Leaver's series at ω, in the precision of C (ComplexF64 or
Complex{BigFloat}; alph and a are the exact Float64 inputs)."""
struct LeaverCoeffs{C<:Complex}
    ω::C; alph::Float64; a::Float64; m::Int; Λ::C
    q::C; σ::C; χ::C
    c0::C; c1::C; c2::C; c3::C; c4::C
end

function leaver_coeffs(ω::Complex{T}, alph, a, m, Λ) where {T<:AbstractFloat}
    al = T(alph); aa = T(a)
    b = sqrt(1 - aa^2); rp, rm = 1 + b, 1 - b
    Λc = Complex{T}(Λ)
    q = -sqrt((al - ω) * (al + ω)); real(q) > 0 && (q = -q)      # μ² - ω² as a product: no cancellation
    ΩH = aa / (2rp)
    σ = 2rp * (ω - ΩH * m) / (rp - rm)
    χ = (al^2 - 2ω^2) / q
    c0 = 1 - 2im * ω - 2im / b * (ω - aa * m / 2)
    c1 = -4 + 4im * (ω - im * q * (1 + b)) + 4im / b * (ω - aa * m / 2) - 2 * (ω^2 + q^2) / q
    c2 = 3 - 2im * ω - 2 * (q^2 - ω^2) / q - 2im / b * (ω - aa * m / 2)
    c3 = 2im * (ω - im * q)^3 / q + 2 * (ω - im * q)^2 * b + q^2 * aa^2 + 2im * q * aa * m - Λc - 1 - (ω - im * q)^2 / q +
         2q * b + 2im / b * ((ω - im * q)^2 / q + 1) * (ω - aa * m / 2)
    c4 = (ω - im * q)^4 / q^2 + 2im * ω * (ω - im * q)^2 / q - 2im / b * (ω - im * q)^2 / q * (ω - aa * m / 2)
    return LeaverCoeffs{Complex{T}}(ω, Float64(alph), Float64(a), m, Λc, q, σ, χ, c0, c1, c2, c3, c4)
end
leaver_coeffs(ω::Number, alph, a, m, Λ) = leaver_coeffs(ComplexF64(ω), alph, a, m, Λ)

"""The coefficients rounded to Float64 (once, from the precision they were computed in)."""
km_round(L::LeaverCoeffs) = LeaverCoeffs{ComplexF64}(ComplexF64(L.ω), L.alph, L.a, L.m, ComplexF64(L.Λ), ComplexF64(L.q),
                                                     ComplexF64(L.σ), ComplexF64(L.χ), ComplexF64(L.c0), ComplexF64(L.c1),
                                                     ComplexF64(L.c2), ComplexF64(L.c3), ComplexF64(L.c4))

lv_α(L::LeaverCoeffs, n) = n^2 + (L.c0 + 1) * n + L.c0
lv_β(L::LeaverCoeffs, n) = -2n^2 + (L.c1 + 2) * n + L.c3
lv_γ(L::LeaverCoeffs, n) = n^2 + (L.c2 - 3) * n + L.c4

"""Continued-fraction ratios t_n = -b_n/b_{n-1} (n = 1..K) and the eigenvalue residual β_0/α_0 - t_1, all in
Float64 (not accurate enough to locate the eigenvalue of high levels at high spin: see km_eigenfrequency)."""
function leaver_cf(L::LeaverCoeffs{ComplexF64}, K::Int)
    t = zeros(ComplexF64, K + 1)
    tmp = 0.0im
    for i in K:-1:1
        tmp = lv_γ(L, i) / (lv_β(L, i) - lv_α(L, i) * tmp)
        t[i] = tmp
    end
    return t, lv_β(L, 0) / lv_α(L, 0) - t[1]
end

"""Float64 tail of the continued fraction: t_{J+1} from the backward recurrence started with t_{K+1} = 0."""
function leaver_cf_tail(L::LeaverCoeffs{ComplexF64}, K::Int, J::Int)
    tmp = 0.0im
    for i in K:-1:(J + 1)
        tmp = lv_γ(L, i) / (lv_β(L, i) - lv_α(L, i) * tmp)
    end
    return tmp
end

"""
    leaver_cf_head(L, J, tJ1; coeffs=false) -> (D, R, b)

Head of the recurrence in the precision of L, by Miller's algorithm: b_J = 1, b_{J+1} = -tJ1 b_J (tJ1 = t_{J+1}
from the tail), b_{n-1} = -(α_n b_{n+1} + β_n b_n)/γ_n down to b_0. D = α_0 b_1 + β_0 b_0 is the eigenvalue
condition without poles, R = β_0/α_0 + b_1/b_0 = β_0/α_0 - t_1 the same condition normalised (independent of the
scale of b, with poles where b_0 = 0). With `coeffs`, b = (b_0, ..., b_{J+1})/b_0, the series coefficients.
"""
function leaver_cf_head(L::LeaverCoeffs{C}, J::Int, tJ1; coeffs::Bool=false) where {C}
    bn = one(C); bp = -C(tJ1)                                   # b_J, b_{J+1}
    bs = coeffs ? Vector{C}(undef, J + 2) : C[]
    coeffs && (bs[J + 1] = bn; bs[J + 2] = bp)
    for i in J:-1:1
        bm = -(lv_α(L, i) * bp + lv_β(L, i) * bn) / lv_γ(L, i)   # b_{i-1}
        bp = bn; bn = bm
        coeffs && (bs[i] = bn)
    end
    D = L.c0 * bp + L.c3 * bn                                   # α_0 = c0, β_0 = c3
    R = L.c3 / L.c0 + bp / bn
    if coeffs
        b0 = bs[1]
        for k in eachindex(bs)
            bs[k] /= b0
        end
    end
    return D, R, bs
end

"""
    km_cf_eval(ω, alph, a, m, Λ, K, J; coeffs=false) -> (D, R, b, L)

Eigenvalue condition of the Leaver recurrence at ω::Complex{BigFloat} in mixed precision: coefficients and the
head n ≤ J in the current BigFloat precision, the tail J < n ≤ K in Float64 from the same coefficients rounded
once. (All ill-conditioning of the recurrence sits in the first ~n_level terms; Float64 arithmetic beyond n = 100
changes the condition by less than 1e-24.)
"""
function km_cf_eval(ω::Complex{BigFloat}, alph, a, m, Λ, K::Int, J::Int; coeffs::Bool=false)
    L = leaver_coeffs(ω, alph, a, m, Λ)
    J = min(J, K)
    tJ1 = J < K ? leaver_cf_tail(km_round(L), K, J) : 0.0im
    D, R, bs = leaver_cf_head(L, J, tJ1; coeffs=coeffs)
    return D, R, bs, L
end

"""Number of series terms for ~1e-17 relative accuracy (b_n ~ exp(-2 sqrt(2κ(r_+-r_-) n)))."""
function leaver_terms(alph, ω, a)
    κ = real(sqrt(alph^2 - real(ω)^2 + 0im))
    rp, rm = 1 + sqrt(1 - a^2), 1 - sqrt(1 - a^2)
    return Int(clamp(ceil((40 / (2 * sqrt(2κ * (rp - rm))))^2), 2000, 4_000_000))
end

"""Hydrogenic estimate of the eigenfrequency: the real part to O(α⁵), the spin term from the hyperfine expansion,
returned as a complex number with zero imaginary part. (Until 2026-11-21 an imaginary part 1e-9 α^(4l+5) was added
as a seed for the secant; for l ≥ 6 it exceeds the level spacing above α ~ 1.2-1.6 and made the root search fail
there.)"""
function km_omega_guess(n, l, m, alph, a)
    spin = l == 0 ? 0.0 : 2a * m * alph^5 / n^3 / (l * (l + 0.5) * (l + 1))
    ωR = alph * (1 - alph^2 / (2n^2) - alph^4 / (8n^4) + alph^4 / n^4 * (2l - 3n + 1) / (l + 0.5) + spin)
    return complex(ωR)
end

"""Distance (in ω, M = 1) from level n to the nearest overtone with the same (l, m), hydrogenic: the unit in which
eigenfrequency errors are measured."""
function km_level_gap(n::Int, l::Int, alph)
    up = 1 / (2n^2) - 1 / (2(n + 1)^2)
    dn = n - 1 > l ? 1 / (2(n - 1)^2) - 1 / (2n^2) : Inf
    return alph^3 * min(up, dn)
end

"""
    km_eigenfrequency(n, l, m, alph, a; ω0, tol=1e-12, steptol=1e-14, maxit=60, bits=256, J=300, K) -> NamedTuple

Eigenfrequency of level (n, l, m): secant iteration on the pole-free mixed-precision condition D (km_cf_eval)
from ω0 (default: the hydrogenic expansion), the spheroidal eigenvalue recomputed at every step. The root is then
validated: R = β_0/α_0 - t_1 is re-evaluated at the root with 2 × bits, with a head of 3J terms and with 2K terms,
and each difference divided by |dR/dω| (finite difference over 1e-6 of the level spacing) estimates the error of
the root in units of the level spacing (km_level_gap). While an estimate exceeds `tol`, the precision is doubled,
the head tripled or the number of terms doubled, and the root recomputed.

Returns (ω, ωbig, err, converged, bits, J, K, gap): ω rounded to ComplexF64, ωbig at `bits` bits, err the sum of
the error estimates in level spacings (Inf if the iteration did not converge). The secant converges to the root
nearest its start, which at large α or from a poor seed can be another overtone: the caller checks that.
"""
function km_eigenfrequency(n::Int, l::Int, m::Int, alph, a; ω0=nothing, tol=1e-12, steptol=1e-14, maxit::Int=60,
                           bits::Int=256, J::Int=300, K=nothing, maxbits::Int=8192, maxJ::Int=30_000)
    gap = km_level_gap(n, l, alph)
    seed = ω0 === nothing ? km_omega_guess(n, l, m, alph, a) : ComplexF64(ω0)
    Kmax = leaver_terms(alph, seed, a)
    Kc = K === nothing ? min(Kmax, max(4000, 400n)) : Int(K)
    Λof(ω) = km_spheroidal(l, m, a^2 * (ComplexF64(ω)^2 - alph^2))[1]
    ωb = Complex{BigFloat}(seed); err = Inf; conv = false; retried = false
    while true
        ωb, conv, R0, slope, Λ = setprecision(BigFloat, bits) do
            f(ω) = km_cf_eval(ω, alph, a, m, Λof(ω), Kc, J)[1]
            ω1 = Complex{BigFloat}(seed); ω2 = ω1 + 1e-3 * gap
            f1 = f(ω1); f2 = f(ω2); cv = false
            for it in 1:maxit
                (f2 == f1 || !isfinite(abs(f2))) && break
                ω3 = ω2 - f2 * (ω2 - ω1) / (f2 - f1)
                ω1, f1 = ω2, f2
                ω2 = ω3; f2 = f(ω2)
                if abs(ω2 - ω1) < steptol * gap
                    cv = true
                    break
                end
            end
            Λr = Λof(ω2)
            Rr = km_cf_eval(ω2, alph, a, m, Λr, Kc, J)[2]
            Rδ = km_cf_eval(ω2 + 1e-6 * gap, alph, a, m, Λr, Kc, J)[2]
            (ω2, cv, Rr, abs(Rδ - Rr) / 1e-6, Λr)
        end
        eprec = ehead = etail = Inf
        if conv && isfinite(slope) && slope > 0
            Rp = setprecision(() -> km_cf_eval(ωb, alph, a, m, Λ, Kc, J)[2], BigFloat, 2bits)
            Rh = setprecision(() -> km_cf_eval(ωb, alph, a, m, Λ, Kc, min(3J, Kc))[2], BigFloat, bits)
            Rt = Kc < Kmax ? setprecision(() -> km_cf_eval(ωb, alph, a, m, Λ, min(2Kc, Kmax), J)[2], BigFloat, bits) : R0
            eprec = Float64(abs(Rp - R0) / slope); ehead = Float64(abs(Rh - R0) / slope)
            etail = Float64(abs(Rt - R0) / slope)
            err = eprec + ehead + etail + Float64(abs(R0) / slope)
        else
            err = Inf
        end
        if isfinite(err) && eprec <= tol && ehead <= tol && etail <= tol
            break
        elseif !isfinite(err)
            (retried || bits >= maxbits) && break            # a failed iteration is retried once at twice the precision
            retried = true; bits *= 2
        elseif eprec > tol && bits < maxbits
            bits *= 2
        elseif ehead > tol && 3J <= maxJ && 3J <= Kc
            J *= 3
        elseif etail > tol && Kc < Kmax
            Kc = min(2Kc, Kmax)
        else
            break
        end
    end
    return (ω=ComplexF64(ωb), ωbig=ωb, err=err, converged=conv, bits=bits, J=min(J, Kc), K=Kc, gap=gap)
end

struct KerrMode
    n::Int; l::Int; m::Int
    alph::Float64; a::Float64
    ω::ComplexF64
    Λ::ComplexF64
    sph::Dict{Int,ComplexF64}
    L::LeaverCoeffs{ComplexF64}
    b::Vector{ComplexF64}
    norm::Float64            # f is multiplied by this (unit Klein-Gordon norm after kerr_mode_normalize!)
    ode::Any                 # nothing (pure series) or (sol, r0, r_end): ODE continuation for r > r0 (radial overtones);
                             # sol(r) = [h, h'] with h = f e^{-iωr} (unnormalised)
end
KerrMode(n, l, m, alph, a, ω, Λ, sph, L, b, norm) = KerrMode(n, l, m, alph, a, ω, Λ, sph, L, b, norm, nothing)

"""Coefficients A, B of the radial ODE f'' + A f' + B f = 0 of the horizon-regular profile (A = (Δ' - 2iK)/Δ,
B = -[iK' + μ²r² + a²ω² - 2amω + Λ]/Δ), and their r-derivatives."""
@inline function km_ode_coeffs(md::KerrMode, r)
    a = md.a; ω = md.ω; m = md.m
    Δ = r^2 - 2r + a^2; dΔ = 2(r - 1)
    K = (r^2 + a^2) * ω - a * m
    NA = dΔ - 2im * K; dNA = 2 - 4im * r * ω
    NB = -(2im * r * ω + md.alph^2 * r^2 + a^2 * ω^2 - 2a * m * ω + md.Λ); dNB = -(2im * ω + 2md.alph^2 * r)
    A = NA / Δ; B = NB / Δ
    return A, B, (dNA * Δ - NA * dΔ) / Δ^2, (dNB * Δ - NB * dΔ) / Δ^2
end

"""Right-hand side of the ODE of h = f e^{-iωr}, u = [h, h']:
h'' + [2(r-1) - 2i(2ωr - am)]/Δ h' + [-q²r² + 2ω²r - 2iω - Λ]/Δ h = 0 (the f-ODE of km_ode_coeffs with
f = h e^{iωr}; q² = μ² - ω² from md.L.q, which was rounded from extended precision, so the weak binding of high
levels is not lost in μ² - ω²)."""
function km_ode_rhs!(du, u, md::KerrMode, r)
    a = md.a; ω = md.ω
    Δ = r^2 - 2r + a^2
    Ah = (2 * (r - 1) - 2im * (2ω * r - a * md.m)) / Δ
    Ch = (-md.L.q^2 * r^2 + 2ω^2 * r - 2im * ω - md.Λ) / Δ
    du[1] = u[2]
    du[2] = -Ah * u[2] - Ch * u[1]
    return nothing
end

"""Radius where the series is handed over to the ODE (u = 0.2) and where the integration stops (outer turning
point of the hydrogenic orbit plus 6/|q|), for a mode with alph, a, n, q."""
function km_ode_range(alph, a, n, q)
    rp, rm = 1 + sqrt(1 - a^2), 1 - sqrt(1 - a^2)
    u0 = 0.2
    r0 = (rp - u0 * rm) / (1 - u0)
    r_end = 2n^2 / alph^2 + 6 / abs(real(q))
    return r0, max(r_end, 2r0)
end

"""ODE continuation of the series: solve the h-ODE from r0, where (f0, f1) = (f, f') of the series (default: the
Float64 series of md; kerr_mode passes the values summed in BigFloat), to r_end; stored as (sol, r0, r_end)."""
function km_ode_continue(md::KerrMode; rtol=1e-13, f0=nothing, f1=nothing)
    r0, r_end = km_ode_range(md.alph, md.a, md.n, md.L.q)
    if f0 === nothing
        f0, f1, _ = kerr_radial_f(KerrMode(md.n, md.l, md.m, md.alph, md.a, md.ω, md.Λ, md.sph, md.L, md.b, md.norm, nothing), r0)
    end
    E = exp(-im * md.ω * r0)
    u0 = ComplexF64[f0 * E, (f1 - im * md.ω * f0) * E]
    prob = ODEProblem(km_ode_rhs!, u0, (r0, r_end), md)
    sol = solve(prob, Vern9(); reltol=rtol, abstol=rtol * abs(u0[1]) * 1e-3, maxiters=10^8)
    return (sol, r0, r_end)
end

"""
    kerr_mode(n, l, m, alph, a; ω0=km_omega_guess(...), K, tol=1e-12, maxit=60, radial=:auto, bits=256, J=300)
        -> (KerrMode, err)

Quasi-bound state (n, l, m) at μM = alph and spin a: eigenfrequency from km_eigenfrequency (mixed-precision
continued fraction, validated a posteriori), then the series coefficients (the first J in BigFloat, rounded). The
returned mode is normalised to unit Klein-Gordon norm (kerr_mode_norm). radial = :auto continues the series by its
ODE for n - l - 1 ≥ 2 (see the file docstring), :ode always, :series never. `err` is the estimated error of the
eigenfrequency in units of the level spacing km_level_gap (Inf if the iteration did not converge); `tol` is the
target for it. K overrides the number of series terms.
"""
function kerr_mode(n::Int, l::Int, m::Int, alph, a; ω0=nothing, K=nothing, tol=1e-12, maxit=60, normalize=true,
                   radial::Symbol=:auto, bits::Int=256, J::Int=300)
    radial in (:auto, :ode, :series) || error("kerr_mode: radial must be :auto, :ode or :series")
    use_ode = radial == :ode || (radial == :auto && n - l - 1 >= 2)
    eig = km_eigenfrequency(n, l, m, alph, a; ω0=ω0, tol=tol, maxit=maxit, bits=bits, J=J, K=use_ode ? nothing : K)
    eig.converged || @warn "kerr_mode: secant iteration did not converge for ($n,$l,$m) α=$alph a=$a"
    ω = eig.ω
    real(ω) < alph || error("kerr_mode: ω = $ω is not bound (Re ω ≥ μ) for ($n,$l,$m) α=$alph a=$a")
    Λ, sph = km_spheroidal(l, m, a^2 * (ω^2 - alph^2))
    Jh = eig.J
    L, b, f0, f1 = setprecision(BigFloat, eig.bits) do
        LB = leaver_coeffs(eig.ωbig, alph, a, m, Λ)
        Lr = km_round(LB)
        if use_ode
            # only u ≤ 0.2 is summed: the head suffices (0.2^J); the handover data come from the BigFloat coefficients
            _, _, bs = leaver_cf_head(LB, Jh, Jh < eig.K ? leaver_cf_tail(Lr, eig.K, Jh) : 0.0im; coeffs=true)
            r0 = km_ode_range(alph, a, n, Lr.q)[1]
            g0, g1, _ = km_series_eval(LB, bs, BigFloat(r0))
            (Lr, ComplexF64.(bs), ComplexF64(g0), ComplexF64(g1))
        else
            # the whole series (it converges slowly at large r): head coefficients in BigFloat, the rest from the
            # Float64 ratios
            Kt = K === nothing ? leaver_terms(alph, ω, a) : Int(K)
            t, _ = leaver_cf(Lr, Kt)
            Jb = min(Jh, Kt - 1)
            _, _, bs = leaver_cf_head(LB, Jb, t[Jb + 1]; coeffs=true)
            bf = zeros(ComplexF64, Kt + 1)
            bf[1:(Jb + 1)] .= ComplexF64.(bs[1:(Jb + 1)])
            nb = Kt + 1
            for i in (Jb + 1):Kt
                bf[i + 1] = -t[i] * bf[i]
                if abs(bf[i + 1]) < 1e-300
                    nb = i + 1
                    break
                end
            end
            (Lr, bf[1:nb], 0.0im, 0.0im)
        end
    end
    md = KerrMode(n, l, m, alph, a, ω, Λ, sph, L, b, 1.0, nothing)
    if use_ode
        md = KerrMode(n, l, m, alph, a, ω, Λ, sph, L, b, 1.0, km_ode_continue(md; f0=f0, f1=f1))
    end
    if normalize
        N = kerr_mode_norm(md)
        isfinite(N) && N > 0 || error("kerr_mode: Klein-Gordon norm $N of ($n,$l,$m) α=$alph a=$a is not finite and positive")
        md = KerrMode(n, l, m, alph, a, ω, Λ, sph, L, b, 1 / sqrt(N), md.ode)
    end
    return md, eig.err
end

"""f, f', f'' at r (the horizon-regular ingoing radial function, including md.norm). The series is
summed in the precision of `r` (BigFloat r avoids rounding accumulation over the ~10^4-10^5 terms).
With an ODE continuation (md.ode, radial overtones) the series is used only for r <= r0; beyond r_end
the leading asymptotic form scales the value at r_end."""
function kerr_radial_f(md::KerrMode, r)
    if md.ode !== nothing && r > md.ode[2]
        f, f1, f2, _ = kerr_radial_ode_jet(md, Float64(r))
        C = Complex{promote_type(typeof(r), Float64)}
        return C(f), C(f1), C(f2)
    end
    return kerr_radial_f_series(md, r)
end

"""f, f', f'', f''' from the ODE continuation (r > r0); beyond r_end the asymptotic prefactor scales f(r_end)."""
function kerr_radial_ode_jet(md::KerrMode, r::Float64)
    sol, r0, r_end = md.ode
    ω = md.ω
    if r <= r_end
        u = sol(r); E = exp(im * ω * r)
        f, f1 = u[1] * E, (u[2] + im * ω * u[1]) * E          # f = h e^{iωr}
        A, B, dA, dB = km_ode_coeffs(md, r)
        f2 = -A * f1 - B * f
        f3 = -dA * f1 - A * f2 - dB * f - B * f1
    else
        a = md.a; rp, rm = 1 + sqrt(1 - a^2), 1 - sqrt(1 - a^2); L = md.L; m = md.m
        p = im * L.σ + L.χ - 1 - im * (2rm * ω - m * a) / (rp - rm)
        fe = sol(r_end)[1] * exp(im * ω * r_end)
        # P(r)/P(r_end) with P(x) = (x - r_-)^p e^{(q+iω)x}, as one exponent: separately both factors overflow for
        # weakly bound high levels ((x - r_-)^p ~ 1e250 for n = 45, α = 0.03), and Inf/Inf gave NaN
        f = fe * exp(p * (log(r - rm) - log(r_end - rm)) + (L.q + im * ω) * (r - r_end))
        lp = p / (r - rm) + (L.q + im * ω)
        f1 = f * lp
        f2 = f * (lp^2 - p / (r - rm)^2)
        f3 = f * (lp^3 - 3lp * p / (r - rm)^2 + 2p / (r - rm)^3)
    end
    return md.norm * f, md.norm * f1, md.norm * f2, md.norm * f3
end

"""f''' of the mode's horizon-regular radial function at r (analytic on the ODE branch, else one ForwardDiff pass through f'')."""
function kerr_radial_f3(md::KerrMode, r::Real)
    if md.ode !== nothing && r > md.ode[2]
        return kerr_radial_ode_jet(md, Float64(r))[4]
    end
    d = ForwardDiff.derivative(t -> (f2 = kerr_radial_f_series(md, t)[3]; [real(f2), imag(f2)]), r)
    return complex(d[1], d[2])
end

kerr_radial_f_series(md::KerrMode, r) = md.norm .* km_series_eval(md.L, md.b, r)

"""f, f', f'' of the Leaver series with coefficients b (b[1] = b_0) at r, unnormalised, in the precision of r and
of the coefficients (BigFloat L, b and r for the handover data of the ODE continuation)."""
function km_series_eval(L::LeaverCoeffs, b::AbstractVector, r)
    a = oftype(float(r), L.a); m = L.m; ω = L.ω
    sq = sqrt(1 - a^2); rp, rm = 1 + sq, 1 - sq
    p = im * L.σ + L.χ - 1 - im * (2rm * ω - m * a) / (rp - rm)
    u = (r - rp) / (r - rm); du = (rp - rm) / (r - rm)^2; d2u = -2 * (rp - rm) / (r - rm)^3
    C = promote_type(Complex{promote_type(typeof(r), Float64)}, eltype(b))
    S0 = zero(C); S1 = zero(C); S2 = zero(C)
    for k in length(b):-1:1
        S2 = S2 * u + 2 * S1
        S1 = S1 * u + S0
        S0 = S0 * u + b[k]
    end
    lp = p / (r - rm) + (L.q + im * ω)
    P = exp(p * log(r - rm) + (L.q + im * ω) * r)
    P1 = P * lp; P2 = P * (lp^2 - p / (r - rm)^2)
    f = P * S0
    f1 = P1 * S0 + P * S1 * du
    f2 = P2 * S0 + 2 * P1 * S1 * du + P * (S2 * du^2 + S1 * d2u)
    return f, f1, f2
end

"""r# = ∫ a/Δ dr = a/(r_+-r_-) ln((r-r_+)/(r-r_-))."""
function kerr_rsharp(r, a)
    a == 0 && return zero(r)
    rp, rm = 1 + sqrt(1 - a^2), 1 - sqrt(1 - a^2)
    return a / (rp - rm) * log((r - rp) / (r - rm))
end

"""
Boyer-Lindquist radial function R and R' at r, ψ = R(r) S(θ) e^{-iωt+imφ}
(R = f e^{-iωr* + imr#}; oscillates as e^{-ikr*} at the horizon).
"""
function kerr_radial_BL(md::KerrMode, r)
    f, f1, _ = kerr_radial_f(md, r)
    a = md.a; ω = md.ω
    Δ = r^2 - 2r + a^2
    K = (r^2 + a^2) * ω - a * md.m
    ph = exp(-im * ω * kerr_rstar(r, a) + im * md.m * kerr_rsharp(r, a))
    return ph * f, ph * (f1 - im * K / Δ * f)
end

"""S(θ) and dS/dθ of the mode's spheroidal harmonic (∫ S² dΩ = 1)."""
function kerr_angular(md::KerrMode, θ)
    S = zero(ComplexF64); dS = zero(ComplexF64)
    for (lp, b) in md.sph
        S += b * km_theta_lm(lp, md.m, θ)
        dS += b * ForwardDiff.derivative(t -> km_theta_lm(lp, md.m, t), θ)
    end
    return S, dS
end

"""
    kerr_mode_norm(md; nr=400, nθ=48, nodes=nothing) -> N

Klein-Gordon norm of ψ/√(2ω) on a v = const slice of ingoing Kerr coordinates
(√-g = Σ sinθ, g^{vv} = a² sin²θ/Σ, g^{vr} = (r²+a²)/Σ, g^{vφ} = a/Σ):

    N = (2π/2ω) ∫dr ∫dθ sinθ |S|² [ 2(r²+a²) Im(f* f') + (2ma - 2ω a² sin²θ) |f|² ],

which tends to ∫|ψ|² d³x in the non-relativistic limit. The r integral uses
Gauss-Legendre panels (km_norm_nodes): for a pure series nr nodes in ln(r - r_+) out to
where the mode has decayed; with an ODE continuation the panels are the steps of the ODE
solution, which resolve the radial nodes of the overtone whatever their number.
`nodes` = (r, weights) replaces the radial quadrature (tests).
"""
function kerr_mode_norm(md::KerrMode; nr::Int=400, nθ::Int=48, nodes=nothing)
    a = md.a; ω = real(md.ω); m = md.m
    rp = 1 + sqrt(1 - a^2)
    xθ, wθ = gauss(nθ)
    θs = π / 2 .* (xθ .+ 1); wθs = π / 2 .* wθ
    ang = sum(wθs .* sin.(θs) .* abs2.(first.(kerr_angular.(Ref(md), θs))))
    ang2 = sum(wθs .* sin.(θs) .^ 3 .* abs2.(first.(kerr_angular.(Ref(md), θs))))
    rn, rw = nodes === nothing ? km_norm_nodes(md, nr) : nodes
    acc = 0.0
    for (r, w) in zip(rn, rw)
        f, f1, _ = kerr_radial_f(md, r)
        acc += w * (2 * (r^2 + a^2) * imag(conj(f) * f1) * ang + (2m * a * ang - 2ω * a^2 * ang2) * abs2(f))
    end
    return 2π * acc / (2ω)
end

"""
    km_radial_nodes(md; per_step=8) -> number of radial nodes, or -1 without an ODE continuation

Zeros of the Boyer-Lindquist radial function R = f e^{-iωr* + imr#} of an ODE-continued mode between half its
hydrogenic inner turning point and the end of the integration: sign changes of Re(R conj(R_peak)) on the ODE steps,
each subdivided per_step times (for a quasi-bound mode R is real up to a constant phase where it is not
evanescent). A level (n, l, m) has n - l - 1 of them, whatever its frequency shift: an independent identification
of the overtone. Not meaningful for l = 0 (no barrier: the ingoing phase near the horizon adds sign changes) or for
resonances as broad as the level spacing.
"""
function km_radial_nodes(md::KerrMode; per_step::Int=8)
    md.ode === nothing && return -1
    sol, r0, r_end = md.ode
    n, l = md.n, md.l
    rin = n^2 / md.alph^2 * (1 - sqrt(max(0.0, 1 - l * (l + 1) / n^2)))
    rlo = max(r0, rin / 2)
    ωr = real(md.ω)
    Rs = ComplexF64[]
    ts = sol.t
    for k in 1:(length(ts) - 1), j in 0:(per_step - 1)
        r = ts[k] + (ts[k + 1] - ts[k]) * j / per_step
        r < rlo && continue
        f = kerr_radial_ode_jet(md, r)[1]
        push!(Rs, f * exp(-im * ωr * kerr_rstar(r, md.a) + im * md.m * kerr_rsharp(r, md.a)))
    end
    isempty(Rs) && return -1
    ref = Rs[argmax(abs.(Rs))]
    g = real.(Rs .* conj(ref))
    cnt = 0; last = 0.0
    for v in g
        v == 0 && continue
        last != 0 && sign(v) != sign(last) && (cnt += 1)
        last = v
    end
    return cnt
end

"""Radial quadrature nodes and weights of the norm integral. Pure series: km_log_panels with nr nodes out to
r_+ + 80 n²/α². ODE continuation: panels in ln(r - r_+) up to the handover radius r0, then two 8-point panels per
step of the ODE solution up to r_end, then the asymptotic tail over 25 decay lengths."""
function km_norm_nodes(md::KerrMode, nr::Int)
    rp = 1 + sqrt(1 - md.a^2)
    md.ode === nothing && return km_log_panels(rp, rp + 80 * md.n^2 / md.alph^2, nr)
    sol, r0, r_end = md.ode
    x, w = gauss(8)
    rn = Float64[]; rw = Float64[]
    function panel(lo, hi)
        append!(rn, (lo + hi) / 2 .+ (hi - lo) / 2 .* x)
        append!(rw, (hi - lo) / 2 .* w)
    end
    edges = rp .+ exp.(range(log(1e-6), log(r0 - rp), length=25))
    edges[end] = r0
    for k in 1:(length(edges) - 1)
        panel(edges[k], edges[k + 1])
    end
    ts = sol.t
    for k in 1:(length(ts) - 1)
        mid = (ts[k] + ts[k + 1]) / 2
        panel(ts[k], mid); panel(mid, ts[k + 1])
    end
    tail = 25 / abs(real(md.L.q))
    for k in 1:12
        panel(r_end + (k - 1) * tail / 12, r_end + k * tail / 12)
    end
    return rn, rw
end

"""Gauss-Legendre nodes and weights on [r_lo, r_hi], uniform panels in ln(r - r_lo) from r_lo + 1e-6."""
function km_log_panels(r_lo, r_hi, n::Int; per_panel::Int=8, eps=1e-6)
    np = max(1, n ÷ per_panel)
    x, w = gauss(per_panel)
    edges = exp.(range(log(eps), log(r_hi - r_lo), length=np + 1))
    rn = Float64[]; rw = Float64[]
    for k in 1:np
        lo, hi = edges[k], edges[k + 1]
        append!(rn, r_lo .+ (lo + hi) / 2 .+ (hi - lo) / 2 .* x)
        append!(rw, (hi - lo) / 2 .* w)
    end
    return rn, rw
end
