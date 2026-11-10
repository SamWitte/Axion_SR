"""
    kerr_modes.jl

Quasi-bound states of a massive scalar on Kerr for the relativistic GW rates
(G = c = M = 1, μ = α). The radial function is Leaver's series (Dolan 2007), as
in solve_sr_rates.jl / find_im_part, in the horizon-regular ingoing form

    ψ = f(r) S(θ) e^{-iωv + imφ̃},   f = R_BL e^{iωr* - imr#},
    R_BL = (r-r_+)^{-iσ} (r-r_-)^{iσ+χ-1} e^{qr} Σ_n b_n u^n,   u = (r-r_+)/(r-r_-),

with v = t + r*, φ̃ = φ + r#, dr#/dr = a/Δ. The minimal solution b_n comes from
the continued-fraction ratios (backward recurrence) and the complex eigenvalue
ω from a secant iteration on the continued fraction, with the exact spheroidal
eigenvalue Λ of (L² - c² cos²θ) S = Λ S, c² = a²(ω² - μ²), from diagonalising
in the Y_lm basis (S = Σ b_l' Y_l'm, so θ-derivatives are exact).

Normalisation: N_j quanta correspond to the field √(N_j/(2ω_j)) (ψ_j + ψ_j*)
with ψ_j of unit Klein-Gordon norm on a v = const slice (finite at the
horizon; the Boyer-Lindquist slice diverges there for real ω). In the
non-relativistic limit this is ∫|ψ|² d³x = 1.

Radial overtones (N = n - l - 1 ≥ 2, kerr_mode radial=:auto): the series Σ b_n uⁿ
behaves as (1-u)^N g(u) with large alternating b_n, so at large r (u → 1) it cancels
to ~r^{-N} of its terms: in Float64 the sum is wrong beyond a few overtones, and since
the Klein-Gordon norm is set there, the normalised profile came out too small by many
orders of magnitude in the inner region that overlaps compact levels. For N ≥ 2 the
series is used only up to u = 0.2, and f is continued outward by integrating its ODE,

    f'' + (Δ' - 2iK)/Δ f' - [iK' + μ²r² + a²ω² - 2amω + Λ]/Δ f = 0,   K = (r²+a²)ω - am,

which is stable up to the outer turning point (the regular solution dominates or
oscillates there) and is stopped 6/|q| beyond it; further out f follows the leading
asymptotic form (r-r_-)^{iσ+χ-1} e^{(q+iω)r}, where it is below e^{-6} of its peak.
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

struct LeaverCoeffs
    ω::ComplexF64; alph::Float64; a::Float64; m::Int; Λ::ComplexF64
    q::ComplexF64; σ::ComplexF64; χ::ComplexF64
    c0::ComplexF64; c1::ComplexF64; c2::ComplexF64; c3::ComplexF64; c4::ComplexF64
end

function leaver_coeffs(ω, alph, a, m, Λ)
    rp, rm = 1 + sqrt(1 - a^2), 1 - sqrt(1 - a^2); b = sqrt(1 - a^2)
    q = -sqrt(alph^2 - ω^2 + 0im); real(q) > 0 && (q = -q)
    ΩH = a / (2rp)
    σ = 2rp * (ω - ΩH * m) / (rp - rm)
    χ = (alph^2 - 2ω^2) / q
    c0 = 1 - 2im * ω - 2im / b * (ω - a * m / 2)
    c1 = -4 + 4im * (ω - im * q * (1 + b)) + 4im / b * (ω - a * m / 2) - 2 * (ω^2 + q^2) / q
    c2 = 3 - 2im * ω - 2 * (q^2 - ω^2) / q - 2im / b * (ω - a * m / 2)
    c3 = 2im * (ω - im * q)^3 / q + 2 * (ω - im * q)^2 * b + q^2 * a^2 + 2im * q * a * m - Λ - 1 - (ω - im * q)^2 / q +
         2q * b + 2im / b * ((ω - im * q)^2 / q + 1) * (ω - a * m / 2)
    c4 = (ω - im * q)^4 / q^2 + 2im * ω * (ω - im * q)^2 / q - 2im / b * (ω - im * q)^2 / q * (ω - a * m / 2)
    return LeaverCoeffs(ω, alph, a, m, Λ, q, σ, χ, c0, c1, c2, c3, c4)
end

lv_α(L::LeaverCoeffs, n) = n^2 + (L.c0 + 1) * n + L.c0
lv_β(L::LeaverCoeffs, n) = -2n^2 + (L.c1 + 2) * n + L.c3
lv_γ(L::LeaverCoeffs, n) = n^2 + (L.c2 - 3) * n + L.c4

"""Continued-fraction ratios t_n = -b_n/b_{n-1} (n = 1..K) and the eigenvalue residual β_0/α_0 - t_1."""
function leaver_cf(L::LeaverCoeffs, K::Int)
    t = zeros(ComplexF64, K + 1)
    tmp = 0.0im
    for i in K:-1:1
        tmp = lv_γ(L, i) / (lv_β(L, i) - lv_α(L, i) * tmp)
        t[i] = tmp
    end
    return t, lv_β(L, 0) / lv_α(L, 0) - t[1]
end

"""Number of series terms for ~1e-17 relative accuracy (b_n ~ exp(-2 sqrt(2κ(r_+-r_-) n)))."""
function leaver_terms(alph, ω, a)
    κ = real(sqrt(alph^2 - real(ω)^2 + 0im))
    rp, rm = 1 + sqrt(1 - a^2), 1 - sqrt(1 - a^2)
    return Int(clamp(ceil((40 / (2 * sqrt(2κ * (rp - rm))))^2), 2000, 4_000_000))
end

"""Hydrogenic estimate of the eigenfrequency (real part to O(α⁵), the spin term from the hyperfine
expansion; a small positive imaginary part as a seed for the secant iteration)."""
function km_omega_guess(n, l, m, alph, a)
    ωR = alph * (1 - alph^2 / (2n^2) - alph^4 / (8n^4) + alph^4 / n^4 * (2l - 3n + 1) / (l + 0.5) +
                 2a * m * alph^5 / n^3 / (l * (l + 0.5) * (l + 1)))
    return ωR + 1e-9im * alph^(4l + 5)
end

struct KerrMode
    n::Int; l::Int; m::Int
    alph::Float64; a::Float64
    ω::ComplexF64
    Λ::ComplexF64
    sph::Dict{Int,ComplexF64}
    L::LeaverCoeffs
    b::Vector{ComplexF64}
    norm::Float64            # f is multiplied by this (unit Klein-Gordon norm after kerr_mode_normalize!)
    ode::Any                 # nothing (pure series) or (sol, r0, r_end): ODE continuation for r > r0 (radial overtones)
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

function km_ode_rhs!(du, u, md::KerrMode, r)
    A, B, _, _ = km_ode_coeffs(md, r)
    du[1] = u[2]
    du[2] = -A * u[2] - B * u[1]
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

"""ODE continuation of the series: solve from r0 (series initial data) to r_end; stored as (sol, r0, r_end)."""
function km_ode_continue(md::KerrMode; rtol=1e-13)
    r0, r_end = km_ode_range(md.alph, md.a, md.n, md.L.q)
    f0, f1, _ = kerr_radial_f(KerrMode(md.n, md.l, md.m, md.alph, md.a, md.ω, md.Λ, md.sph, md.L, md.b, md.norm, nothing), r0)
    prob = ODEProblem(km_ode_rhs!, ComplexF64[f0, f1], (r0, r_end), md)
    sol = solve(prob, Vern9(); reltol=rtol, abstol=rtol * abs(f0) * 1e-3, maxiters=10^8)
    return (sol, r0, r_end)
end

"""
    kerr_mode(n, l, m, alph, a; ω0=km_omega_guess(...), K, tol=1e-13, maxit=60, radial=:auto) -> (KerrMode, residual)

Quasi-bound state (n, l, m) at μM = alph and spin a: secant iteration on the
continued fraction from ω0, then the series coefficients. The returned mode is
normalised to unit Klein-Gordon norm (kerr_mode_norm). radial = :auto continues
the series by its ODE for n - l - 1 ≥ 2 (see the file docstring), :ode always, :series never.
"""
function kerr_mode(n::Int, l::Int, m::Int, alph, a; ω0=nothing, K=nothing, tol=1e-13, maxit=60, normalize=true,
                   radial::Symbol=:auto)
    ω0 === nothing && (ω0 = km_omega_guess(n, l, m, alph, a))
    Kt = K === nothing ? leaver_terms(alph, ω0, a) : K
    function resid(ω)
        Λ, _ = km_spheroidal(l, m, a^2 * (ω^2 - alph^2))
        _, r = leaver_cf(leaver_coeffs(ω, alph, a, m, Λ), Kt)
        return r
    end
    ω1 = complex(ω0); ω2 = ω1 * (1 + 1e-7) + 1e-12im
    f1 = resid(ω1); f2 = resid(ω2)
    converged = false
    for it in 1:maxit
        ω3 = ω2 - f2 * (ω2 - ω1) / (f2 - f1)
        ω1, f1 = ω2, f2
        ω2 = ω3; f2 = resid(ω2)
        if abs(ω2 - ω1) < tol * abs(ω2)
            converged = true
            break
        end
    end
    converged || @warn "kerr_mode: secant iteration did not converge for ($n,$l,$m) α=$alph a=$a; |res| = $(abs(f2))"
    ω = ω2
    real(ω) < alph || error("kerr_mode: ω = $ω is not bound (Re ω ≥ μ) for ($n,$l,$m) α=$alph a=$a")
    Λ, sph = km_spheroidal(l, m, a^2 * (ω^2 - alph^2))
    L = leaver_coeffs(ω, alph, a, m, Λ)
    t, res = leaver_cf(L, Kt)
    b = zeros(ComplexF64, Kt + 1); b[1] = 1
    nb = Kt + 1
    for i in 1:Kt
        b[i + 1] = -t[i] * b[i]
        if abs(b[i + 1]) < 1e-300
            nb = i + 1
            break
        end
    end
    md = KerrMode(n, l, m, alph, a, ω, Λ, sph, L, b[1:nb], 1.0, nothing)
    use_ode = radial == :ode || (radial == :auto && n - l - 1 >= 2)
    radial in (:auto, :ode, :series) || error("kerr_mode: radial must be :auto, :ode or :series")
    if use_ode
        md = KerrMode(n, l, m, alph, a, ω, Λ, sph, L, b[1:nb], 1.0, km_ode_continue(md))
    end
    if normalize
        md = KerrMode(n, l, m, alph, a, ω, Λ, sph, L, b[1:nb], 1 / sqrt(kerr_mode_norm(md)), md.ode)
    end
    return md, res
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
    if r <= r_end
        u = sol(r)
        f, f1 = u[1], u[2]
        A, B, dA, dB = km_ode_coeffs(md, r)
        f2 = -A * f1 - B * f
        f3 = -dA * f1 - A * f2 - dB * f - B * f1
    else
        a = md.a; rp, rm = 1 + sqrt(1 - a^2), 1 - sqrt(1 - a^2); L = md.L; ω = md.ω; m = md.m
        p = im * L.σ + L.χ - 1 - im * (2rm * ω - m * a) / (rp - rm)
        fe = sol(r_end)[1]
        P(x) = exp(p * log(x - rm) + (L.q + im * ω) * x)
        f = fe * P(r) / P(r_end)
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

function kerr_radial_f_series(md::KerrMode, r)
    a = md.a; rp, rm = 1 + sqrt(1 - a^2), 1 - sqrt(1 - a^2); L = md.L; ω = md.ω; m = md.m
    p = im * L.σ + L.χ - 1 - im * (2rm * ω - m * a) / (rp - rm)
    u = (r - rp) / (r - rm); du = (rp - rm) / (r - rm)^2; d2u = -2 * (rp - rm) / (r - rm)^3
    C = Complex{promote_type(typeof(r), Float64)}
    S0 = zero(C); S1 = zero(C); S2 = zero(C)
    for k in length(md.b):-1:1
        S2 = S2 * u + 2 * S1
        S1 = S1 * u + S0
        S0 = S0 * u + md.b[k]
    end
    lp = p / (r - rm) + (L.q + im * ω)
    P = exp(p * log(r - rm) + (L.q + im * ω) * r)
    P1 = P * lp; P2 = P * (lp^2 - p / (r - rm)^2)
    f = P * S0
    f1 = P1 * S0 + P * S1 * du
    f2 = P2 * S0 + 2 * P1 * S1 * du + P * (S2 * du^2 + S1 * d2u)
    return md.norm * f, md.norm * f1, md.norm * f2
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
    kerr_mode_norm(md; nr=400, nθ=48) -> N

Klein-Gordon norm of ψ/√(2ω) on a v = const slice of ingoing Kerr coordinates
(√-g = Σ sinθ, g^{vv} = a² sin²θ/Σ, g^{vr} = (r²+a²)/Σ, g^{vφ} = a/Σ):

    N = (2π/2ω) ∫dr ∫dθ sinθ |S|² [ 2(r²+a²) Im(f* f') + (2ma - 2ω a² sin²θ) |f|² ],

which tends to ∫|ψ|² d³x in the non-relativistic limit. The r integral uses
Gauss-Legendre panels in ln(r - r_+) out to where the mode has decayed.
"""
function kerr_mode_norm(md::KerrMode; nr::Int=400, nθ::Int=48)
    a = md.a; ω = real(md.ω); m = md.m
    rp = 1 + sqrt(1 - a^2)
    xθ, wθ = gauss(nθ)
    θs = π / 2 .* (xθ .+ 1); wθs = π / 2 .* wθ
    ang = sum(wθs .* sin.(θs) .* abs2.(first.(kerr_angular.(Ref(md), θs))))
    ang2 = sum(wθs .* sin.(θs) .^ 3 .* abs2.(first.(kerr_angular.(Ref(md), θs))))
    rmax = rp + 80 * md.n^2 / md.alph^2
    rn, rw = km_log_panels(rp, rmax, nr)
    acc = 0.0
    for (r, w) in zip(rn, rw)
        f, f1, _ = kerr_radial_f(md, r)
        acc += w * (2 * (r^2 + a^2) * imag(conj(f) * f1) * ang + (2m * a * ang - 2ω * a^2 * ang2) * abs2(f))
    end
    return 2π * acc / (2ω)
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
