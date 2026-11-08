"""
    gw_kerr.jl

Relativistic gravitational-wave emission of a scalar cloud on Kerr: annihilation
a × b -> graviton (ω̃ = ω_a + ω_b, m̃ = m_a + m_b) and transition a -> b + graviton
(ω̃ = ω_a - ω_b, m̃ = m_a - m_b), from the Teukolsky amplitude at infinity
(teukolsky.jl) sourced by the quasi-bound modes of kerr_modes.jl. G = c = M = 1,
μ = α; the levels are evaluated at Re ω_j.

Field: Φ = Σ_j √(N_j/(2ω_j)) (ψ_j + ψ_j*), ψ_j = R_j(r) S_j(θ) e^{-iω_j t + im_j φ}
of unit Klein-Gordon norm. Since n and m̄ are null with n·m̄ = 0, the source
components are bilinear in first derivatives, T_nn = (n·∂Φ)², T_nm̄ = (n·∂Φ)(m̄·∂Φ),
T_m̄m̄ = (m̄·∂Φ)², and the component with the process's (ω̃, m̃) is

    T̂_nn = c n_a n_b,   T̂_nm̄ = (c/2)(n_a m_b + n_b m_a),   T̂_m̄m̄ = c m_a m_b,

n_j = n·∂ψ_j, m_j = m̄·∂ψ_j (ψ_b -> ψ_b* for transitions), c = √(N_aN_b/(ω_aω_b))
for a ≠ b and N_a/(2ω_a) for a = b. With N_a = N_b = 1 the power to infinity is
P = Σ_ℓ |Z_ℓ|²/(2πω̃²) (both helicities, by the parity of the levels) and the
rate in the solver's convention (du_a/dτ = -Γ u_a u_b, τ = μt, u = N/(GM²)) is

    Γ(α, a) = P/(α ω̃)      (G = c = M = 1 units of P and ω̃).
"""

"""Tetrad components of the pair on the grid: (Tnn, Tnm, Tmm, ω̃, m̃)."""
function gw_kerr_tetrad(ma::KerrMode, mb::KerrMode, kind::Symbol, rn, θn)
    kind in (:ann, :tr) || error("gw_kerr_tetrad: kind must be :ann or :tr")
    a = ma.a
    a == mb.a && ma.alph == mb.alph || error("gw_kerr_tetrad: the two modes must share (α, a)")
    ωa, ωb = real(ma.ω), real(mb.ω)
    same = (ma.n, ma.l, ma.m) == (mb.n, mb.l, mb.m)
    conjb = kind == :tr
    conjb && same && error("gw_kerr_tetrad: a transition needs two distinct levels")
    ω̃ = conjb ? ωa - ωb : ωa + ωb
    m̃ = conjb ? ma.m - mb.m : ma.m + mb.m
    c = same ? 1 / (2ωa) : 1 / sqrt(ωa * ωb)
    # radial factors: n·∂ψ = -[ΔR' + iKR] S/(2Σ) = -Δ f' e^{-iωr*+imr#} S/(2Σ) (the iKR cancels exactly
    # against the phase derivative, so use f' directly),  m̄·∂ψ = (ρ/√2)[S' + (m/sinθ - aω sinθ) S] R
    C = Complex{promote_type(eltype(rn), Float64)}
    nr_a = zeros(C, length(rn)); nr_b = similar(nr_a); Ra = similar(nr_a); Rb = similar(nr_a)
    for (i, r) in enumerate(rn)
        Δ = r^2 - 2r + a^2
        for (md, nr_, R_) in ((ma, nr_a, Ra), (mb, nr_b, Rb))
            f, f1, _ = kerr_radial_f(md, r)
            ph = exp(-im * md.ω * kerr_rstar(r, a) + im * md.m * kerr_rsharp(r, a))
            nr_[i] = -Δ * f1 * ph; R_[i] = f * ph
        end
        if conjb
            nr_b[i] = conj(nr_b[i]); Rb[i] = conj(Rb[i])
        end
    end
    Sa = similar(θn, ComplexF64); La = similar(Sa); Sb = similar(Sa); Lb = similar(Sa)
    for (j, θ) in enumerate(θn)
        S, dS = kerr_angular(ma, θ)
        Sa[j] = S; La[j] = dS + (ma.m / sin(θ) - a * ωa * sin(θ)) * S
        S, dS = kerr_angular(mb, θ)
        if conjb
            Sb[j] = conj(S); Lb[j] = conj(dS - (mb.m / sin(θ)) * S + a * ωb * sin(θ) * S)   # conj(m·∂ψ_b)
        else
            Sb[j] = S; Lb[j] = dS + (mb.m / sin(θ) - a * ωb * sin(θ)) * S
        end
    end
    Tnn = zeros(C, length(rn), length(θn)); Tnm = similar(Tnn); Tmm = similar(Tnn)
    for (i, r) in enumerate(rn), (j, θ) in enumerate(θn)
        Σ = r^2 + a^2 * cos(θ)^2
        ρ = 1 / (r - im * a * cos(θ))
        na = nr_a[i] * Sa[j] / (2Σ); nb = nr_b[i] * Sb[j] / (2Σ)
        mA = ρ / sqrt(2) * La[j] * Ra[i]; mB = ρ / sqrt(2) * Lb[j] * Rb[i]
        Tnn[i, j] = c * na * nb
        Tnm[i, j] = c / 2 * (na * mB + nb * mA)
        Tmm[i, j] = c * mA * mB
    end
    return Tnn, Tnm, Tmm, ω̃, m̃
end

"""
    gw_kerr_grid(ma, mb, ω̃; eps=1e-6, ppr=8, ppw=8, decay=40, nθ=48, per_panel=8)

Quadrature grid for the pair: Gauss-Legendre panels in r from r_+ + eps out to where
R_a R_b has decayed to e^{-decay} (R^in grows as r^3 in the wave zone and annihilation
amplitudes are ~1e-9 of their bulk partial sums, so the tail must be below ~1e-12 of the
peak: decay = 45 is the minimum there, 60 the default), each panel of width
min((r - r_+)/ppr, (2π/ω̃)/ppw), so
that both the radial structure (which lives on the scale of r - r_+) and the oscillation
e^{iω̃ r*} of R^in (annihilations have ω̃ ≈ 2μ, with the cloud in the wave zone) are
resolved with spectral accuracy. The by-parts (kernel) form of the amplitude differs from
the direct form by total derivatives of bulk-sized terms, so the panel quadrature error
sets the floor of the strongly cancelling annihilation integrals; ppr matters for that.
Gauss-Legendre in θ.
"""
function gw_kerr_grid(ma::KerrMode, mb::KerrMode, ω̃; eps=1e-6, ppr::Real=8.0, ppw::Real=8.0, decay::Real=60.0,
                      nθ::Int=48, per_panel::Int=8)
    rp = 1 + sqrt(1 - ma.a^2)
    qsum = abs(real(ma.L.q)) + abs(real(mb.L.q))
    xmax = decay / qsum
    hosc = 2π / ω̃ / ppw
    x, w = gauss(per_panel)
    rn = Float64[]; rw = Float64[]
    lo = eps
    while lo < xmax
        h = min(lo / ppr, hosc)
        hi = min(lo + h, xmax)
        append!(rn, rp .+ (lo + hi) / 2 .+ (hi - lo) / 2 .* x)
        append!(rw, (hi - lo) / 2 .* w)
        lo = hi
    end
    xθ, wθ = gauss(nθ)
    return rn, rw, π / 2 .* (xθ .+ 1), π / 2 .* wθ
end
gw_kerr_grid(ma::KerrMode, mb::KerrMode; kwargs...) = gw_kerr_grid(ma, mb, abs(real(ma.ω) - real(mb.ω)) > 0 ? abs(real(ma.ω) - real(mb.ω)) : real(ma.ω) + real(mb.ω); kwargs...)

"""
    gw_kerr_rate(ma, mb, kind; lmax=12, ltol=1e-8, nθ=48, eps=1e-6, ppr=8, ppw=8, decay=60,
                 T=Float64, verbose=false, teuk_kwargs...) -> (Γ, P, Dict ℓ => |Z_ℓ|², ω̃, m̃)

Rate of the annihilation (kind = :ann) or transition a -> b (kind = :tr,
needs Re ω_a > Re ω_b) in the solver's convention, summing ℓ from max(2, |m̃|)
until two consecutive ℓ contribute less than ltol of the total.
"""
function gw_kerr_rate(ma::KerrMode, mb::KerrMode, kind::Symbol; lmax::Int=12, ltol=1e-8, nθ::Int=48, eps=1e-6,
                      ppr::Real=8.0, ppw::Real=8.0, decay::Real=60.0, T=Float64, verbose::Bool=false, teuk_kwargs...)
    ω̃0 = kind == :tr ? real(ma.ω) - real(mb.ω) : real(ma.ω) + real(mb.ω)
    ω̃0 > 0 || error("gw_kerr_rate: ω̃ = $ω̃0 ≤ 0 (order a transition as a -> b with Re ω_a > Re ω_b)")
    rn, rw, θn, θw = gw_kerr_grid(ma, mb, ω̃0; eps=eps, ppr=ppr, ppw=ppw, decay=decay, nθ=nθ)
    rn = T.(rn); rw = T.(rw)          # the r-dependent part of the integral is done in precision T
    Tnn, Tnm, Tmm, ω̃, m̃ = gw_kerr_tetrad(ma, mb, kind, rn, θn)
    modes = Dict{Int,Float64}()
    P = 0.0
    nsmall = 0
    # R^in must be available over the whole source grid and is stepped through its nodes
    kw = merge((r_in_max=1.05 * rn[end], tstops_in=rn, T=T), values(teuk_kwargs))
    for l in max(2, abs(m̃)):lmax
        tr = teuk_radial(-2, l, m̃, T(ω̃), T(ma.a); kw...)
        Z = teuk_Z_grid(tr, Tnn, Tnm, Tmm, rn, rw, θn, θw)
        modes[l] = Float64(abs2(Z) / (2π * ω̃^2))
        P += modes[l]
        verbose && println("  ℓ = $l: |Z|²/(2πω̃²) = $(modes[l])   W_spread = $(tr.W_spread)")
        nsmall = modes[l] < ltol * P ? nsmall + 1 : 0
        nsmall >= 2 && break
    end
    return P / (ma.alph * Float64(ω̃)), P, modes, Float64(ω̃), m̃
end
