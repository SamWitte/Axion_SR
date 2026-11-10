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
    # against the phase derivative, so use f' directly),  m̄·∂ψ = (ρ/√2)[S' + (m/sinθ - aω sinθ) S] R.
    # The source is quasi-stationary at Re ω: the profile is the horizon-regular f(r) with the phase
    # e^{-i Re ω r* + imr#}. (With the complex ω in the phase, R = f e^{-iωr*+imr#} is the Boyer-Lindquist
    # eigenfunction, which carries e^{Im ω r*}: harmless for Z^∞, where R^in ∝ Δ² suppresses the horizon,
    # but its derivatives ~ Im ω/(r-r_+)² spoil the horizon amplitude.)
    C = Complex{promote_type(eltype(rn), Float64)}
    nr_a = zeros(C, length(rn)); nr_b = similar(nr_a); Ra = similar(nr_a); Rb = similar(nr_a)
    for (i, r) in enumerate(rn)
        Δ = r^2 - 2r + a^2
        for (md, nr_, R_) in ((ma, nr_a, Ra), (mb, nr_b, Rb))
            f, f1, _ = kerr_radial_f(md, r)
            ph = exp(-im * real(md.ω) * kerr_rstar(r, a) + im * md.m * kerr_rsharp(r, a))
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

# kerr_radial_f3 (f''' of the mode) is in kerr_modes.jl.

"""
(S, S', S'', L, L', L'') of mode md at θ, L = S' + q S with q = sgn (m/sinθ - aω sinθ);
S'' and S''' from the spheroidal equation (exact). sgn = -1 with conj = true gives the
factor of m̄·∂ψ* (transitions): conj(S' - (m/sinθ - aω sinθ) S) and its derivatives.
"""
function kerr_angular_jet(md::KerrMode, θ; conjugate::Bool=false)
    a = md.a; ω = real(md.ω); m = md.m
    S, S1 = kerr_angular(md, θ)
    c2 = a^2 * (md.ω^2 - md.alph^2)
    sθ, cθ = sin(θ), cos(θ)
    V = c2 * cθ^2 - m^2 / sθ^2 + md.Λ
    dV = -2c2 * cθ * sθ + 2m^2 * cθ / sθ^3
    S2 = -cθ / sθ * S1 - V * S
    S3 = S1 / sθ^2 - cθ / sθ * S2 - dV * S - V * S1
    sg = conjugate ? -1 : 1
    q = sg * (m / sθ - a * ω * sθ)
    dq = sg * (-m * cθ / sθ^2 - a * ω * cθ)
    d2q = sg * (m * (sθ^2 + 2cθ^2) / sθ^3 + a * ω * sθ)
    L = S1 + q * S; L1 = S2 + dq * S + q * S1; L2 = S3 + d2q * S + 2dq * S1 + q * S2
    vals = (S, S1, S2, L, L1, L2)
    return conjugate ? conj.(vals) : vals
end

"""
    gw_kerr_source_hat(ma, mb, kind, rn, θn) -> (Shat, Ph, ω̃, m̃)

Direct (operator) form of the Teukolsky source of the pair on the (rn, θn) grid, for
teuk_Z_hat (Z^∞ with R^in, the horizon amplitude Z^H with R^up):

    Shat[i, j] = ρ^{-5} ρ̄^{-1} (B'_2 + B'*_2)/Ph_i,   Ph = e^{-iω̃ r* + i m̃ r#}   (ST03 eqs. 9-12).

The tetrad components are those of gw_kerr_tetrad, written as Ph × (horizon-regular
part): n·∂ψ_j = -Δ f_j' ph_j S_j/(2Σ), m̄·∂ψ_j = (ρ/√2) L_j f_j ph_j, with
ph_j = e^{-i Re ω_j r* + i m_j r#} (quasi-stationary source at Re ω), so that ph_a ph_b
(ph_a ph_b* for transitions) = Ph exactly. Their second-order (r, θ) jets come from
f, f', f'' (Leaver series), f''' (ForwardDiff) and the spheroidal equation, and the
source operators act on them with J_+ -> ∂_r (teuk_source_B2_jet phase_factored). No
integration by parts, so nothing is dropped at r_+.

`window` (tests): a function r -> (W, W', W'') multiplying all three components.
"""
function gw_kerr_source_hat(ma::KerrMode, mb::KerrMode, kind::Symbol, rn, θn; window=nothing)
    kind in (:ann, :tr) || error("gw_kerr_source_hat: kind must be :ann or :tr")
    a = ma.a
    a == mb.a && ma.alph == mb.alph || error("gw_kerr_source_hat: the two modes must share (α, a)")
    ωa, ωb = real(ma.ω), real(mb.ω)
    same = (ma.n, ma.l, ma.m) == (mb.n, mb.l, mb.m)
    conjb = kind == :tr
    conjb && same && error("gw_kerr_source_hat: a transition needs two distinct levels")
    ω̃ = conjb ? ωa - ωb : ωa + ωb
    m̃ = conjb ? ma.m - mb.m : ma.m + mb.m
    c = same ? 1 / (2ωa) : 1 / sqrt(ωa * ωb)
    TR = promote_type(eltype(rn), Float64)
    # radial jets: G_j = -Δ f_j' (n·∂ψ_j factor), F_j = f_j (m̄·∂ψ_j factor)
    radj = Vector{NTuple{4,Jet2{TR}}}(undef, length(rn))
    Ph = Vector{Complex{TR}}(undef, length(rn))
    for (i, r) in enumerate(rn)
        Δ = r^2 - 2r + a^2; dΔ = 2(r - 1)
        js = map((ma, mb)) do md
            f, f1, f2 = kerr_radial_f(md, r)
            f3 = kerr_radial_f3(md, r)
            G = -(Δ * f1); G1 = -(dΔ * f1 + Δ * f2); G2 = -(2f1 + 2dΔ * f2 + Δ * f3)
            (G, G1, G2, f, f1, f2)
        end
        (Ga, G1a, G2a, Fa, F1a, F2a), (Gb, G1b, G2b, Fb, F1b, F2b) = js
        if conjb
            Gb, G1b, G2b, Fb, F1b, F2b = conj.((Gb, G1b, G2b, Fb, F1b, F2b))
        end
        radj[i] = (jet_radial(Ga, G1a, G2a, TR), jet_radial(Fa, F1a, F2a, TR), jet_radial(Gb, G1b, G2b, TR),
                   jet_radial(Fb, F1b, F2b, TR))
        Ph[i] = exp(-im * ω̃ * kerr_rstar(r, a) + im * m̃ * kerr_rsharp(r, a))
    end
    angj = map(θn) do θ
        Sa, S1a, S2a, La, L1a, L2a = kerr_angular_jet(ma, θ)
        Sb, S1b, S2b, Lb, L1b, L2b = kerr_angular_jet(mb, θ; conjugate=conjb)
        (jet_angular(Sa, S1a, S2a, TR), jet_angular(La, L1a, L2a, TR), jet_angular(Sb, S1b, S2b, TR), jet_angular(Lb, L1b, L2b, TR))
    end
    Shat = zeros(Complex{TR}, length(rn), length(θn))
    is2 = 1 / sqrt(TR(2))
    for (i, r) in enumerate(rn)
        jGa, jFa, jGb, jFb = radj[i]
        jr = jet_var_r(TR(r))
        if window !== nothing
            jGa = jGa * jet_radial(window(r)..., TR)        # n_a and m_a both carry one factor of W
            jFa = jFa * jet_radial(window(r)..., TR)
        end
        for (j, θ) in enumerate(θn)
            jSa, jLa, jSb, jLb = angj[j]
            jθ = jet_var_t(TR(θ))
            cθ = jet_cos(jθ)
            inv2Σ = inv(2 * (jr * jr + TR(a)^2 * (cθ * cθ)))
            ρ = inv(jr + (-im * TR(a)) * cθ)
            na = jGa * jSa * inv2Σ; nb = jGb * jSb * inv2Σ
            mA = is2 * (ρ * jLa * jFa); mB = is2 * (ρ * jLb * jFb)
            Tnn = c * (na * nb)
            Tnm = (c / 2) * (na * mB + nb * mA)
            Tmm = c * (mA * mB)
            ρv = ρ.v
            Shat[i, j] = ρv^-5 * conj(ρv)^-1 * teuk_source_B2_jet(Tnn, Tnm, Tmm, a, ω̃, m̃, r, θ; phase_factored=true)
        end
    end
    return Shat, Ph, ω̃, m̃
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
    # out to where the product has decayed: decay/qsum, or, for radial overtones (whose decay only starts at the
    # outer turning point ~2n²/α²), the end of the more compact mode's ODE range (km_ode_range)
    xmax = max(decay / qsum, min(km_ode_range(ma.alph, ma.a, ma.n, ma.L.q)[2], km_ode_range(mb.alph, mb.a, mb.n, mb.L.q)[2]))
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
                 T=Float64, verbose=false, horizon=false, teuk_kwargs...)
        -> (Γ, P_∞, Dict ℓ => P_∞ℓ, ω̃, m̃, Γ_H, Dict ℓ => P_Hℓ, chk)

Rate of the annihilation (kind = :ann) or transition a -> b (kind = :tr,
needs Re ω_a > Re ω_b) in the solver's convention, summing ℓ from max(2, |m̃|)
until two consecutive ℓ contribute less than ltol of the total. Γ and P_∞ are the
flux to infinity (kernel form, teuk_Z_grid).

With horizon = true the flux into the horizon is computed too, from the direct
(phase-factored) form of the source (gw_kerr_source_hat, teuk_Z_hat with R^up; the
kernel form cannot be used there: its integrations by parts drop boundary terms that
do not vanish at r_+ for R^up): P_H = Σ_ℓ teuk_horizon_alpha |Z^H_ℓ|²/(2πω̃²) (both
helicities) and Γ_H = P_H/(α ω̃), negative for superradiant gravitons
(k = ω̃ - m̃ Ω_H < 0). The total rate of the process is Γ + Γ_H, and P_H enters the BH
mass and angular-momentum budget. The same direct form gives P_∞ again, and
chk = P_∞(direct)/P_∞(kernel) - 1 is returned as a consistency check of the source
jets and the quadrature. Γ_H and chk are NaN without horizon.
"""
function gw_kerr_rate(ma::KerrMode, mb::KerrMode, kind::Symbol; lmax::Int=12, ltol=1e-8, nθ::Int=48, eps=1e-6,
                      ppr::Real=8.0, ppw::Real=8.0, decay::Real=60.0, T=Float64, verbose::Bool=false,
                      horizon::Bool=false, teuk_kwargs...)
    ω̃0 = kind == :tr ? real(ma.ω) - real(mb.ω) : real(ma.ω) + real(mb.ω)
    ω̃0 > 0 || error("gw_kerr_rate: ω̃ = $ω̃0 ≤ 0 (order a transition as a -> b with Re ω_a > Re ω_b)")
    rn, rw, θn, θw = gw_kerr_grid(ma, mb, ω̃0; eps=eps, ppr=ppr, ppw=ppw, decay=decay, nθ=nθ)
    rn = T.(rn); rw = T.(rw)          # the r-dependent part of the integral is done in precision T
    Tnn, Tnm, Tmm, ω̃, m̃ = gw_kerr_tetrad(ma, mb, kind, rn, θn)
    horizon && ((Shat, Ph) = gw_kerr_source_hat(ma, mb, kind, rn, θn)[1:2])
    modes = Dict{Int,Float64}()
    modesH = Dict{Int,Float64}()
    P = 0.0; PH = 0.0; Pd = 0.0
    nsmall = 0
    # R^in must be available over the whole source grid and is stepped through its nodes; for the
    # horizon flux R^up is carried in to the first node as well
    kw = merge((r_in_max=1.05 * rn[end], tstops_in=rn, T=T), values(teuk_kwargs))
    horizon && (kw = merge(kw, (r_up_min=rn[1], tstops_up=rn)))
    for l in max(2, abs(m̃)):lmax
        tr = teuk_radial(-2, l, m̃, T(ω̃), T(ma.a); kw...)
        Z = teuk_Z_grid(tr, Tnn, Tnm, Tmm, rn, rw, θn, θw)
        modes[l] = Float64(abs2(Z) / (2π * ω̃^2))
        P += modes[l]
        if horizon
            ZH = teuk_Z_hat(tr, Shat, Ph, rn, rw, θn, θw; horizon=true)
            modesH[l] = Float64(teuk_horizon_alpha(tr) * abs2(ZH) / (2π * ω̃^2))
            PH += modesH[l]
            Pd += Float64(abs2(teuk_Z_hat(tr, Shat, Ph, rn, rw, θn, θw)) / (2π * ω̃^2))
        end
        verbose && println("  ℓ = $l: |Z|²/(2πω̃²) = $(modes[l])" * (horizon ? "   P_H = $(modesH[l])" : "") * "   W_spread = $(tr.W_spread)")
        # converged when two consecutive ℓ add less than ltol of the total (horizon part included: near-horizon
        # sources radiate into the hole at higher ℓ than to infinity)
        small = modes[l] < ltol * P && (!horizon || abs(modesH[l]) < ltol * (P + abs(PH)))
        nsmall = small ? nsmall + 1 : 0
        nsmall >= 2 && break
    end
    ΓH = horizon ? PH / (ma.alph * Float64(ω̃)) : NaN
    chk = horizon ? Pd / P - 1 : NaN
    return P / (ma.alph * Float64(ω̃)), P, modes, Float64(ω̃), m̃, ΓH, modesH, chk
end
