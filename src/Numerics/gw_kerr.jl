"""
    gw_kerr.jl

Relativistic GW annihilation (a × b) and transition (a -> b) rates of a scalar cloud on Kerr: Teukolsky amplitude
(teukolsky.jl) sourced by the kerr_modes.jl modes at Re ω; G = c = M = 1, μ = α, Γ = P/(α ω̃) in solver units.
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
    # n.dpsi = -Delta f' e^{-iw r* + im r#} S/(2 Sigma), mbar.dpsi = (rho/sqrt2)[S' + (m/sin - a w sin) S] R;
    # the source is taken at Re w (with complex w the horizon derivatives carry spurious e^{Im w r*}).
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
Direct (no integration by parts, valid at r_+) Teukolsky source for teuk_Z_hat: Shat = ρ^{-5} ρ̄^{-1} (B'_2 + B'*_2)/Ph,
Ph = e^{-iω̃ r* + i m̃ r#} (ST03 eqs. 9-12). `window` (tests): r -> (W, W', W'') multiplying the source.
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

"""Outer turning point of the hydrogenic orbit of a mode's level, n²/α² (1 + sqrt(1 - l(l+1)/n²)) (M = 1)."""
gw_kerr_turning(md::KerrMode) = md.n^2 / md.alph^2 * (1 + sqrt(max(0.0, 1 - md.l * (md.l + 1) / md.n^2)))

"""
    gw_kerr_rmax(ma, mb, decay) -> x_max (= r_max - r_+)
Where R_a R_b has decayed by e^{-decay}; each decays as e^{-|q| r} only beyond its turning point (gw_kerr_turning).
"""
function gw_kerr_rmax(ma::KerrMode, mb::KerrMode, decay::Real)
    qa, qb = abs(real(ma.L.q)), abs(real(mb.L.q))
    ra, rb = gw_kerr_turning(ma), gw_kerr_turning(mb)
    (r1, q1), (r2, q2) = ra <= rb ? ((ra, qa), (rb, qb)) : ((rb, qb), (ra, qa))
    D1 = q1 * (r2 - r1)
    x = D1 >= decay ? r1 + decay / q1 : r2 + (decay - D1) / (q1 + q2)
    return max(x, decay / (qa + qb))
end

"""
    gw_kerr_grid(ma, mb, ω̃; eps=1e-6, ppr=8, ppw=8, decay=60, nθ=48, per_panel=8)
Gauss-Legendre panels in r from r_+ + eps to gw_kerr_rmax(decay) (annihilations need decay ≥ 45), each of width
min((r - r_+)/ppr, (2π/ω̃)/ppw); Gauss-Legendre in θ.
"""
function gw_kerr_grid(ma::KerrMode, mb::KerrMode, ω̃; eps=1e-6, ppr::Real=8.0, ppw::Real=8.0, decay::Real=60.0,
                      nθ::Int=48, per_panel::Int=8)
    rp = 1 + sqrt(1 - ma.a^2)
    xmax = gw_kerr_rmax(ma, mb, decay)
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
    gw_kerr_rate(ma, mb, kind; lmax=12, ltol=1e-8, horizon=false, kw...) -> (Γ, P_∞, P_∞ℓ, ω̃, m̃, Γ_H, P_Hℓ, chk)
Annihilation (:ann) or transition a -> b (:tr, Re ω_a > Re ω_b) rate to infinity in solver units, ℓ summed to ltol;
horizon = true adds Γ_H (< 0 if superradiant) and chk = P_∞(direct)/P_∞(kernel) - 1 (else NaN).
"""
function gw_kerr_rate(ma::KerrMode, mb::KerrMode, kind::Symbol; lmax::Int=12, ltol=1e-8, nθ::Int=48, eps=1e-6,
                      ppr::Real=8.0, ppw::Real=8.0, decay::Real=60.0, T=Float64, verbose::Bool=false,
                      horizon::Bool=false, ltrunc=1e-6, teuk_kwargs...)
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
    # R^in (and R^up for the horizon flux) stepped through the source nodes
    kw = merge((r_in_max=1.05 * rn[end], tstops_in=rn, T=T), values(teuk_kwargs))
    horizon && (kw = merge(kw, (r_up_min=rn[1], tstops_up=rn)))
    l0 = max(2, abs(m̃))
    for l in l0:lmax
        tr = try
            teuk_radial(-2, l, m̃, T(ω̃), T(ma.a); kw...)
        catch
            # Teukolsky failure at high l ends the sum if the previous l added < ltrunc of the total
            if l > l0 && modes[l - 1] < ltrunc * P && (!horizon || abs(modesH[l - 1]) < ltrunc * (P + abs(PH)))
                verbose && println("  ℓ = $l: Teukolsky solution failed; sum ended at ℓ = $(l - 1), whose term is $(modes[l - 1] / P) of the total")
                break
            end
            rethrow()
        end
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
        # converged when two consecutive l add < ltol of the total (horizon part included)
        small = modes[l] < ltol * P && (!horizon || abs(modesH[l]) < ltol * (P + abs(PH)))
        nsmall = small ? nsmall + 1 : 0
        nsmall >= 2 && break
    end
    ΓH = horizon ? PH / (ma.alph * Float64(ω̃)) : NaN
    chk = horizon ? Pd / P - 1 : NaN
    return P / (ma.alph * Float64(ω̃)), P, modes, Float64(ω̃), m̃, ΓH, modesH, chk
end
