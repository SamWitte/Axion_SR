# Quasi-bound Kerr modes of a massive scalar in ingoing Kerr coordinates (G = M = 1, mu = alpha):
#     phi = f(r) S(theta) e^{-i w v + i m phi~},     f(r) = R_BL(r) e^{i w r* - i m r#}  (smooth at r+)
# R_BL from Leaver's series (Dolan 2007; the recurrence used by find_im_part in solve_sr_rates.jl):
#     R = (r-r+)^{-i sig} (r-r-)^{i sig + chi - 1} e^{q r} sum_n b_n u^n,   u = (r-r+)/(r-r-)
# The minimal solution b_n is obtained from the continued-fraction ratios (backward, stable), and the eigenvalue is
# refined with the same continued fraction, using the exact spheroidal eigenvalue (Y_lm-basis diagonalisation)
# instead of its small-c series.  S(theta) = sum_l' b_l' Y_l'm, so derivatives are exact.

# ---------------------------------------------------------------- spheroidal harmonics (s = 0)
"""Spheroidal eigenfunction of (L^2 - c2 cos^2) S = Lam S continuously connected to Y_lm.  Returns (Lam, coeffs::Dict l' => b)."""
function spheroidal_coeffs(l, m, c2; nextra=8, nth=60)
    ls = [lp for lp in abs(m):(l + 2nextra) if iseven(lp - l)]
    T = SWSHTable(m, maximum(ls), nth; smax=0)
    n = length(ls)
    A = zeros(ComplexF64, n, n)
    for (i, li) in enumerate(ls), (j, lj) in enumerate(ls)
        cc = 2pi * sum(T.w .* T.Y[(0, li)] .* T.Y[(0, lj)] .* T.x .^ 2)
        A[i, j] = (i == j ? li * (li + 1) : 0.0) - c2 * cc
    end
    E = eigen(A)
    i0 = findfirst(==(l), ls)
    k = argmax(abs.(E.vectors[i0, :]))
    v = E.vectors[:, k]
    v ./= sqrt(sum(v .^ 2))                    # ∑ b^2 = 1 (complex-symmetric normalisation)
    v .*= sign(real(v[i0]))
    return E.values[k], Dict(ls[i] => v[i] for i in 1:n)
end

"""S, dS/dθ, d²S/dθ² at the nodes of table T (which must have T.m == m and contain the needed l')."""
function spheroidal_on_nodes(coeffs, T::SWSHTable)
    S = zeros(ComplexF64, length(T.th)); dS = zeros(ComplexF64, length(T.th)); d2S = zeros(ComplexF64, length(T.th))
    for (lp, b) in coeffs
        lp > T.Lmax && (abs(b) > 1e-14 ? error("spheroidal expansion needs l' = $lp > table Lmax") : continue)
        S .+= b .* T.Y[(0, lp)]; dS .+= b .* T.dY[(0, lp)]; d2S .+= b .* T.d2Y[(0, lp)]
    end
    return S, dS, d2S
end

# ---------------------------------------------------------------- Leaver radial functions
struct LeaverCoeffs
    w::ComplexF64; alph::Float64; a::Float64; m::Int; Lam::ComplexF64
    q::ComplexF64; sig::ComplexF64; chi::ComplexF64
    al::Function; be::Function; ga::Function
end

function leaver_coeffs(w, alph, a, m, Lam)
    rp, rm = rplus(a), rminus(a); b = sqrt(1 - a^2)
    q = -sqrt(alph^2 - w^2 + 0im); real(q) > 0 && (q = -q)
    OmH = a / (2rp)
    sig = 2rp * (w - OmH * m) / (rp - rm)
    chi = (alph^2 - 2w^2) / q
    c0 = 1 - 2im * w - 2im / b * (w - a * m / 2)
    c1 = -4 + 4im * (w - im * q * (1 + b)) + 4im / b * (w - a * m / 2) - 2 * (w^2 + q^2) / q
    c2 = 3 - 2im * w - 2 * (q^2 - w^2) / q - 2im / b * (w - a * m / 2)
    c3 = 2im * (w - im * q)^3 / q + 2 * (w - im * q)^2 * b + q^2 * a^2 + 2im * q * a * m - Lam - 1 - (w - im * q)^2 / q +
         2q * b + 2im / b * ((w - im * q)^2 / q + 1) * (w - a * m / 2)
    c4 = (w - im * q)^4 / q^2 + 2im * w * (w - im * q)^2 / q - 2im / b * (w - im * q)^2 / q * (w - a * m / 2)
    al(n) = n^2 + (c0 + 1) * n + c0
    be(n) = -2n^2 + (c1 + 2) * n + c3
    ga(n) = n^2 + (c2 - 3) * n + c4
    return LeaverCoeffs(w, alph, a, m, Lam, q, sig, chi, al, be, ga)
end

"""Continued-fraction ratios t_n = -b_n/b_{n-1}, n = 1..K, and the eigenvalue residual β0/α0 - t_1."""
function leaver_cf(L::LeaverCoeffs, K)
    t = zeros(ComplexF64, K + 1)
    tmp = 0.0im
    for i in K:-1:1
        tmp = L.ga(i) / (L.be(i) - L.al(i) * tmp)
        t[i] = tmp
    end
    return t, L.be(0) / L.al(0) - t[1]
end

"""Number of series terms: b_n ~ exp(-2 sqrt(2 kappa (r+ - r-) n)); ask for ~1e-17 relative."""
cf_terms(alph, w, a) = (kap = real(sqrt(alph^2 - real(w)^2 + 0im)); Int(clamp(ceil((40 / (2 * sqrt(2kap * (rplus(a) - rminus(a)))))^2), 2000, 4_000_000)))

struct KerrMode
    n::Int; l::Int; m::Int
    alph::Float64; a::Float64
    w::ComplexF64
    Lam::ComplexF64
    sph::Dict{Int,ComplexF64}
    L::LeaverCoeffs
    b::Vector{ComplexF64}
    norm::Float64
end

"""Refine the eigenfrequency (secant on the continued fraction) starting from w0; build the mode."""
function kerr_mode(n, l, m, alph, a, w0; K=nothing, tol=1e-13, maxit=60, verbose=false)
    Kt = isnothing(K) ? cf_terms(alph, w0, a) : K
    function resid(w)
        Lam, _ = spheroidal_coeffs(l, m, a^2 * (w^2 - alph^2))
        _, r = leaver_cf(leaver_coeffs(w, alph, a, m, Lam), Kt)
        return r
    end
    w1 = complex(w0); w2 = w1 * (1 + 1e-7) + 1e-12im
    f1 = resid(w1); f2 = resid(w2)
    for it in 1:maxit
        w3 = w2 - f2 * (w2 - w1) / (f2 - f1)
        w1, f1 = w2, f2
        w2 = w3; f2 = resid(w2)
        verbose && println("  it $it  w = $w2  |res| = $(abs(f2))")
        abs(w2 - w1) < tol * abs(w2) && break
    end
    w = w2
    Lam, sph = spheroidal_coeffs(l, m, a^2 * (w^2 - alph^2))
    L = leaver_coeffs(w, alph, a, m, Lam)
    t, res = leaver_cf(L, Kt)
    b = zeros(ComplexF64, Kt + 1); b[1] = 1
    for i in 1:Kt
        b[i+1] = -t[i] * b[i]
        abs(b[i+1]) < 1e-300 && break
    end
    mode = KerrMode(n, l, m, alph, a, w, Lam, sph, L, b, 1.0)
    # normalise ∫|R|^2 r^2 dr = 1 (crudely on a log grid; only sets the scale)
    rr = 10 .^ range(log10(rplus(a) * (1 + 1e-6)), log10(60 * n^2 / alph^2), 20000)
    fv = [radial_f(mode, r)[1] for r in rr]
    nm = sum(abs2.(fv[1:end-1]) .* rr[1:end-1] .^ 2 .* diff(rr))
    return KerrMode(n, l, m, alph, a, w, Lam, sph, L, b, 1 / sqrt(nm)), res
end

"""f, f', f'' at (possibly complex) r for the smooth ingoing radial function."""
function radial_f(md::KerrMode, r)
    a = md.a; rp, rm = rplus(a), rminus(a); L = md.L; w = md.w; m = md.m
    p = im * L.sig + L.chi - 1 - im * (2rm * w - m * a) / (rp - rm)
    u = (r - rp) / (r - rm); du = (rp - rm) / (r - rm)^2; d2u = -2 * (rp - rm) / (r - rm)^3
    S0 = 0.0im; S1 = 0.0im; S2 = 0.0im
    # Horner for sum b_n u^n and derivatives
    for k in length(md.b):-1:1
        S2 = S2 * u + 2 * S1
        S1 = S1 * u + S0
        S0 = S0 * u + md.b[k]
    end
    lp = p / (r - rm) + (L.q + im * w)                  # P'/P
    P = exp(p * log(r - rm) + (L.q + im * w) * r)
    P1 = P * lp; P2 = P * (lp^2 - p / (r - rm)^2)
    f = P * S0
    f1 = P1 * S0 + P * S1 * du
    f2 = P2 * S0 + 2 * P1 * S1 * du + P * (S2 * du^2 + S1 * d2u)
    return md.norm * f, md.norm * f1, md.norm * f2
end
