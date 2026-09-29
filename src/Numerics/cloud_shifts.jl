"""
    cloud_shifts.jl  (included by gw_rates.jl)

Frequency shifts of the levels from the cloud's own gravity and from the axion
quartic self-interaction, at first order, for an arbitrary set of occupied
hydrogenic levels (Schrödinger–Poisson / Gross–Pitaevskii with
ψ = Σ_j sqrt(N_j) ψ_j e^{-iω_j t}). Keeping every term resonant with level i:

  gravity:  δω_i/μ = -α^3 Σ_j u_j [ D_ij + (1 - δ_ij) E_ij ]
     D_ij = ∫∫ |ψ_i(x)|^2 |ψ_j(x')|^2 / |x - x'|          (direct)
     E_ij = ∫∫ ψ_i*ψ_j(x) ψ_j*ψ_i(x') / |x - x'|          (exchange)
  quartic:  δω_i/μ = -(α^5/8) (M_pl/f_a)^2 Σ_j u_j (2 - δ_ij) χ_ij,
     χ_ij = ∫ |ψ_i|^2 |ψ_j|^2 d^3x
     (from V ⊃ -μ^2 φ^4/(24 f_a^2); contact direct + exchange give the 2)

with u_j = N_j/(G M^2) and all integrals in Bohr units (a0 = 1/(μα)).
Baryakhtar et al. 2021 App. G keep only the direct gravity term; App. B (B14)
has the same quartic combination as here.
"""

# ∫_0^∞ r^p e^{-β r} dr
gw_gint(p, β) = factorial(big(p)) / β^(p + 1)

# I_L = ∫∫ f(r) g(r') r^2 r'^2 r_<^L / r_>^{L+1} dr dr'  for f, g = Σ c r^p e^{-β r}
function gw_coulomb_radial(f::GWRadial{T}, g::GWRadial{T}, L::Int) where {T}
    βf, βg = f.beta, g.beta
    tot = zero(T)
    for (pf, cf) in zip(f.pows, f.coefs), (pg, cg) in zip(g.pows, g.coefs)
        # r' < r: r^{-(L+1)} ∫_0^r x^{n1} e^{-βg x} = n1!/βg^{n1+1} [1 - e^{-βg r} Σ_k (βg r)^k/k!]
        n1 = pg + 2 + L
        p1 = pf + 1 - L
        s1 = gw_gint(p1, βf)
        for k in 0:n1
            s1 -= βg^k / factorial(big(k)) * gw_gint(p1 + k, βf + βg)
        end
        t1 = factorial(big(n1)) / βg^(n1 + 1) * s1
        # r' > r: r^L ∫_r^∞ x^{n2} e^{-βg x} = n2!/βg^{n2+1} e^{-βg r} Σ_k (βg r)^k/k!
        n2 = pg + 1 - L
        p2 = pf + 2 + L
        s2 = zero(T)
        for k in 0:n2
            s2 += βg^k / factorial(big(k)) * gw_gint(p2 + k, βf + βg)
        end
        t2 = factorial(big(n2)) / βg^(n2 + 1) * s2
        tot += cf * cg * (t1 + t2)
    end
    return tot
end

function gw_angular3(l1, m1, l2, m2, L, M)
    # ∫ Y*_{l1 m1} Y_{l2 m2} Y_{L M} dΩ  (nonzero for M = m1 - m2)
    M == m1 - m2 || return 0.0
    xs, ws = gauss((l1 + l2 + L) ÷ 2 + 3)
    return 2π * sum(w * gw_theta_lm(l1, m1, acos(x)) * gw_theta_lm(l2, m2, acos(x)) *
                    gw_theta_lm(L, M, acos(x)) for (x, w) in zip(xs, ws))
end

"""
    gw_shift_kernels(si, sj; radcache=Dict()) -> (KG, χ)

KG = D_ij + (i != j) E_ij and χ = χ_ij (Bohr units) for levels si, sj = (n,l,m).
"""
function gw_shift_kernels(si::NTuple{3,Int}, sj::NTuple{3,Int}; radcache=Dict{Any,BigFloat}())
    (ni, li, mi), (nj, lj, mj) = si, sj
    Ri = gw_hydrogen_radial(ni, li; T=BigFloat); Rj = gw_hydrogen_radial(nj, lj; T=BigFloat)
    rad(key, f, g, L) = get!(() -> gw_coulomb_radial(f, g, L), radcache, (key, L))
    fii = gw_radial_product(Ri, Ri, 0); fjj = gw_radial_product(Rj, Rj, 0)
    D = 0.0
    for L in 0:2:(2 * min(li, lj))
        ai = gw_angular3(li, mi, li, mi, L, 0)
        aj = gw_angular3(lj, mj, lj, mj, L, 0)
        (ai == 0 || aj == 0) && continue
        D += 4π / (2L + 1) * ai * aj * Float64(rad((:D, ni, li, nj, lj), fii, fjj, L))
    end
    E = 0.0
    if si != sj
        fij = gw_radial_product(Ri, Rj, 0)
        for L in max(abs(li - lj), abs(mi - mj)):(li + lj)
            isodd(li + lj + L) && continue
            h = gw_angular3(li, mi, lj, mj, L, mi - mj)
            h == 0 && continue
            E += 4π / (2L + 1) * h^2 * Float64(rad((:E, ni, li, nj, lj), fij, fij, L))
        end
    end
    # contact overlap
    xs, ws = gauss(li + lj + 3)
    angc = 2π * sum(w * gw_theta_lm(li, mi, acos(x))^2 * gw_theta_lm(lj, mj, acos(x))^2 for (x, w) in zip(xs, ws))
    χ = Float64(gw_radial_moment(gw_radial_product(fii, fjj, 2))) * angc
    return D + E, χ, D, E
end

"""
    gw_level_shift_matrices(I, J) -> (KG, χ)

Kernel matrices for target levels I (whose frequencies are wanted) and source
levels J (occupied), each a vector of (n,l,m).
"""
function gw_level_shift_matrices(I, J)
    KG = zeros(length(I), length(J)); X = zeros(length(I), length(J))
    radcache = Dict{Any,BigFloat}()
    for (a, si) in enumerate(I), (b, sj) in enumerate(J)
        KG[a, b], X[a, b] = gw_shift_kernels(si, sj; radcache=radcache)[1:2]
        si == sj && (X[a, b] *= 0.5)       # (2 - δ_ij) / 2 folded in below
    end
    return KG, X
end

"""
    gw_level_shifts(KG, X, uJ, α; fa=nothing) -> δω/μ for the target levels

uJ: occupations of the source levels at one time. fa (GeV) = nothing switches the
self-interaction part off.
"""
function gw_level_shifts(KG, X, uJ, α; fa=nothing)
    δ = -α^3 .* (KG * uJ)
    if fa !== nothing
        δ .-= (α^5 / 4) * (M_pl / fa)^2 .* (X * uJ)
    end
    return δ
end
