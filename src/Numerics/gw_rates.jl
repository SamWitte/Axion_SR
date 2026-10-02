"""
    gw_rates.jl

Gravitational-wave annihilation and transition rates for arbitrary pairs of
hydrogenic scalar levels |nlm>, in the non-relativistic (NR) limit.

Conventions (solve_system_unified.jl)
-------------------------------------
Occupations are u_i = N_i / (G M^2) and time is tau = mu t, so that every rate
here is a dimensionless number multiplying u_a u_b:

  annihilation a x b -> g : du_a/dtau = du_b/dtau = -rate u_a u_b
                            (a == b: du_a/dtau = -2 rate u_a^2)
  transition   a -> b + g : du_a/dtau = -rate u_a u_b,  du_b/dtau = +rate u_a u_b
                            (omega_a > omega_b; stimulated emission)

The solver builds a GWCache once (gw_build_cache) and calls gw_rhs! on every
right-hand-side evaluation, so the rates follow the evolving BH mass (α) and
spin (hyperfine part of Δω).

Physics (NR limit)
------------------
* Annihilation (omega_GW ~ 2 mu, graviton wavelength << cloud): flat-space
  Weinberg formula with the TT part of the kinetic stress T_ij = d_i phi d_j phi
  (Arvanitaki, Baryakhtar & Huang 2015, App. A). The Fourier transform is done
  exactly: gradient formula -> Gaunt coefficients -> closed-form radial
  integrals of r^p e^{-beta r} j_J(kr). We return the leading term C alpha^p of
  the small-alpha expansion. Caveat: for 2p x 2p the flat-space leading term
  cancels for every l=1 pair (alpha^16 instead of the alpha^14 found once the BH potential is
  included); see `gw_literature_annihilation`.
* Transition (omega_GW = Delta omega ~ mu alpha^2, wavelength >> cloud): leading
  mass (I_L) and current (J_L) multipoles of the cross density/current, summed
  over every allowed L with the Thorne/Blanchet flux
      P = G sum_L [ (L+1)(L+2)/((L-1)L L!(2L+1)!!) |I_L^(L+1)|^2
                  + 4L(L+2)/((L-1)(L+1)!(2L+1)!!) |J_L^(L+1)|^2 ].
  Coefficients are alpha-independent; Delta omega enters at load time.

All lengths are in Bohr units a0 = 1/(mu alpha).

Hook for relativistic rates: `gw_rates_relativistic` (not implemented yet).
"""

using WignerSymbols: clebschgordan
using SpecialFunctions: loggamma
using QuadGK: gauss
using ForwardDiff

# ----------------------------------------------------------------------------
# Spherical harmonics (orthonormal, Condon-Shortley phase)
# ----------------------------------------------------------------------------

"""
    gw_theta_lm(l, m, θ)

theta part of Y_lm, so that Y_lm(θ,φ) = gw_theta_lm(l,m,θ) e^{imφ}. Works for
negative m and for ForwardDiff duals.
"""
function gw_theta_lm(l::Int, m::Int, θ)
    if abs(m) > l
        return zero(cos(θ))
    end
    if m < 0
        return (isodd(m) ? -1 : 1) * gw_theta_lm(l, -m, θ)
    end
    x = cos(θ)
    s = sin(θ)
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
        p_next = a * (x * p_cur - b * p_prev)
        p_prev, p_cur = p_cur, p_next
    end
    return p_cur
end

gw_dtheta_lm(l, m, θ) = ForwardDiff.derivative(t -> gw_theta_lm(l, m, t), θ)

# ----------------------------------------------------------------------------
# Hydrogenic radial functions as sum_k c_k r^{p_k} e^{-beta r}
# ----------------------------------------------------------------------------

struct GWRadial{T}
    pows::Vector{Int}
    coefs::Vector{T}
    beta::T
end

"""
    gw_hydrogen_radial(n, l; T=Float64)

R_nl in Bohr units, normalised to ∫ R^2 r^2 dr = 1, positive at small r.
"""
function gw_hydrogen_radial(n::Int, l::Int; T=Float64)
    N = n - l - 1
    norm = sqrt(T(2) / n)^3 * sqrt(T(factorial(big(N))) / (2n * T(factorial(big(n + l)))))
    pows = Int[]
    coefs = T[]
    for j in 0:N
        c = norm * (T(2) / n)^(l + j) * (isodd(j) ? -1 : 1) *
            T(binomial(big(n + l), N - j)) / T(factorial(big(j)))
        push!(pows, l + j)
        push!(coefs, c)
    end
    return GWRadial{T}(pows, coefs, T(1) / n)
end

function gw_simplify(pows, coefs, beta::T) where {T}
    d = Dict{Int, T}()
    for (p, c) in zip(pows, coefs)
        d[p] = get(d, p, zero(T)) + c
    end
    ks = sort(collect(keys(d)))
    return GWRadial{T}(ks, [d[k] for k in ks], beta)
end

# g_up = R' - l R / r  (couples to L = l+1),  g_dn = R' + (l+1) R / r  (L = l-1)
function gw_grad_radial(R::GWRadial{T}, l::Int, up::Bool) where {T}
    pows = Int[]
    coefs = T[]
    shift = up ? -l : (l + 1)
    for (p, c) in zip(R.pows, R.coefs)
        push!(pows, p - 1); push!(coefs, c * (p + shift))
        push!(pows, p);     push!(coefs, -R.beta * c)
    end
    keep = coefs .!= 0
    return gw_simplify(pows[keep], coefs[keep], R.beta)
end

function gw_radial_product(A::GWRadial{T}, B::GWRadial{T}, extra_pow::Int) where {T}
    pows = Int[]
    coefs = T[]
    for (pa, ca) in zip(A.pows, A.coefs), (pb, cb) in zip(B.pows, B.coefs)
        push!(pows, pa + pb + extra_pow)
        push!(coefs, ca * cb)
    end
    return gw_simplify(pows, coefs, A.beta + B.beta)
end

# ∫_0^∞ R(r) dr for R = sum c r^p e^{-beta r}
function gw_radial_moment(R::GWRadial{T}) where {T}
    s = zero(T)
    for (p, c) in zip(R.pows, R.coefs)
        s += c * T(factorial(big(p))) / R.beta^(p + 1)
    end
    return s
end

# ----------------------------------------------------------------------------
# Closed-form  k^{p+1} ∫_0^∞ r^p e^{-beta r} j_J(k r) dr,  as a function of b = beta/k.
# Uses h_J^(1) = (-i)^{J+1} e^{iz}/z sum_s i^s (J+s)!/(s!(J-s)!(2z)^s). The single
# log-divergent s = J = p term is purely imaginary; its finite real part is
# (2J)!/(J! 2^J) arctan(1/b).
# ----------------------------------------------------------------------------

function gw_Itilde(p::Int, J::Int, b::Float64)
    p >= J || error("gw_Itilde requires p >= J (p=$p, J=$J)")
    λ = complex(b, -1.0)
    tot = zero(ComplexF64)
    for s in 0:J
        ν = p - s
        ν == 0 && continue
        la = loggamma(J + s + 1) - loggamma(s + 1) - loggamma(J - s + 1) - s * log(2.0) + loggamma(ν)
        tot += im^s * exp(la) * λ^(-ν)
    end
    res = real((-im)^(J + 1) * tot)
    if p == J
        res += exp(loggamma(2J + 1) - loggamma(J + 1) - J * log(2.0)) * atan(1 / b)
    end
    return res
end

# ----------------------------------------------------------------------------
# Gradient of psi_nlm in the spherical basis:  ∇ψ = sum_q c_q e_q,
#   c_q = sum_{L=l±1} A_{L,q} g_L(r) Y_{L, m-q}
# e_{+1} = -(x+iy)/√2, e_0 = z, e_{-1} = (x-iy)/√2
# ----------------------------------------------------------------------------

# WignerSymbols keeps a global, non-thread-safe cache: memoise per thread and
# take a lock on a miss.
const GW_CG_LOCK = ReentrantLock()
const GW_CG_CACHES = [Dict{NTuple{6,Int}, Float64}() for _ in 1:Threads.nthreads()]
function gw_cg(j1, m1, j2, m2, j3, m3)
    c = GW_CG_CACHES[Threads.threadid()]
    get!(c, (j1, m1, j2, m2, j3, m3)) do
        lock(GW_CG_LOCK) do
            clebschgordan(Float64, j1, m1, j2, m2, j3, m3)
        end
    end
end

const GW_EQ = Dict(
    1  => ComplexF64[-1 / sqrt(2), -im / sqrt(2), 0],
    0  => ComplexF64[0, 0, 1],
    -1 => ComplexF64[1 / sqrt(2), -im / sqrt(2), 0],
)

struct GWGradTerm
    q::Int
    L::Int
    mu::Int          # azimuthal index of Y_{L,mu}, mu = m - q
    A::Float64
end

function gw_grad_terms(l::Int, m::Int)
    terms = GWGradTerm[]
    for q in -1:1
        μ = m - q
        L = l + 1
        if abs(μ) <= L
            A = -sqrt((l + 1) / (2l + 1)) * gw_cg(L, μ, 1, q, l, m)
            A != 0 && push!(terms, GWGradTerm(q, L, μ, A))
        end
        L = l - 1
        if L >= 0 && abs(μ) <= L
            A = sqrt(l / (2l + 1)) * gw_cg(L, μ, 1, q, l, m)
            A != 0 && push!(terms, GWGradTerm(q, L, μ, A))
        end
    end
    return terms
end

"""
    gw_grad_psi(n, l, m, x) -> Vector{ComplexF64}

Cartesian gradient of the normalised hydrogenic psi_nlm at Cartesian point x
(Bohr units), built from the spherical-basis gradient formula. Used in tests.
"""
function gw_grad_psi(n, l, m, x)
    r = sqrt(sum(abs2, x)); θ = acos(x[3] / r); φ = atan(x[2], x[1])
    R = gw_hydrogen_radial(n, l)
    out = zeros(ComplexF64, 3)
    for t in gw_grad_terms(l, m)
        g = gw_grad_radial(R, l, t.L == l + 1)
        gval = sum(c * r^p for (p, c) in zip(g.pows, g.coefs)) * exp(-g.beta * r)
        out .+= t.A * gval * gw_theta_lm(t.L, t.mu, θ) * exp(im * t.mu * φ) .* GW_EQ[t.q]
    end
    return out
end

function gw_psi(n, l, m, x)
    r = sqrt(sum(abs2, x)); θ = acos(x[3] / r); φ = atan(x[2], x[1])
    R = gw_hydrogen_radial(n, l)
    return sum(c * r^p for (p, c) in zip(R.pows, R.coefs)) * exp(-R.beta * r) *
           gw_theta_lm(l, m, θ) * exp(im * m * φ)
end

gw_gaunt(L1, m1, L2, m2, J, M) =
    sqrt((2L1 + 1) * (2L2 + 1) / (4π * (2J + 1))) *
    gw_cg(L1, 0, L2, 0, J, 0) * gw_cg(L1, m1, L2, m2, J, M)

# ----------------------------------------------------------------------------
# Annihilation: flat-space power
# ----------------------------------------------------------------------------

"""
    gw_annihilation_F(sa, sb, α) -> (Fhat, Pmin)

Angular-integrated TT power functional F = ∫dΩ Λ_{ij,kl} T*_ij T_kl for the
source S_ij = ∂_iψ_a ∂_jψ_b (+ i<->j if a != b), evaluated at κ = k a0 = 2/α.
Returned as F = κ^{-2 Pmin} Fhat to avoid underflow.
"""
function gw_annihilation_F(sa::NTuple{3,Int}, sb::NTuple{3,Int}, α::Float64)
    (na, la, ma), (nb, lb, mb) = sa, sb
    κ = 2 / α
    Pmin = la + lb - 1
    Ra = gw_hydrogen_radial(na, la); Rb = gw_hydrogen_radial(nb, lb)
    ta = gw_grad_terms(la, ma); tb = gw_grad_terms(lb, mb)
    ga = Dict(L => gw_grad_radial(Ra, la, L == la + 1) for L in (la - 1, la + 1) if L >= 0)
    gb = Dict(L => gw_grad_radial(Rb, lb, L == lb + 1) for L in (lb - 1, lb + 1) if L >= 0)

    nθ = la + lb + 6
    xs, ws = gauss(nθ)
    θs = acos.(xs)
    T = zeros(ComplexF64, 3, 3, nθ)

    radcache = Dict{Tuple{Int,Int,Int}, Float64}()
    function radJ(La, Lb, J)
        get!(radcache, (La, Lb, J)) do
            P = gw_radial_product(ga[La], gb[Lb], 2)
            b = P.beta / κ
            s = 0.0
            for (p, c) in zip(P.pows, P.coefs)
                s += c * gw_Itilde(p, J, b) * κ^(-(p + 1 - Pmin))
            end
            s
        end
    end

    Yk = Dict{Tuple{Int,Int}, Vector{Float64}}()
    for a in ta, b in tb
        M = a.mu + b.mu
        eij = GW_EQ[a.q] * transpose(GW_EQ[b.q])
        for J in abs(a.L - b.L):(a.L + b.L)
            (isodd(a.L + b.L + J) || abs(M) > J) && continue
            G = gw_gaunt(a.L, a.mu, b.L, b.mu, J, M)
            G == 0 && continue
            amp = a.A * b.A * G * 4π * (-im)^J * radJ(a.L, b.L, J)
            y = get!(Yk, (J, M)) do
                [gw_theta_lm(J, M, θ) for θ in θs]
            end
            for k in 1:nθ
                @views T[:, :, k] .+= (amp * y[k]) .* eij
            end
        end
    end

    Fhat = 0.0
    for k in 1:nθ
        Tk = T[:, :, k]
        if sa != sb
            Tk = Tk + transpose(Tk)
        end
        n̂ = [sin(θs[k]), 0.0, cos(θs[k])]
        Pp = [1.0 0 0; 0 1.0 0; 0 0 1.0] - n̂ * n̂'
        A = Pp * Tk * Pp
        tr = A[1, 1] + A[2, 2] + A[3, 3]
        val = real(sum(conj.(Tk) .* (A .- 0.5 .* tr .* Pp)))
        Fhat += ws[k] * val
    end
    return 2π * Fhat, Pmin
end

"""
    gw_annihilation_rate_flat(sa, sb, α)

Full flat-space annihilation rate (all orders in α within the flat-space,
ω_GW = 2μ approximation), in the solver's convention (see file header).
"""
function gw_annihilation_rate_flat(sa, sb, α)
    Fhat, Pmin = gw_annihilation_F(sa, sb, α)
    return α^6 * (α / 2)^(2Pmin) * Fhat / (2π)
end

"""
    gw_annihilation_leading(sa, sb; α1=1e-3) -> (C, p)

Leading NR term, rate ≈ C α^p, extracted from the exact flat-space expression at
α1 and 2α1 (Richardson-corrected for the O(α) term).
"""
function gw_annihilation_leading(sa, sb; α1=1e-3)
    α2 = 2α1
    F1, Pmin = gw_annihilation_F(sa, sb, α1)
    F2, _ = gw_annihilation_F(sa, sb, α2)
    # rate/(α^{6+2Pmin} 2^{-2Pmin}/(2π)) = Fhat, which scales as α^{p - 6 - 2Pmin}
    (F1 == 0 || F2 == 0) && return (0.0, 0)
    q = round(Int, log(abs(F2 / F1)) / log(2.0))
    p = 6 + 2Pmin + q
    c1 = F1 / α1^q
    c2 = F2 / α2^q
    Cf = 2c1 - c2
    return Cf * 2.0^(-2Pmin) / (2π), p
end

# ----------------------------------------------------------------------------
# Transitions: multipole moments (α-independent)
# ----------------------------------------------------------------------------

gw_dfact(n) = n <= 0 ? 1.0 : prod(Float64, n:-2:1)
gw_cmass(L) = 4π * (L + 1) * (L + 2) / ((L - 1) * L * gw_dfact(2L + 1)^2)
gw_ccurr(L) = 16π * L * (L + 2) / ((L - 1) * (L + 1) * gw_dfact(2L + 1)^2)

"""
    gw_transition_channels(sa, sb; radcache=nothing, angcache=nothing) -> Vector{(L, kind, K)}

Multipole channels for a <-> b. kind = :mass or :current. With δ = Δω/μ, the rate
of each channel is K α^{2-2L} δ^{2L+1} (mass) or K α^{4-2L} δ^{2L+1} (current).
Optional Dict caches reuse the radial (m-independent) and angular (n-independent)
factors across pairs.
"""
function gw_transition_channels(sa::NTuple{3,Int}, sb::NTuple{3,Int}; radcache=nothing, angcache=nothing)
    (na, la, ma), (nb, lb, mb) = sa, sb
    m = ma - mb
    out = Tuple{Int, Symbol, Float64}[]
    Ls = max(2, abs(m), abs(la - lb)):(la + lb)
    isempty(Ls) && return out
    radf = function (extra)
        key = (na, la, nb, lb, extra)
        f() = Float64(gw_radial_moment(gw_radial_product(gw_hydrogen_radial(na, la; T=BigFloat),
                                                         gw_hydrogen_radial(nb, lb; T=BigFloat), extra)))
        radcache === nothing ? f() : get!(f, radcache, key)
    end
    angf = function (L)
        key = (la, ma, lb, mb, L)
        angcache === nothing ? gw_transition_angular(la, ma, lb, mb, L) :
                               get!(() -> gw_transition_angular(la, ma, lb, mb, L), angcache, key)
    end
    for L in Ls
        ang = angf(L)
        abs(ang) < 1e-14 && continue
        if iseven(la + lb + L)
            K = 2 * gw_cmass(L) * (radf(L + 2) * ang)^2
            K > 0 && push!(out, (L, :mass, K))
        else
            K = 2 * gw_ccurr(L) * ((π / L) * radf(L + 1) * ang)^2
            K > 0 && push!(out, (L, :current, K))
        end
    end
    return out
end

# Angular factor of the mass (even parity) or current (odd parity) moment.
function gw_transition_angular(la, ma, lb, mb, L)
    m = ma - mb
    if iseven(la + lb + L)
        xs, ws = gauss((la + lb + L) ÷ 2 + 2)
        return 2π * sum(w * gw_theta_lm(la, ma, acos(x)) * gw_theta_lm(lb, mb, acos(x)) *
                        gw_theta_lm(L, m, acos(x)) for (x, w) in zip(xs, ws))
    else
        θg, wθ = gauss(la + lb + L + 30, 0, π)
        return sum(w * (m * (gw_theta_lm(la, ma, θ) * gw_dtheta_lm(lb, mb, θ) -
                             gw_theta_lm(lb, mb, θ) * gw_dtheta_lm(la, ma, θ)) * gw_theta_lm(L, m, θ) -
                        (ma + mb) * gw_theta_lm(la, ma, θ) * gw_theta_lm(lb, mb, θ) * gw_dtheta_lm(L, m, θ))
                   for (θ, w) in zip(θg, wθ))
    end
end

gw_transition_rate(channels, α, δ) =
    sum((kind == :mass ? K * α^(2 - 2L) : K * α^(4 - 2L)) * δ^(2L + 1) for (L, kind, K) in channels; init=0.0)

# ----------------------------------------------------------------------------
# Level energies (same spectrum as ergL in solve_sr_rates.jl), in units of mu
# ----------------------------------------------------------------------------

gw_omega(n, l, m, α, a) = 1 - α^2 / (2n^2) - α^4 / (8n^4) + α^4 / n^4 * (2l - 3n + 1) / (l + 0.5) +
                          2a * m * α^5 / n^3 / (l * (l + 0.5) * (l + 1))

# ----------------------------------------------------------------------------
# Tables for a whole state list
# ----------------------------------------------------------------------------

"""
    gw_state_list(Nmax) -> Vector{NTuple{3,Int}}

The solver's level list: (n,l,m) with 1 <= m <= l < n <= Nmax, then the
truncation modes (2l+1, 2l, 2m) that fall outside Nmax (same order as
setup_quantum_levels_standard).
"""
function gw_state_list(Nmax::Int)
    st = NTuple{3,Int}[]
    for n in 1:Nmax, l in 1:(n - 1), m in 1:l
        push!(st, (n, l, m))
    end
    seen = Set{NTuple{3,Int}}()
    for n in 1:Nmax, l in 1:(n - 1), m in 1:l
        t = (2l + 1, 2l, 2m)
        if t[1] > Nmax && !(t in seen)
            push!(st, t); push!(seen, t)
        end
    end
    return st
end

gw_table_path(Nmax) = joinpath(@__DIR__, "..", "rate_sve", "gw_nr_rates_Nmax_$(Nmax).txt")

# |Δω/μ| at α = 1, maximised over BH spin: an upper bound used only for pruning.
gw_delta_bound(sa, sb) = maximum(abs(gw_omega(sa..., 1.0, a) - gw_omega(sb..., 1.0, a)) for a in (0.0, 1.0))

"""
    gw_build_nr_table(Nmax; path=gw_table_path(Nmax), floor_at_alpha1=1e-60, rel_channel=1e-8)

Compute leading NR annihilation (C, p) and transition multipole coefficients for
every unordered pair of levels in `gw_state_list(Nmax)` and write

  ann  na la ma nb lb mb  C   p   0                       rate = C α^p
  tr   na la ma nb lb mb  K   L   kind(1=mass,2=current)  rate = K α^{2-2L | 4-2L} δ^{2L+1}

Rates only fall with decreasing α, so entries are pruned at α = 1: an annihilation
is dropped if C < floor_at_alpha1; a transition channel is dropped if its α = 1
rate (with the largest Δω the spectrum allows) is below floor_at_alpha1 or below
rel_channel times the pair's dominant channel. `threaded=true` spreads pairs over
Julia threads (can deadlock on Julia 1.8 with many threads; serial is the default).
"""
function gw_build_nr_table(Nmax::Int; path=gw_table_path(Nmax), floor_at_alpha1=1e-60,
                           rel_channel=1e-8, verbose=true, threaded=false)
    states = gw_state_list(Nmax)
    pairs = [(states[i], states[j]) for i in 1:length(states) for j in i:length(states)]
    rows = [Vector{Any}[] for _ in 1:length(pairs)]
    nt = Threads.nthreads()
    radcaches = [Dict{NTuple{5,Int}, Float64}() for _ in 1:nt]
    angcaches = [Dict{NTuple{5,Int}, Float64}() for _ in 1:nt]
    done = Threads.Atomic{Int}(0)
    t0 = time()
    verbose && println("GW NR table: $(length(states)) states, $(length(pairs)) pairs, $(nt) threads")
    # Compile every code path serially first: Julia 1.8 can deadlock when many
    # threads JIT-compile the same methods at once.
    gw_annihilation_leading((2, 1, 1), (3, 2, 1))
    gw_transition_channels((3, 2, 2), (2, 1, 1))
    gw_transition_channels((3, 2, 2), (3, 1, 1); radcache=radcaches[1], angcache=angcaches[1])
    gw_delta_bound((3, 2, 2), (2, 1, 1))
    verbose && flush(stdout)
    work = function (ip)
        tid = Threads.threadid()
        sa, sb = pairs[ip]
        C, p = gw_annihilation_leading(sa, sb)
        C >= floor_at_alpha1 && push!(rows[ip], Any["ann", sa..., sb..., C, p, 0])
        if sa != sb
            ch = gw_transition_channels(sa, sb; radcache=radcaches[tid], angcache=angcaches[tid])
            if !isempty(ch)
                δ1 = gw_delta_bound(sa, sb)
                r1 = [K * δ1^(2L + 1) for (L, kind, K) in ch]
                keep = max(floor_at_alpha1, rel_channel * maximum(r1))
                for ((L, kind, K), r) in zip(ch, r1)
                    r >= keep && push!(rows[ip], Any["tr", sa..., sb..., K, L, kind == :mass ? 1 : 2])
                end
            end
        end
        c = Threads.atomic_add!(done, 1) + 1
        if verbose && (c % 20000 == 0)
            println("  $(c)/$(length(pairs)) pairs  ($(round(time() - t0, digits=1)) s)")
            flush(stdout)
        end
    end
    if threaded
        Threads.@threads :static for ip in 1:length(pairs)
            work(ip)
        end
    else
        foreach(work, 1:length(pairs))
    end
    open(path, "w") do io
        println(io, "# GW NR rates, Nmax=$(Nmax), solver convention (Numerics/gw_rates.jl).")
        println(io, "# ann: rate = C alpha^p.  tr: rate = K alpha^(2-2L mass | 4-2L current) (dw/mu)^(2L+1)")
        println(io, "# type na la ma nb lb mb coeff power kind")
        for rs in rows, r in rs
            println(io, join(r, " "))
        end
    end
    verbose && println("wrote $(path) in $(round(time() - t0, digits=1)) s")
    return path
end

# ----------------------------------------------------------------------------
# Literature overrides and the relativistic hook
# ----------------------------------------------------------------------------

# |R_nl(r)/r^l| at r -> 0 (Bohr units), squared.
gw_origin_amp2(n, l) = (2 / n)^(2l + 3) * exp(loggamma(n + l + 1) - loggamma(n - l) - 2 * loggamma(2l + 2)) / (2n)

"""
    gw_literature_annihilation(sa, sb) -> (C, p) or nothing

Leading-order annihilation terms that supersede the flat-space value. For every
l = 1 pair (n p x n' p) the flat-space O(α^14) term cancels, leaving α^16; the BH
potential restores α^14 with dE/dt = (484+9π^2)/23040 (M_c/M)^2 α^14 for 2p x 2p
(Yoshino & Kodama 2014; Brito et al. 2017; Baryakhtar et al. 2021 Table IV quotes
1e-2 α^14). The large-k tail is set by the wavefunctions at the origin, so other
n, n' follow from |R_n1(0)|^2 |R_n'1(0)|^2, times 4 for distinct levels (the
same scaling the flat-space coefficients obey exactly).
"""
function gw_literature_annihilation(sa, sb)
    (sa[2] == 1 && sb[2] == 1) || return nothing
    C211 = (484 + 9π^2) / 46080
    w(n) = gw_origin_amp2(n, 1) / gw_origin_amp2(2, 1)
    return (C211 * w(sa[1]) * w(sb[1]) * (sa == sb ? 1 : 4), 14)
end

function gw_rates_relativistic(args...; kwargs...)
    # Should return a GWCache-compatible object for (Nmax, modes, mu, M, a).
    error("Relativistic GW rates (gw_model=:rel) are not implemented yet. " *
          "This is the hook for Teukolsky-based annihilation/transition rates.")
end

# ----------------------------------------------------------------------------
# Solver-facing: rate cache evaluated inside the ODE right-hand side
# ----------------------------------------------------------------------------

const GW_TABLE_CACHE = Dict{String, Any}()

function gw_load_table(path)
    get!(GW_TABLE_CACHE, path) do
        ann = Tuple{NTuple{3,Int}, NTuple{3,Int}, Float64, Int}[]
        trans = Dict{Tuple{NTuple{3,Int},NTuple{3,Int}}, Vector{Tuple{Int,Symbol,Float64}}}()
        for line in eachline(path)
            startswith(line, "#") && continue
            f = split(line)
            sa = (parse(Int, f[2]), parse(Int, f[3]), parse(Int, f[4]))
            sb = (parse(Int, f[5]), parse(Int, f[6]), parse(Int, f[7]))
            if f[1] == "ann"
                push!(ann, (sa, sb, parse(Float64, f[8]), parse(Int, f[9])))
            else
                push!(get!(trans, (sa, sb), Tuple{Int,Symbol,Float64}[]),
                      (parse(Int, f[9]), f[10] == "1" ? :mass : :current, parse(Float64, f[8])))
            end
        end
        (ann, trans)
    end
end

"""
GW channels resolved to solver level indices. Annihilation k: rate C α^p.
Transition t: rate sum_q K_q α^{apow_q} |δ|^{dpow_q} over channels
q in off[t]:off[t+1]-1 (stored as ln K), with δ = ω_i - ω_j recomputed from the current α and
spin (the emission direction follows the sign of δ).
"""
struct GWCache
    n::Vector{Int}; l::Vector{Int}; m::Vector{Int}; ω::Vector{Float64}
    ann_i::Vector{Int}; ann_j::Vector{Int}; ann_C::Vector{Float64}; ann_p::Vector{Int}
    tr_i::Vector{Int}; tr_j::Vector{Int}; tr_off::Vector{Int}
    ch_lnK::Vector{Float64}; ch_apow::Vector{Int}; ch_dpow::Vector{Int}
end

gw_empty_cache(modes) = GWCache([md[1] for md in modes], [md[2] for md in modes], [md[3] for md in modes],
                                zeros(length(modes)), Int[], Int[], Float64[], Int[], Int[], Int[], [1],
                                Float64[], Int[], Int[])

"""
    gw_build_cache(Nmax, modes, mu, M, a; gw_model=:nonrel, min_rate_per_yr=1e-10,
                   literature_overrides=true) -> GWCache

Resolve the GW table for `Nmax` onto the solver's level list `modes`
((n,l,m,...) tuples, as from setup_quantum_levels_standard).

gw_model: :nonrel (NR tables, Numerics/gw_rates.jl), :rel (relativistic hook,
not implemented) or :off (no GW emission).

Channels with rate * mu/hbar * yr < min_rate_per_yr at the initial (α, a) are
dropped: with the default 1e-10 their e-fold time exceeds a Hubble time even at
u = 1 (α only decreases during spin-down, so this stays conservative). At
α ≳ 0.2 and Nmax ≳ 15 this still leaves O(10^5) mostly high-n transitions.
`literature_overrides` replaces the flat-space l=1 x l=1 annihilations by the
BH-potential-corrected α^14 result (gw_literature_annihilation).
With accretion the BH can grow: pass the largest mass reached as `M_cut`, and
a channel is kept if it passes the cut at either M or M_cut.
"""
function gw_build_cache(Nmax, modes, mu, M, a; gw_model=:nonrel, min_rate_per_yr=1e-10,
                        literature_overrides=true, M_cut=M)
    c = gw_empty_cache(modes)
    gw_model == :off && return c
    gw_model == :rel && return gw_rates_relativistic(Nmax, modes, mu, M, a)
    gw_model == :nonrel || error("unknown gw_model $(gw_model); use :nonrel, :rel or :off")
    path = gw_table_path(Nmax)
    isfile(path) || error("GW table $(path) missing; build it with `julia src/scripts/build_gw_tables.jl $(Nmax)`")
    ann, trans = gw_load_table(path)
    idx = Dict((md[1], md[2], md[3]) => i for (i, md) in enumerate(modes))
    α = GNew * M * mu
    αc = GNew * max(M, M_cut) * mu
    floor_rate = min_rate_per_yr / (mu / hbar * YEAR_IN_SECONDS)
    for (sa, sb, C, p) in ann
        (haskey(idx, sa) && haskey(idx, sb)) || continue
        if literature_overrides
            lit = gw_literature_annihilation(sa, sb)
            lit === nothing || ((C, p) = lit)
        end
        C * max(α^p, αc^p) >= floor_rate || continue
        push!(c.ann_i, idx[sa]); push!(c.ann_j, idx[sb]); push!(c.ann_C, C); push!(c.ann_p, p)
    end
    for ((sa, sb), ch) in trans
        (haskey(idx, sa) && haskey(idx, sb)) || continue
        δ = gw_omega(sa..., α, a) - gw_omega(sb..., α, a)
        δc = gw_omega(sa..., αc, a) - gw_omega(sb..., αc, a)
        ((δ != 0 && gw_transition_rate(ch, α, abs(δ)) >= floor_rate) ||
         (δc != 0 && gw_transition_rate(ch, αc, abs(δc)) >= floor_rate)) || continue
        push!(c.tr_i, idx[sa]); push!(c.tr_j, idx[sb])
        for (L, kind, K) in ch
            push!(c.ch_lnK, log(K)); push!(c.ch_apow, kind == :mass ? 2 - 2L : 4 - 2L); push!(c.ch_dpow, 2L + 1)
        end
        push!(c.tr_off, length(c.ch_lnK) + 1)
    end
    return c
end

@inline function gw_transition_rate_cached(c::GWCache, t, α, ad)
    lα, lδ = log(α), log(ad)
    r = 0.0
    @inbounds for q in c.tr_off[t]:(c.tr_off[t + 1] - 1)
        r += exp(c.ch_lnK[q] + c.ch_apow[q] * lα + c.ch_dpow[q] * lδ)
    end
    return r
end

"""
    gw_rhs!(du, u, c::GWCache, α, a)

Add GW annihilation and transition terms to du (solver units: rates in units of
mu multiplying u_i u_j, before the solver's mu/hbar*yr conversion), using the
current α = G M mu and BH spin a. Mass and spin are untouched: the gravitons
leave the system.
"""
function gw_rhs!(du, u, c::GWCache, α, a)
    @inbounds for i in eachindex(c.n)
        c.ω[i] = gw_omega(c.n[i], c.l[i], c.m[i], α, a)
    end
    @inbounds for k in eachindex(c.ann_i)
        i, j = c.ann_i[k], c.ann_j[k]
        r = c.ann_C[k] * α^c.ann_p[k] * u[i] * u[j]
        du[i] -= r
        du[j] -= r        # i == j: -2r, two quanta per annihilation
    end
    @inbounds for t in eachindex(c.tr_i)
        i, j = c.tr_i[t], c.tr_j[t]
        δ = c.ω[i] - c.ω[j]
        δ == 0 && continue
        r = gw_transition_rate_cached(c, t, α, abs(δ)) * u[i] * u[j]
        hi, lo = δ > 0 ? (i, j) : (j, i)
        du[hi] -= r
        du[lo] += r
    end
    return du
end

"""
    gw_channel_rates(c::GWCache, α, a) -> (ann, tr)

Instantaneous rates for inspection / signal modelling. ann: (i, j, rate) with
ω_GW = ω_i + ω_j; tr: (hi, lo, rate) for hi -> lo with ω_GW = ω_hi - ω_lo
(frequencies in units of mu, from gw_omega).
"""
function gw_channel_rates(c::GWCache, α, a)
    ann = [(c.ann_i[k], c.ann_j[k], c.ann_C[k] * α^c.ann_p[k]) for k in eachindex(c.ann_i)]
    tr = Tuple{Int,Int,Float64}[]
    for t in eachindex(c.tr_i)
        i, j = c.tr_i[t], c.tr_j[t]
        δ = gw_omega(c.n[i], c.l[i], c.m[i], α, a) - gw_omega(c.n[j], c.l[j], c.m[j], α, a)
        δ == 0 && continue
        hi, lo = δ > 0 ? (i, j) : (j, i)
        push!(tr, (hi, lo, gw_transition_rate_cached(c, t, α, abs(δ))))
    end
    return ann, tr
end

# ----------------------------------------------------------------------------
# Signal post-processing: frequency, power and strain of every GW line
# ----------------------------------------------------------------------------

const GW_ERG_PER_EV = 1.602176634e-12
const GW_G_CGS = 6.674e-8
const GW_C_CGS = 2.99792458e10
const GW_KPC_CM = 3.0857e21

"""
    gw_lines(timeT, states, modes, spin, mass, mu, Nmax; d_kpc=1.0, rel_floor=1e-8,
             n_peak_samples=400, gw_model=:nonrel, min_rate_per_yr=1e-10,
             literature_overrides=true)

GW lines of an evolved cloud. `states[i, k]` is u_i at time `timeT[k]` (years) for
level `modes[i]` = (n, l, m, ...) (rows of a Modes_ file; pruned lists are fine),
`spin`/`mass` are the BH spin and mass (M_sun) at the same times, `mu` in eV.

Each annihilation a x b and transition a -> b is one line with, at every time,
  f  [Hz]    = ω_GW μ / (2π ħ),   ω_GW = ω_a + ω_b  or  |ω_a - ω_b|   (gw_omega)
  P  [erg/s] = (emission events per s) × ħ ω_GW,
               events/s = rate(α(t), a(t)) u_a u_b G M_0^2 μ/ħ
  h0         = sqrt(8 G P / (c^3 (2πf)^2 d^2)),  i.e. sqrt(A_+^2 + A_x^2) averaged
               over source orientation, for a source at distance d_kpc.
Lines whose peak power (on `n_peak_samples` sampled times) is below
rel_floor × the brightest line are dropped. Output is on `n_out` log-spaced
times. Returns (t, lines): t [yr] and a vector of NamedTuples (kind, a, b, f,
df_cloud, P, h0) sorted by peak power.

Frequencies include the BH drift (α(t) = G M(t) μ, spin in the spectrum) and,
with `self_gravity=true` and/or `fa` (GeV) given, the cloud's self-gravity and
quartic self-interaction level shifts (cloud_shifts.jl) from every level
occupied above occ_rel × the largest occupation; `df_cloud` is that cloud part
of f. Powers use the unshifted spectrum, as in the evolution. `M0` is the BH mass
that normalises u (the initial mass of the run); pass it when `timeT` etc. are a
late-time slice of the run.
"""
function gw_lines(timeT, states, modes, spin, mass, mu, Nmax; d_kpc=1.0, rel_floor=1e-8,
                  n_peak_samples=400, n_out=2000, gw_model=:nonrel, min_rate_per_yr=1e-10,
                  literature_overrides=true, self_gravity=true, fa=nothing, occ_rel=1e-8,
                  M0=mass[1])
    nt = length(timeT)
    c = gw_build_cache(Nmax, modes, mu, mass[1], spin[1]; gw_model=gw_model,
                       min_rate_per_yr=min_rate_per_yr, literature_overrides=literature_overrides,
                       M_cut=maximum(mass))
    nlines = length(c.ann_i) + length(c.tr_i)
    ω = zeros(length(c.n))
    peak = zeros(nlines)
    for k in gw_time_samples(timeT, n_peak_samples)
        α = gw_set_omegas!(ω, c, mass[k], spin[k], mu)
        for q in 1:nlines
            peak[q] = max(peak[q], gw_line_at(c, q, ω, α, states, k, M0, mu)[2])
        end
    end
    lines = NamedTuple[]
    (nlines == 0 || maximum(peak) == 0) && return (t = Float64[], lines = lines)
    keep = findall(peak .>= rel_floor * maximum(peak))
    sort!(keep, by=q -> -peak[q])
    # cloud self-gravity / self-interaction shifts: target levels = those in kept
    # lines, sources = levels ever occupied above occ_rel × the largest occupation
    na = length(c.ann_i)
    lvl(q) = q <= na ? (c.ann_i[q], c.ann_j[q]) : (c.tr_i[q - na], c.tr_j[q - na])
    Iset = unique(vcat([collect(lvl(q)) for q in keep]...))
    umax = vec(maximum(states, dims=2))
    Jset = findall(umax .>= occ_rel * maximum(umax))
    shifts = self_gravity || fa !== nothing
    if shifts
        st(i) = (c.n[i], c.l[i], c.m[i])
        KG, X = gw_level_shift_matrices(st.(Iset), st.(Jset))
        self_gravity || (KG .= 0)
    end
    ks = gw_time_samples(timeT, n_out)
    F = zeros(length(keep), length(ks)); Pw = zeros(length(keep), length(ks)); dF = zeros(length(keep), length(ks))
    ωs = similar(ω)
    for (col, k) in enumerate(ks)
        α = gw_set_omegas!(ω, c, mass[k], spin[k], mu)
        ωs .= ω
        if shifts
            ωs[Iset] .+= gw_level_shifts(KG, X, states[Jset, k], α; fa=fa)
        end
        for (row, q) in enumerate(keep)
            F0, Pw[row, col] = gw_line_at(c, q, ω, α, states, k, M0, mu)
            i, j = lvl(q)
            F[row, col] = (q <= na ? ωs[i] + ωs[j] : abs(ωs[i] - ωs[j])) * mu / hbar / (2π)
            dF[row, col] = F[row, col] - F0
        end
    end
    d = d_kpc * GW_KPC_CM
    for (row, q) in enumerate(keep)
        isann = q <= na
        ia, ib = isann ? (c.ann_i[q], c.ann_j[q]) : (c.tr_i[q - na], c.tr_j[q - na])
        if !isann && ω[ia] < ω[ib]      # label transitions hi -> lo (at the last output time)
            ia, ib = ib, ia
        end
        f = F[row, :]
        h0 = [f[k] > 0 ? sqrt(8GW_G_CGS * Pw[row, k] / GW_C_CGS^3) / (2π * f[k] * d) : 0.0 for k in eachindex(f)]
        push!(lines, (kind = isann ? :annihilation : :transition,
                      a = (c.n[ia], c.l[ia], c.m[ia]), b = (c.n[ib], c.l[ib], c.m[ib]),
                      f = f, df_cloud = dF[row, :], P = Pw[row, :], h0 = h0))
    end
    return (t = timeT[ks], lines = lines)
end

# Up to n indices into timeT, log-spaced in time (plus the first and last sample).
function gw_time_samples(timeT, n)
    length(timeT) <= n && return collect(eachindex(timeT))
    tpos = timeT[timeT .> 0]
    targets = exp10.(range(log10(minimum(tpos)), log10(timeT[end]), length=n))
    ks = [clamp(searchsortedfirst(timeT, t), 1, length(timeT)) for t in targets]
    return unique(vcat(1, ks, length(timeT)))
end

function gw_set_omegas!(ω, c::GWCache, M, a, mu)
    α = GNew * M * mu
    @inbounds for i in eachindex(c.n)
        ω[i] = gw_omega(c.n[i], c.l[i], c.m[i], α, a)
    end
    return α
end

# (f [Hz], P [erg/s]) of line q at output time index k; q <= n_ann: annihilation.
function gw_line_at(c::GWCache, q, ω, α, states, k, M0, mu)
    na = length(c.ann_i)
    if q <= na
        i, j = c.ann_i[q], c.ann_j[q]
        wgw = ω[i] + ω[j]
        r = c.ann_C[q] * α^c.ann_p[q]
    else
        t = q - na
        i, j = c.tr_i[t], c.tr_j[t]
        wgw = abs(ω[i] - ω[j])
        r = wgw == 0 ? 0.0 : gw_transition_rate_cached(c, t, α, wgw)
    end
    events_per_s = r * states[i, k] * states[j, k] * (GNew * M0^2 * M_to_eV) * mu / hbar
    return wgw * mu / hbar / (2π), events_per_s * wgw * mu * GW_ERG_PER_EV
end

include(joinpath(@__DIR__, "cloud_shifts.jl"))
