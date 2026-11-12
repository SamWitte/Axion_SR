# Spin-weighted spherical harmonics sY_lm(theta) (Goldberg et al. 1967 convention, phi = 0) and their first and
# second theta-derivatives on Gauss-Legendre nodes in cos(theta).  Evaluated in BigFloat and rounded.
#   eth  sY_lm = +sqrt((l-s)(l+s+1)) s+1Y_lm,     ethbar sY_lm = -sqrt((l+s)(l-s+1)) s-1Y_lm,
#   eth eta = -(d_theta + i csc d_phi - s cot) eta,  ethbar eta = -(d_theta - i csc d_phi + s cot) eta,
# hence d_theta = -(eth + ethbar)/2 on any spin-weighted quantity.

using LinearAlgebra

function gauss_legendre(n)
    beta = [k / sqrt(4k^2 - 1) for k in 1:n-1]
    E = eigen(SymTridiagonal(zeros(n), beta))
    return E.values, 2 .* E.vectors[1, :] .^ 2
end

# Goldberg: sY_lm = (-1)^m sqrt((l+m)!(l-m)!(2l+1)/(4pi (l+s)!(l-s)!)) sin^{2l}(th/2)
#                   * sum_r C(l-s,r) C(l+s,r+s-m) (-1)^{l-r-s} cot^{2r+s-m}(th/2)
function sYlm_big(s, l, m, th::BigFloat)
    (l < abs(s) || l < abs(m)) && return big(0.0)
    pre = (-1)^m * sqrt(factorial(big(l + m)) * factorial(big(l - m)) * (2l + 1) /
                        (4 * big(pi) * factorial(big(l + s)) * factorial(big(l - s))))
    tot = big(0.0)
    ct = cot(th / 2)
    for r in 0:(l - s)
        k = r + s - m
        (k < 0 || k > l + s) && continue
        tot += binomial(big(l - s), big(r)) * binomial(big(l + s), big(k)) * (-1)^(l - r - s) * ct^(2r + s - m)
    end
    return pre * sin(th / 2)^(2l) * tot
end

eth_up(s, l) = (l < abs(s) || l < abs(s + 1)) ? 0.0 : sqrt((l - s) * (l + s + 1))       # coefficient of s+1Y
eth_dn(s, l) = (l < abs(s) || l < abs(s - 1)) ? 0.0 : -sqrt((l + s) * (l - s + 1))      # coefficient of s-1Y

"""
    SWSHTable(m, Lmax, nth; smax=2)

Values Y[s][l] (vector over nodes), first and second theta-derivatives dY, d2Y, for s = -smax..smax, l = 0..Lmax.
Derivatives use d_theta = -(eth + ethbar)/2, which needs spin weights up to smax + 2.
"""
struct SWSHTable
    m::Int
    Lmax::Int
    x::Vector{Float64}      # cos(theta) nodes
    w::Vector{Float64}      # Gauss weights (∫ dOmega = 2pi sum w)
    th::Vector{Float64}
    Y::Dict{Tuple{Int,Int},Vector{Float64}}
    dY::Dict{Tuple{Int,Int},Vector{Float64}}
    d2Y::Dict{Tuple{Int,Int},Vector{Float64}}
end

function SWSHTable(m, Lmax, nth; smax=2, prec=256)
    x, w = gauss_legendre(nth)
    th = acos.(x)
    Yb = Dict{Tuple{Int,Int},Vector{BigFloat}}()
    setprecision(BigFloat, prec) do
        thb = [acos(big(xi)) for xi in x]
        # x from eigen() is Float64; refine theta for consistency: use acos of the Float64 node
        for s in -(smax + 2):(smax + 2), l in 0:Lmax
            Yb[(s, l)] = [sYlm_big(s, l, m, t) for t in thb]
        end
    end
    zb(n) = zeros(BigFloat, n)
    getb(s, l) = get(Yb, (s, l), zb(nth))
    dYb(s, l) = -(eth_up(s, l) .* getb(s + 1, l) .+ eth_dn(s, l) .* getb(s - 1, l)) ./ 2
    Y = Dict{Tuple{Int,Int},Vector{Float64}}(); dY = empty(Y); d2Y = empty(Y)
    for s in -smax:smax, l in 0:Lmax
        Y[(s, l)] = Float64.(getb(s, l))
        dY[(s, l)] = Float64.(dYb(s, l))
        d2Y[(s, l)] = Float64.(-(eth_up(s, l) .* dYb(s + 1, l) .+ eth_dn(s, l) .* dYb(s - 1, l)) ./ 2)
    end
    return SWSHTable(m, Lmax, x, w, th, Y, dY, d2Y)
end

lmin_of(s, m) = max(abs(s), abs(m))
