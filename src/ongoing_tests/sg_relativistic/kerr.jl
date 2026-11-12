# Kerr background in ingoing Kerr coordinates x = (v, r, theta, phi~), G = M = 1:
#   ds^2 = -(1 - 2r/Sig) dv^2 + 2 dv dr - (4 a r s^2/Sig) dv dphi - 2 a s^2 dr dphi + Sig dth^2
#          + ((r^2+a^2)^2 - a^2 Del s^2) s^2/Sig dphi^2
# Regular at the future horizon r = r+.  r may be complex (exterior complex scaling); every background quantity
# is analytic in r, so radial derivatives are taken along the real direction.

using ForwardDiff

function gcomp(i, j, a, r, th)
    i > j && return gcomp(j, i, a, r, th)
    s = sin(th); c = cos(th)
    Sig = r^2 + a^2 * c^2
    if (i, j) == (1, 1)
        return -(1 - 2r / Sig)
    elseif (i, j) == (1, 2)
        return one(Sig)
    elseif (i, j) == (1, 4)
        return -2a * r * s^2 / Sig
    elseif (i, j) == (2, 4)
        return -a * s^2 * one(Sig)
    elseif (i, j) == (3, 3)
        return Sig
    elseif (i, j) == (4, 4)
        Del = r^2 - 2r + a^2
        return ((r^2 + a^2)^2 - a^2 * Del * s^2) * s^2 / Sig
    end
    return zero(Sig)
end

"""Background at one point: g, ginv, Gam[l,m,n] = Γ^l_{mn}, dGam[k][l,m,n] = ∂_k Γ^l_{mn} (k = 1..4, zero for v, phi)."""
struct BGPoint
    g::Matrix{ComplexF64}
    gi::Matrix{ComplexF64}
    Gam::Array{ComplexF64,3}
    dGam::Vector{Array{ComplexF64,3}}
end

function bg_point(a, r::Number, th::Real)
    D = ForwardDiff.derivative
    g = zeros(ComplexF64, 4, 4)
    dg = [zeros(ComplexF64, 4, 4) for _ in 1:4]                 # dg[k][i,j] = ∂_k g_ij
    ddg = [zeros(ComplexF64, 4, 4) for _ in 1:4, _ in 1:4]        # ddg[k,l][i,j]
    for i in 1:4, j in i:4
        f(t, u) = gcomp(i, j, a, r + t, th + u)
        v = f(0.0, 0.0)
        gr = D(t -> f(t, 0.0), 0.0)
        gt = D(u -> f(0.0, u), 0.0)
        grr = D(t -> D(t2 -> f(t + t2, 0.0), 0.0), 0.0)
        grt = D(t -> D(u -> f(t, u), 0.0), 0.0)
        gtt = D(u -> D(u2 -> f(0.0, u + u2), 0.0), 0.0)
        for (p, q) in ((i, j), (j, i))
            g[p, q] = v
            dg[2][p, q] = gr; dg[3][p, q] = gt
            ddg[2, 2][p, q] = grr; ddg[2, 3][p, q] = grt; ddg[3, 2][p, q] = grt; ddg[3, 3][p, q] = gtt
        end
    end
    gi = inv(g)
    dgi = [-gi * dg[k] * gi for k in 1:4]
    Gam = zeros(ComplexF64, 4, 4, 4)
    for l in 1:4, m in 1:4, n in 1:4
        Gam[l, m, n] = 0.5 * sum(gi[l, k] * (dg[m][k, n] + dg[n][k, m] - dg[k][m, n]) for k in 1:4)
    end
    dGam = [zeros(ComplexF64, 4, 4, 4) for _ in 1:4]
    for p in (2, 3), l in 1:4, m in 1:4, n in 1:4
        dGam[p][l, m, n] = 0.5 * sum(dgi[p][l, k] * (dg[m][k, n] + dg[n][k, m] - dg[k][m, n]) +
                                     gi[l, k] * (ddg[p, m][k, n] + ddg[p, n][k, m] - ddg[p, k][m, n]) for k in 1:4)
    end
    return BGPoint(g, gi, Gam, dGam)
end

rplus(a) = 1 + sqrt(1 - a^2)
rminus(a) = 1 - sqrt(1 - a^2)
OmegaH(a) = a / (2 * rplus(a))
# tortoise-like coordinates: dr*/dr = (r^2+a^2)/Del, dr#/dr = a/Del  (r > r+)
function rstar(a, r)
    rp, rm = rplus(a), rminus(a)
    return r + 2rp / (rp - rm) * log((r - rp) / 2) - 2rm / (rp - rm) * log((r - rm) / 2)
end
rsharp(a, r) = (rp = rplus(a); rm = rminus(a); a / (rp - rm) * log((r - rp) / (r - rm)))

# future-directed unit normal to the Kerr-Schild slices t_KS = v - r (timelike everywhere, regular at the horizon):
#   n_μ = -N ∂_μ t_KS = -N (1, -1, 0, 0),   N^-2 = (r^2 + 2r + a^2 cos^2)/Sigma
function ks_normal(a, r, th)
    c = cos(th); Sig = r^2 + a^2 * c^2
    N = sqrt(Sig / (r^2 + 2r + a^2 * c^2))
    return ComplexF64[-N, N, 0, 0]
end
