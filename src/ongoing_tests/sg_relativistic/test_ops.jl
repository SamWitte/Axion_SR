# Checks of the background and the pointwise operators:
#  (1) dG[h] = 0 for h = d g / dM (Kerr mass perturbation, stationary)
#  (2) dG[h] = 0 for pure gauge h = L_xi g with xi ∝ e^{-i w v + i m phi}
#  (3) -δ□[L_xi g] phi = □(xi·∂phi) - xi·∂(□phi)   (Lie derivative of the wave operator)
include("kerr.jl"); include("jets.jl")
using ForwardDiff
const D = ForwardDiff.derivative

# jet values of a tensor field hfun(i,j,r,th) (complex), as a 60-vector
function jetvals(hfun, r, th)
    J = zeros(ComplexF64, NJ)
    for (c, (i, j)) in enumerate(PAIRS)
        f(t, u) = hfun(i, j, r + t, th + u)
        J[jet(c, 1)] = f(0.0, 0.0)
        J[jet(c, 2)] = D(t -> f(t, 0.0), 0.0)
        J[jet(c, 3)] = D(u -> f(0.0, u), 0.0)
        J[jet(c, 4)] = D(t -> D(t2 -> f(t + t2, 0.0), 0.0), 0.0)
        J[jet(c, 5)] = D(t -> D(u -> f(t, u), 0.0), 0.0)
        J[jet(c, 6)] = D(u -> D(u2 -> f(0.0, u + u2), 0.0), 0.0)
    end
    return J
end
contract(Ejet, J) = [sum(Ejet[μ, ν, :] .* J) for μ in 1:4, ν in 1:4]

a = 0.7
for (r, th) in ((5.3, 0.9), (2.1 + 0.0im, 2.3), (40.0 + 25.0im, 1.4))
    bg = bg_point(a, r, th)
    # (1) mass perturbation, w = m = 0
    function hM(i, j, r, th)
        i > j && return hM(j, i, r, th)
        s = sin(th); c = cos(th); Sig = r^2 + a^2 * c^2
        (i, j) == (1, 1) && return 2r / Sig
        (i, j) == (1, 4) && return -2a * r * s^2 / Sig
        (i, j) == (4, 4) && return 2a^2 * r * s^4 / Sig
        return zero(Sig)
    end
    H, dH, ddH = jet_arrays(0.0, 0)
    D1, D2 = cov_derivs(bg, H, dH, ddH)
    E, dG, Z = einstein_jets(bg, D1, D2)
    J = jetvals(hM, r, th)
    println("r=$r th=$th  (1) max|dG[dM g]| = ", maximum(abs.(contract(dG, J))), "   (scale ", maximum(abs.(J)), ")")

    # (2) pure gauge with w, m
    w = 0.013; m = -1
    xi(λ, r, th) = (λ == 1 ? (1 + 0.3cos(th)) / r : λ == 2 ? sin(th)^2 * r / (r + 3) : λ == 3 ? 0.2 * sin(th) * cos(th) / r^2 : 0.5 / (r^2 + 1)) * (1 + 0.1im)
    dxi(k, λ, r, th) = k == 1 ? -im * w * xi(λ, r, th) : k == 4 ? im * m * xi(λ, r, th) :
                       k == 2 ? D(t -> xi(λ, r + t, th), 0.0) : D(u -> xi(λ, r, th + u), 0.0)
    dgc(k, i, j, r, th) = k == 2 ? D(t -> gcomp(i, j, a, r + t, th), 0.0) : k == 3 ? D(u -> gcomp(i, j, a, r, th + u), 0.0) : 0.0
    hgauge(i, j, r, th) = sum(xi(l, r, th) * dgc(l, i, j, r, th) + gcomp(l, j, a, r, th) * dxi(i, l, r, th) +
                              gcomp(i, l, a, r, th) * dxi(j, l, r, th) for l in 1:4)
    H, dH, ddH = jet_arrays(w, m)
    D1, D2 = cov_derivs(bg, H, dH, ddH)
    E, dG, Z = einstein_jets(bg, D1, D2)
    J = jetvals(hgauge, r, th)
    println("              (2) max|dG[L_xi g]| = ", maximum(abs.(contract(dG, J))), "   (|h| ~ ", maximum(abs.(J)), ", |E| ~ ", maximum(abs.(contract(E, J))), ")")

    # (3) scalar coupling: phi = f(r,th) e^{-i w1 v + i m1 phi}; check  -δ□[L_xi g]phi = □(L_xi phi) - L_xi(□phi)
    w1 = 0.29; m1 = 1
    f(r, th) = exp(-0.1r) * r * sin(th) * (1 + 0.2im * cos(th))
    # derivatives of a scalar F(r,th) e^{-i W v + i M phi}
    Frr_dummy(F, r, th) = zero(D(t -> F(r + t, th), 0.0))
    function sder(F, W, M, r, th)
        d1 = [-im * W * F(r, th), D(t -> F(r + t, th), 0.0), D(u -> F(r, th + u), 0.0), im * M * F(r, th)]
        k = [-im * W, 0, 0, im * M]
        dd = fill(zero(d1[1] + Frr_dummy(F, r, th)), 4, 4)
        Frr = D(t -> D(t2 -> F(r + t + t2, th), 0.0), 0.0); Frt = D(t -> D(u -> F(r + t, th + u), 0.0), 0.0)
        Ftt = D(u -> D(u2 -> F(r, th + u + u2), 0.0), 0.0)
        for i in 1:4, j in 1:4
            if i in (1, 4) && j in (1, 4)
                dd[i, j] = k[i] * k[j] * F(r, th)
            elseif i in (1, 4)
                dd[i, j] = k[i] * d1[j]
            elseif j in (1, 4)
                dd[i, j] = k[j] * d1[i]
            else
                dd[i, j] = (i, j) == (2, 2) ? Frr : (i, j) == (3, 3) ? Ftt : Frt
            end
        end
        return d1, dd
    end
    function boxs(F, W, M, r, th)   # □ of F e^{...}, divided by the exponential
        bgp = bg_point(a, r, th)
        d1, dd = sder(F, W, M, r, th)
        return sum(bgp.gi[μ, ν] * (dd[μ, ν] - sum(bgp.Gam[λ, μ, ν] * d1[λ] for λ in 1:4)) for μ in 1:4, ν in 1:4)
    end
    dphi, ddphi = sder(f, w1, m1, r, th)
    K = box_coupling(bg, H, D1, Z, dphi, ddphi)
    lhs = sum(K .* J)
    # xi·∂phi has frequency w + w1, azimuthal m + m1
    Lphi(r, th) = sum(xi(l, r, th) * sder(f, w1, m1, r, th)[1][l] for l in 1:4)
    boxf(r, th) = boxs(f, w1, m1, r, th)
    fd(F; h=1e-3) = (-F(2h) + 8F(h) - 8F(-h) + F(-2h)) / (12h)
    rhs = boxs(Lphi, w + w1, m + m1, r, th) - sum(xi(l, r, th) * (l == 1 ? -im * w1 * boxf(r, th) : l == 4 ? im * m1 * boxf(r, th) :
                                                               l == 2 ? fd(t -> boxf(r + t, th)) : fd(u -> boxf(r, th + u))) for l in 1:4)
    println("              (3) -δ□[L_xi g]phi = $(round(lhs, sigdigits=10))   [□,L_xi]phi = $(round(rhs, sigdigits=10))   rel diff = ", abs(lhs - rhs) / abs(rhs))
end
