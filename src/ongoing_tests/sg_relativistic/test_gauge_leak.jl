# Does a pure-gauge metric perturbation h = L_zeta g give zero Inf amplitude through gf_amplitudes?
include("exchange.jl")
const SRC = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(SRC, "ongoing_tests", "sg_vs_si", "sg_vs_si.jl"))
using Printf
const FD = ForwardDiff.derivative
alph = 0.2; a = 0.95; mu = alph / GNew
g1 = find_im_part(mu, 1.0, a, 3, 2, 2; Ntot_force=5000, for_s_rates=true, return_both=true)
g3 = find_im_part(mu, 1.0, a, 2, 1, 1; Ntot_force=5000, for_s_rates=true, return_both=true)
md1, _ = kerr_mode(3, 2, 2, alph, a, ComplexF64(Float64(g1[1]) + im * Float64(g1[2])))
md3, _ = kerr_mode(2, 1, 1, alph, a, ComplexF64(Float64(g3[1]) + im * Float64(g3[2])))
es = setup(md1, md3; Lmax=9, drho=0.03)
w = es.wg; m = es.mg
# zeta: translation-like near the hole (x+iy direction ~ sin θ e^{iφ}), windowed at large r
win(r) = exp(-(r / 400)^2)
zeta(λ, r, th) = (λ == 2 ? sin(th) : λ == 3 ? cos(th) / r : λ == 4 ? im / r : 0.0) * win(r)   # ∂_x+i∂_y in flat space
dz(k, λ, r, th) = k == 1 ? -im * w * zeta(λ, r, th) : k == 4 ? im * m * zeta(λ, r, th) :
                  k == 2 ? FD(t -> zeta(λ, r + t, th), 0.0) : FD(u -> zeta(λ, r, th + u), 0.0)
dgc(k, i, j, r, th) = k == 2 ? FD(t -> gcomp(i, j, a, r + t, th), 0.0) : k == 3 ? FD(u -> gcomp(i, j, a, r, th + u), 0.0) : 0.0
hg(i, j, r, th) = sum(zeta(l, r, th) * dgc(l, i, j, r, th) + gcomp(l, j, a, r, th) * dz(i, l, r, th) + gcomp(i, l, a, r, th) * dz(j, l, r, th) for l in 1:4)
function jetvals(hfun, r, th)
    J = zeros(ComplexF64, NJ)
    for (c, (i, j)) in enumerate(PAIRS)
        f(t, u) = hfun(i, j, r + t, th + u)
        J[jet(c, 1)] = f(0.0, 0.0); J[jet(c, 2)] = FD(t -> f(t, 0.0), 0.0); J[jet(c, 3)] = FD(u -> f(0.0, u), 0.0)
        J[jet(c, 4)] = FD(t -> FD(t2 -> f(t + t2, 0.0), 0.0), 0.0); J[jet(c, 5)] = FD(t -> FD(u -> f(t, u), 0.0), 0.0)
        J[jet(c, 6)] = FD(u -> FD(u2 -> f(0.0, u + u2), 0.0), 0.0)
    end
    J
end
Hj = jet_arrays(w, m)
w1 = real(md1.w); S1, dS1, d2S1 = es.S1
nth = length(es.ab.T.th)
S4 = zeros(ComplexF64, es.ireal, nth)
for i in 1:es.ireal
    r = real(es.rg.r[i]); f1 = radial_f(md1, r)
    for (j, th) in enumerate(es.ab.T.th)
        bg = bg_point(a, r, th)
        D1, D2 = cov_derivs(bg, Hj...)
        _, _, Z = einstein_jets(bg, D1, D2)
        d1, dd = scalar_jets(f1..., S1[j], dS1[j], d2S1[j], w1, md1.m)
        K = box_coupling(bg, Hj[1], D1, Z, d1, dd)
        S4[i, j] = sum(K .* jetvals(hg, r, th))
    end
end
w4 = 2w1 - real(md3.w); m4 = 2md1.m - md3.m; l4 = 3
rp = rplus(a); rmax = 30.0 / sqrt(w4^2 - alph^2)
rlist = collect(10 .^ range(log10(rp * 1.001), log10(rmax), 6000))
T_SI, T_N, T_NR, sph4 = newtonian_sources(md1, md3, l4, m4, w4, rlist)
rs, Tg = gr_radial_source(es, S4, l4, m4, w4, sph4)
Tgauge = interp_rho(es, rs, Tg, rlist)
maxN = 3
rmax_1 = (2.0^(2.0 * 3 - 2 * (1 + 3)) * gamma(2 + 2 * 3) / (0.03^2) / factorial(2 * 3 - 1) * 7.0)
rmax_ratio = (2.0^(2.0 * maxN - 2 * (1 + maxN)) * gamma(2 + 2 * maxN) / alph^2 / factorial(big(2 * maxN - 1)) * 7.0) / rmax_1
h_mve = Float64(0.2 / rmax_ratio * 10.0)
for (hm, ef) in ((h_mve, 1e-3), (4h_mve, 1e-3), (h_mve, 1e-5), (4h_mve, 1e-6))
    rl = collect(10 .^ range(log10(rp * (1 + ef)), log10(rmax), 6000))
    Tg2 = interp_rho(es, rs, Tg, rl)
    Tn2 = newtonian_sources(md1, md3, l4, m4, w4, rl)[2]
    amps, _, _ = gf_amplitudes(complex(w4), a, alph, l4, m4, rl, [Tn2, Tg2], rmax, true; h_mve=hm, eps_fac=ef)
    @printf("h_mve=%.3g eps=%.0e:  A_gauge/A_N = %.3e   (|source| near horizon / in cloud: %.2e)\n", hm, ef, abs(amps[2] / amps[1]),
            maximum(abs.(Tg[1:10])) / maximum(abs.(Tg)))
end
