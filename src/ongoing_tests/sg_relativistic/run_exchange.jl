# Driver: exact one-graviton exchange amplitude vs Newtonian (fixed BH), Newtonian + BH recoil, and contact (SI).
# usage: julia run_exchange.jl <S1> <S3> <BH|Inf> <alpha> [spin] [Lextra] [drho]
include("exchange.jl"); include("hybrid.jl"); include("pinned.jl")
const SRC = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(SRC, "ongoing_tests", "sg_vs_si", "sg_vs_si.jl"))      # gf_amplitudes, find_im_part, constants
using Printf, DelimitedFiles
LinearAlgebra.BLAS.set_num_threads(2)

function run_channel(st1, st3, to_inf, alph; a=0.95, Lextra=6, drho=0.02, rpts=6000, verbose=true, kap0=0.0, kappow=0,
                     method=:pinned, R1fac=1.1, R2fac=3.5, frame_scale=1.0, nstrut=3, dkind=:d, odd=:auto, ord=8, csfac=4.0,
                     setup_kw...)
    n1, l1, m1 = st1; n3, l3, m3 = st3
    mu = alph / GNew
    g1 = find_im_part(mu, 1.0, a, n1, l1, m1; Ntot_force=5000, for_s_rates=true, return_both=true)
    g3 = find_im_part(mu, 1.0, a, n3, l3, m3; Ntot_force=5000, for_s_rates=true, return_both=true)
    md1, r1 = kerr_mode(n1, l1, m1, alph, a, ComplexF64(Float64(g1[1]) + im * Float64(g1[2])))
    md3, r3 = kerr_mode(n3, l3, m3, alph, a, ComplexF64(Float64(g3[1]) + im * Float64(g3[2])))
    verbose && @printf("modes: w1 = %.12f%+.2ei (res %.1e)   w3 = %.12f%+.2ei (res %.1e)\n", real(md1.w), imag(md1.w), abs(r1), real(md3.w), imag(md3.w), abs(r3))
    w4 = 2 * real(md1.w) - real(md3.w); m4 = 2m1 - m3
    @assert (w4 > alph) == to_inf
    lmin = max(abs(2l1 - l3), abs(m4)); iseven(lmin + 2l1 + l3) || (lmin += 1)
    l4 = lmin

    es = setup(md1, md3; Lmax=l1 + l3 + Lextra, drho=drho, ord=ord, csfac=csfac, setup_kw...)
    rc = real(es.rg.r[es.ireal])
    if abs(es.mg) == 1 && method == :pinned
        X, hA, acoef = solve_h_pinned(es; R1=R1fac * rc, R2=R2fac * rc, verbose=verbose, frame_scale=frame_scale, nstrut=nstrut, dkind=dkind, odd=odd)
        S4 = mode4_source_pinned(es, X, hA); Zrel = zeros(es.ireal)
    elseif abs(es.mg) == 1 && method == :hybrid
        hg = HybridGauge(a, es.wg, es.mg; R1=R1fac * rc, R2=R2fac * rc)
        X, beta = solve_h_hybrid(es, hg; verbose=verbose)
        S4, Zrel = mode4_source(es, X)
    else
        kap = r -> kap0 / (1 + abs(r))^kappow
        X, A, b = solve_h(es; verbose=verbose, kap=kap)
        S4, Zrel = mode4_source(es, X)
    end
    # gauge-condition check, weighted to where the source lives
    wts = [maximum(abs.(S4[i, :])) for i in 1:es.ireal]
    verbose && method != :pinned && @printf("  Lorenz condition |Z|/|∇h|: max over cloud region %.2e (median %.2e)\n", maximum(Zrel[wts .> 1e-3 * maximum(wts)]), sort(Zrel)[length(Zrel) ÷ 2])

    # fine BL grid and Green's function (as sg_vs_si.jl)
    rp = rplus(a); eps_fac = 1e-3; maxN = max(n1, n3)
    rmax = to_inf ? 30.0 / sqrt(w4^2 - alph^2) :
           2.0^(2.0 * maxN - 2 * (1 + maxN)) * gamma(2 + 2 * maxN) / alph^2 / factorial(2 * maxN - 1) * 10.0
    rlist = collect(10 .^ range(log10(rp * (1.0 + eps_fac)), log10(rmax), rpts))
    rmax_1 = (2.0^(2.0 * 3 - 2 * (1 + 3)) * gamma(2 + 2 * 3) / (0.03^2) / factorial(2 * 3 - 1) * 7.0)
    rmax_ratio = (2.0^(2.0 * maxN - 2 * (1 + maxN)) * gamma(2 + 2 * maxN) / alph^2 / factorial(big(2 * maxN - 1)) * 7.0) / rmax_1
    h_mve = Float64(0.2 / rmax_ratio * 10.0)

    T_SI, T_N, T_NR, sph4 = newtonian_sources(md1, md3, l4, m4, w4, rlist)
    rs, Tg = gr_radial_source(es, S4, l4, m4, w4, sph4)
    T_GR = interp_rho(es, rs, Tg, rlist)
    amps, _, _ = gf_amplitudes(complex(w4), a, alph, l4, m4, rlist, [T_SI, T_N, T_NR, T_GR], rmax, to_inf; h_mve=h_mve, eps_fac=eps_fac)
    aSI, aN, aNR, aGR = amps
    C_N = 64 * alph^4 * abs2(0.5 * aN) / abs2(0.5 * aSI)
    res = (; alph, a, l4, m4, w4, aSI, aN, aNR, aGR, gr_over_n=aGR / aN, nr_over_n=aNR / aN, gr_over_nr=aGR / aNR, C_N, si_over_n=aSI / aN,
             zmax=maximum(Zrel[wts .> 1e-3 * maximum(wts)]))
    verbose && @printf("  alpha=%.4f l4=%d: A_GR/A_N = %.6f%+.6fi (|.|^2 = %.5f)   A_NR/A_N = %.6f%+.6fi (|.|^2 = %.5f)   A_GR/A_NR |.|^2 = %.5f   C_N = %.4e\n",
                       alph, l4, real(res.gr_over_n), imag(res.gr_over_n), abs2(res.gr_over_n), real(res.nr_over_n), imag(res.nr_over_n),
                       abs2(res.nr_over_n), abs2(res.gr_over_nr), C_N)
    return res
end

if abspath(PROGRAM_FILE) == @__FILE__() && !isdefined(Main, :NO_MAIN)
    st(s) = (parse(Int, s[1]), parse(Int, s[2]), parse(Int, s[3]))
    a = length(ARGS) >= 5 ? parse(Float64, ARGS[5]) : 0.95
    Lx = length(ARGS) >= 6 ? parse(Int, ARGS[6]) : 6
    dr = length(ARGS) >= 7 ? parse(Float64, ARGS[7]) : 0.03
    k0 = length(ARGS) >= 8 ? parse(Float64, ARGS[8]) : 0.0
    kp = length(ARGS) >= 9 ? parse(Int, ARGS[9]) : 0
    t0 = time()
    run_channel(st(ARGS[1]), st(ARGS[2]), ARGS[3] == "Inf", parse(Float64, ARGS[4]); a=a, Lextra=Lx, drho=dr, kap0=k0, kappow=kp)
    println("  total time ", round(time() - t0, digits=1), " s")
end
