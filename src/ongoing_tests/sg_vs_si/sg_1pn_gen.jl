# First post-Newtonian corrections to the self-gravity amplitude for UNEQUAL incoming states,  N1 x N2 -> N3 x {BH, Inf}.
#
# Generalisation of sg_1pn.jl (N1 x N1 -> N3).  Every piece of the source of mode 4 is bilinear in a "spectator" state
# psi_s (the incoming state that turns into mode 4) and the transition density / current of the other incoming state
# psi_d with psi_3:  rho = conj(psi3) psi_d.  For two different incoming states the source is the sum of the two
# assignments,
#     s_X = s_X[spectator = 1, density = (3,2)] + s_X[spectator = 2, density = (3,1)],
# with Delta = (omega_d - omega_3) G M of the pair that forms the density.  For N1 = N2 this is twice the source of
# sg_1pn.jl, so the ratios delta_X = A_X / A_N are identical.
#
# NOTE (2026-11-19): the delta_X are corrections to a Newtonian amplitude written with canonical Schroedinger
# wavefunctions in harmonic coordinates.  The rate tables (and the relativistic exchange) use the Klein-Gordon mode
# functions in the Boyer-Lindquist radius, and that Newtonian amplitude differs at relative order alpha^2.  The sum of
# the delta_X is therefore NOT the correction of a tabulated rate; it also holds the hole fixed and has no pair
# (s-channel) term.  The complete 1PN correction in the variables of the tables: src/scripts/sg_1pn.
#
# Usage:  julia sg_1pn_gen.jl <S1> <S2> <S3> <BH|Inf> <alpha> <l4[,l4,...]> <outdir>
#             -> <outdir>/<S1>_<S2>_<S3>_<dest>_<alpha>.dat, one row per l4:
#                l4  Re,Im A_N (arbitrary common normalisation)  Re,Im of delta_GM, delta_hij, delta_T00, delta_ret, delta_BHc
#                (total 1PN correction of the amplitude: delta = GM + 2 hij + T00 + ret + BHc, as in plot_1pn.py)
#         julia sg_1pn_gen.jl test      integration-by-parts checks with unequal states
include(joinpath(@__DIR__, "sg_1pn.jl"))

function all_sources2(sp::State, sd::State, s3::State, alph, Delta, r, Lmax)
    c3 = cstate(s3)
    rho = c3.psi * sd.psi                                          # conj(psi3) psi_d
    Phi, Phir, Phith = coulomb(rho, r, Lmax)
    sN = (-alph^2) * (sp.psi * Phi)
    jr = (-0.5im) * (c3.psi * sd.gr - sd.psi * c3.gr)
    jth = (-0.5im) * (c3.psi * sd.gth - sd.psi * c3.gth)
    jph = (-0.5im) * (c3.psi * sd.gph - sd.psi * c3.gph)
    jp, jm, jz = cart(jr, jth, jph)
    Ap, _, _ = coulomb(jp, r, Lmax + 1); Am, _, _ = coulomb(jm, r, Lmax + 1); Az, _, _ = coulomb(jz, r, Lmax + 1)
    gp, gm, gz = cart(sp.gr, sp.gth, sp.gph)
    Adotg = 0.5 * (Ap * gm) + 0.5 * (Am * gp) + Az * gz
    sGM = (-4im) * Adotg + (2 * alph * Delta) * (sp.psi * Phi)
    gdg = graddot(c3, sd)
    PhiK, _, _ = coulomb(gdg, r, Lmax + 2)
    gradPhi_dot_g = Phir * sp.gr + Phith * sp.gth + Fld((im * Phi.m) .* (Phi.v ./ r) ./ SINTH', Phi.m) * sp.gph
    contact = sp.psi * sd.psi * c3.psi
    sS = (-1.0) * (sp.psi * PhiK) + gradPhi_dot_g + Phi * sp.lap - (6pi) * contact
    sT00 = 0.5 * (sS + (6pi) * contact)
    Psi = linear_kernel(rho, r, Lmax)
    sret = (alph^2 * Delta^2 / 2) * (sp.psi * Psi)
    Phi_rho_r, _, _ = coulomb(rmul(rho, 1 ./ r), r, Lmax)
    s3b = alph^2 * (rmul(sp.psi * Phi, 1 ./ r) + sp.psi * Phi_rho_r)
    if rho.m == 0
        s3b = s3b + alph^2 * rmul(volint(rmul(rho, 1 ./ r), r) * sp.psi, 1 ./ r)
    end
    return Dict(:N => sN, :GM => sGM, :hij => 0.5 * sS, :S => sS, :T00 => sT00, :ret => sret, :BHc => s3b),
           (; rho, Phi, jp, jm, jz, Ap, Am, Az, c3, gdg)
end

const KEYS1PN = [:N, :GM, :hij, :T00, :ret, :BHc]

function point_1pn_gen(st1, st2, st3, to_inf, alph, l4s; a=0.95, rpts=3000)
    n1, l1, m1 = st1; n2, l2, m2 = st2; n3, l3, m3 = st3
    mu = alph / GNew; M = 1.0
    m4 = m1 + m2 - m3; maxN = max(n1, n2, n3)
    en(n, l, m) = complex(ergL(n, l, m, mu, M, a; full=true) * GNew * M)
    erg = en(n1, l1, m1) + en(n2, l2, m2) - en(n3, l3, m3)
    @assert (real(erg) > alph) == to_inf
    rp = 1.0 + sqrt(1.0 - a^2); eps_fac = 1e-3
    rmax = to_inf ? 30.0 / real(sqrt(erg^2 - alph^2)) :
           2.0^(2.0 * maxN - 2 * (1 + maxN)) * gamma(2 + 2 * maxN) / alph^2 / factorial(2 * maxN - 1) * 10.0
    r = collect(10 .^ range(log10(rp * (1.0 + eps_fac)), log10(rmax), rpts))
    rmax_1 = (2.0^(2.0 * 3 - 2 * (1 + 3)) * gamma(2 + 2 * 3) / (0.03^2) / factorial(2 * 3 - 1) * 7.0)
    rmax_ratio = (2.0^(2.0 * maxN - 2 * (1 + maxN)) * gamma(2 + 2 * maxN) / alph^2 / factorial(big(2 * maxN - 1)) * 7.0) / rmax_1
    h_mve = Float64(0.2 / rmax_ratio * 10.0)

    sa = hyd_state(n1, l1, m1, alph, r); sb = hyd_state(n2, l2, m2, alph, r); s3 = hyd_state(n3, l3, m3, alph, r)
    srcA, _ = all_sources2(sa, sb, s3, alph, bohr(n2, alph) - bohr(n3, alph), r, l2 + l3 + 3)   # spectator 1, density (3,2)
    srcB, _ = all_sources2(sb, sa, s3, alph, bohr(n1, alph) - bohr(n3, alph), r, l1 + l3 + 3)   # spectator 2, density (3,1)
    out = Dict{Int,Any}()
    for l4 in l4s
        T = [project(srcA[k] + srcB[k], r, a, l4, m4) for k in KEYS1PN]
        amps, _, _ = gf_amplitudes(erg, a, alph, l4, m4, r, T, rmax, to_inf; h_mve=h_mve, eps_fac=eps_fac)
        out[l4] = (amps[1], Dict(k => amps[i] / amps[1] for (i, k) in enumerate(KEYS1PN)))
    end
    return out
end

# ---------------------------------------------------------------- self-test with unequal states (all bound, hydrogenic)
# Mode 4 is replaced by a bound hydrogenic state; for each assignment (spectator s, density (3,d)) the projection of
# the source,  ∫ conj(psi4) s_X d^3x,  is compared with the symmetric double integral evaluated directly.
function selftest_gen(sts, std, st3, st4; alph=0.2)
    nmx = maximum(s[1] for s in (sts, std, st3, st4))
    r = collect(10 .^ range(log10(1e-3), log10(40 * nmx^2 / alph^2), 6000))
    sp = hyd_state(sts..., alph, r); sd = hyd_state(std..., alph, r); s3 = hyd_state(st3..., alph, r); s4 = hyd_state(st4..., alph, r)
    Delta = bohr(std[1], alph) - bohr(st3[1], alph); Lmax = std[2] + st3[2] + 3
    src, aux = all_sources2(sp, sd, s3, alph, Delta, r, Lmax)
    c4 = cstate(s4)
    via_source(k) = volint(c4.psi * src[k], r)
    rho4s = c4.psi * sp.psi
    N_dir = -alph^2 * volint(rho4s * aux.Phi, r)
    jr = (-0.5im) * (c4.psi * sp.gr - sp.psi * c4.gr); jth = (-0.5im) * (c4.psi * sp.gth - sp.psi * c4.gth); jph = (-0.5im) * (c4.psi * sp.gph - sp.psi * c4.gph)
    jp4, jm4, jz4 = cart(jr, jth, jph)
    GM_dir = 4 * volint(0.5 * (jp4 * aux.Am) + 0.5 * (jm4 * aux.Ap) + jz4 * aux.Az, r)
    g4s = graddot(c4, sp)
    lap4s = c4.psi * sp.lap + sp.psi * c4.lap + 2.0 * g4s
    Tkk4s = g4s - 0.75 * lap4s
    Tkk3d_pot, _, _ = coulomb(aux.gdg - 0.75 * (aux.c3.psi * sd.lap + sd.psi * aux.c3.lap + 2.0 * aux.gdg), r, Lmax + 2)
    S_dir = -(volint(rho4s * Tkk3d_pot, r) + volint(Tkk4s * aux.Phi, r))
    println("selftest alpha=$alph: spectator $(sts), density ($(st3), $(std)), mode 4 -> bound $(st4)")
    for (k, dir) in [(:N, N_dir), (:GM, GM_dir), (:S, S_dir)]
        v = via_source(k)
        println("  $k  via source = $(round(v, sigdigits=8))   direct = $(round(dir, sigdigits=8))   rel.diff = $(round(abs(v - dir) / abs(dir), sigdigits=3))")
    end
end

parse_st(s) = occursin('-', s) ? Tuple(parse.(Int, split(s, '-'))) : (parse(Int, s[1:1]), parse(Int, s[2:2]), parse(Int, s[3:3]))

if abspath(PROGRAM_FILE) == @__FILE__
    if ARGS[1] == "test"
        selftest(alph=0.2)                                                   # the original equal-state check, for reference
        selftest_gen((2, 1, 1), (2, 1, 1), (3, 2, 2), (3, 2, 0); alph=0.2)   # must reproduce it
        selftest_gen((3, 2, 2), (5, 4, 4), (2, 1, 1), (6, 5, 5); alph=0.2)   # 322 x 544 -> 211, both assignments
        selftest_gen((5, 4, 4), (3, 2, 2), (2, 1, 1), (6, 5, 5); alph=0.2)
        selftest_gen((3, 2, 2), (4, 2, 2), (2, 1, 1), (4, 3, 3); alph=0.1)   # 322 x 422 -> 211
        selftest_gen((4, 2, 2), (3, 2, 2), (2, 1, 1), (4, 3, 3); alph=0.1)
    else
        S1, S2, S3, dest = ARGS[1:4]; alph = parse(Float64, ARGS[5]); l4s = parse.(Int, split(ARGS[6], ",")); outdir = ARGS[7]
        mkpath(outdir); t0 = time()
        d = point_1pn_gen(parse_st(S1), parse_st(S2), parse_st(S3), dest == "Inf", alph, l4s)
        open(joinpath(outdir, "$(S1)_$(S2)_$(S3)_$(dest)_$(ARGS[5]).dat"), "w") do io
            for l4 in l4s
                aN, dl = d[l4]
                println(io, join([l4, real(aN), imag(aN), vcat([[real(dl[k]), imag(dl[k])] for k in (:GM, :hij, :T00, :ret, :BHc)]...)...], " "))
                tot = dl[:GM] + 2 * dl[:hij] + dl[:T00] + dl[:ret] + dl[:BHc]
                println("$(S1) x $(S2) -> $(S3) x $(dest)  alpha=$(alph) l4=$(l4)  GM=$(round(dl[:GM], sigdigits=3)) hij=$(round(dl[:hij], sigdigits=3)) T00=$(round(dl[:T00], sigdigits=3)) ret=$(round(dl[:ret], sigdigits=3)) 3b=$(round(dl[:BHc], sigdigits=3))  total=$(round(tot, sigdigits=4))  [$(round(time() - t0, digits=1)) s]")
            end
        end
    end
end
