# Recoil term of the self-gravity amplitude N1 x N2 -> N3 x {BH, Inf} at 1PN: the pieces that need the transition
# current (hydrogenic states, same machinery, units and normalisation as sg_1pn_gen.jl: sources s = S / (2 mu G),
# lengths in GM, the sum of the two exchanges "spectator a, density (3, b)").
#
# In the frame of the hole the field feels the fictitious potential a_BH . x.  At 1PN (note, Section 3.7) the
# acceleration of the hole is
#     a_BH = grad v(0) + 4 d/dt v_i(0),
# the second term being gravitomagnetic induction by the transition current,  -4 i w_g int j(y)/|y| d^3y.
# For each exchange with |m_b - m_3| = 1:
#   R   : s_R   = +alpha^2 (a . x) psi_a,             a  = int rho(y) y/|y|^3 d^3y      (Newtonian recoil; check:
#                                                                                         A_R/A_N = A_NR/A_N - 1)
#   ind : s_ind = -4 i alpha Delta (A0 . x) psi_a,    A0 = int j(y)/|y| d^3y            (induction),  Delta = (w_b - w_3) GM
# Output: delta_R = A_R/A_N and delta_ind = A_ind/A_N, A_N the Newtonian amplitude with the hole held fixed.
#
# usage: julia sg_1pn_recoil.jl <S1> <S2> <S3> <BH|Inf> <alpha> <l4[,l4,...]> <outdir>
#        -> <outdir>/<S1>_<S2>_<S3>_<dest>_<alpha>.dat, one row per l4:  l4  Re,Im delta_R  Re,Im delta_ind
include(joinpath(@__DIR__, "sg_1pn_gen.jl"))
using Printf

cvolint(F::Fld, r) = (@assert F.m == 0; 2pi * trapz(r .^ 2 .* (F.v * THw), r))          # complex volume integral

# recoil and induction sources of one exchange: spectator sp, transition density / current of (3, d) in aux
function recoil_sources(sp::State, aux, alph, Delta, r)
    Mt = aux.rho.m                                                           # m_d - m_3
    abs(Mt) == 1 || return nothing
    # V . x = 1/2 (V+ x- + V- x+) + Vz z  with  V+- = Vx +- i Vy,  x+- = r sin(theta) e^{+-i phi}.
    # rho and j have azimuthal number Mt = -+1, so only the component V_(-Mt) of the integrals survives.
    jc = Mt == -1 ? aux.jp : aux.jm                                          # the Cartesian-spherical component with m = 0
    acomp = cvolint(Fld(aux.rho.v .* transpose(SINTH) ./ r .^ 2, aux.rho.m - Mt), r)     # a_(-Mt) = int rho y_(-Mt) / |y|^3
    A0comp = cvolint(rmul(jc, 1 ./ r), r)                                                # A0_(-Mt) = int j_(-Mt) / |y|
    xfac = Fld(r * transpose(SINTH), Mt) * sp.psi                                        # x_(Mt) psi_a
    return (alph^2 * 0.5 * acomp) * xfac, (-4im * alph * Delta * 0.5 * A0comp) * xfac
end

function point_recoil(st1, st2, st3, to_inf, alph, l4s; a=0.95, rpts=3000)
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
    DA = bohr(n2, alph) - bohr(n3, alph); DB = bohr(n1, alph) - bohr(n3, alph)
    srcA, auxA = all_sources2(sa, sb, s3, alph, DA, r, l2 + l3 + 3)              # spectator 1, density (3,2)
    srcB, auxB = all_sources2(sb, sa, s3, alph, DB, r, l1 + l3 + 3)              # spectator 2, density (3,1)
    zero_ = 0.0 * srcA[:N]
    recA = recoil_sources(sa, auxA, alph, DA, r); recB = recoil_sources(sb, auxB, alph, DB, r)
    sR = (recA === nothing ? zero_ : recA[1]) + (recB === nothing ? zero_ : recB[1])
    sind = (recA === nothing ? zero_ : recA[2]) + (recB === nothing ? zero_ : recB[2])
    out = Dict{Int,Any}()
    for l4 in l4s
        T = [project(srcA[:N] + srcB[:N], r, a, l4, m4), project(sR, r, a, l4, m4), project(sind, r, a, l4, m4)]
        amps, _, _ = gf_amplitudes(erg, a, alph, l4, m4, r, T, rmax, to_inf; h_mve=h_mve, eps_fac=eps_fac)
        out[l4] = (amps[2] / amps[1], amps[3] / amps[1])
    end
    return out, recA !== nothing || recB !== nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    S1, S2, S3, dest = ARGS[1:4]; alph = parse(Float64, ARGS[5]); l4s = parse.(Int, split(ARGS[6], ",")); outdir = ARGS[7]
    mkpath(outdir)
    d, dip = point_recoil(parse_st(S1), parse_st(S2), parse_st(S3), dest == "Inf", alph, l4s)
    open(joinpath(outdir, "$(S1)_$(S2)_$(S3)_$(dest)_$(ARGS[5]).dat"), "w") do io
        for l4 in l4s
            dR, dind = d[l4]
            @printf(io, "%d %.8e %.8e %.8e %.8e\n", l4, real(dR), imag(dR), real(dind), imag(dind))
            @printf("%s x %s -> %s %s  alpha %s l4 %d:  A_R/A_N = %+.5f%+.5fi   A_ind/A_N = %+.5e%+.5ei   (/alpha^2: %+.4f)%s\n", S1, S2, S3, dest, ARGS[5], l4,
                    real(dR), imag(dR), real(dind), imag(dind), real(dind) / alph^2, dip ? "" : "   [no dipole exchange]")
        end
    end
end
