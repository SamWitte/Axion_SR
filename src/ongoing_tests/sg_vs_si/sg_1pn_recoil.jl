# Recoil term of the self-gravity amplitude N1 x N1 -> N3 x {BH, Inf} at 1PN: the pieces that need the transition
# current (hydrogenic states, same machinery and units as sg_1pn.jl: sources s = S / (2 mu G), lengths in GM).
#
# In the frame of the hole the field feels the fictitious potential a_BH . x.  At 1PN (note, Section 3.7) the
# acceleration of the hole is
#     a_BH = grad v(0) + 4 d/dt v_i(0),
# the second term being gravitomagnetic induction by the transition current,  -4 i w_g int j(y)/|y| d^3y.
#   R   : s_R   = +alpha^2 (a . x) psi1,            a  = int rho(y) y/|y|^3 d^3y        (Newtonian recoil; check:
#                                                                                         A_R/A_N = A_NR/A_N - 1)
#   ind : s_ind = -4 i alpha Delta (A0 . x) psi1,   A0 = int j(y)/|y| d^3y              (induction)
# Only transitions with |m1 - m3| = 1 have a recoil.  Output: delta_R = A_R/A_N and delta_ind = A_ind/A_N.
#
# usage: julia sg_1pn_recoil.jl <S1> <S3> <BH|Inf> <alpha,alpha,...> <outdir>
#        -> <outdir>/<S1>_<S1>_<S3>_<dest>_<alpha>.dat :  alpha  Re,Im delta_R  Re,Im delta_ind
include(joinpath(@__DIR__, "sg_1pn.jl"))
using Printf

cvolint(F::Fld, r) = (@assert F.m == 0; 2pi * trapz(r .^ 2 .* (F.v * THw), r))          # complex volume integral

function point_recoil(st1, st3, to_inf, alph; a=0.95, rpts=3000)
    n1, l1, m1 = st1; n3, l3, m3 = st3
    mu = alph / GNew; M = 1.0
    m4 = 2m1 - m3; maxN = max(n1, n3); Mt = m1 - m3
    @assert abs(Mt) == 1 "no recoil: |m1 - m3| != 1"
    erg1 = complex(ergL(n1, l1, m1, mu, M, a; full=true) * GNew * M)
    erg3 = complex(ergL(n3, l3, m3, mu, M, a; full=true) * GNew * M)
    erg = 2erg1 - erg3
    @assert (real(erg) > alph) == to_inf
    rp = 1.0 + sqrt(1.0 - a^2); eps_fac = 1e-3
    rmax = to_inf ? 30.0 / real(sqrt(erg^2 - alph^2)) :
           2.0^(2.0 * maxN - 2 * (1 + maxN)) * gamma(2 + 2 * maxN) / alph^2 / factorial(2 * maxN - 1) * 10.0
    r = collect(10 .^ range(log10(rp * (1.0 + eps_fac)), log10(rmax), rpts))
    rmax_1 = (2.0^(2.0 * 3 - 2 * (1 + 3)) * gamma(2 + 2 * 3) / (0.03^2) / factorial(2 * 3 - 1) * 7.0)
    rmax_ratio = (2.0^(2.0 * maxN - 2 * (1 + maxN)) * gamma(2 + 2 * maxN) / alph^2 / factorial(big(2 * maxN - 1)) * 7.0) / rmax_1
    h_mve = Float64(0.2 / rmax_ratio * 10.0)

    s1 = hyd_state(n1, l1, m1, alph, r); s3 = hyd_state(n3, l3, m3, alph, r)
    Delta = bohr(n1, alph) - bohr(n3, alph)
    Lmax = l1 + l3 + 3
    src, aux = all_sources(s1, s3, alph, Delta, r, Lmax)
    # V . x = 1/2 (V+ x- + V- x+) + Vz z  with  V+- = Vx +- i Vy,  x+- = r sin(theta) e^{+-i phi}.
    # rho and j have azimuthal number Mt = -+1, so only the component V_(-Mt) of the integrals survives.
    xs_(s) = Fld(r * transpose(SINTH), s)                                    # r sin(theta) with azimuthal number s
    wgt(F, s) = Fld(F.v .* transpose(SINTH) ./ r .^ 2, F.m + s)              # F * y_(s) / |y|^3
    jc = Mt == -1 ? aux.jp : aux.jm                                          # the Cartesian-spherical component with m = 0
    @assert jc.m == 0
    acomp = cvolint(wgt(aux.rho, -Mt), r)                                    # a_(-Mt) = int rho y_(-Mt) / |y|^3
    A0comp = cvolint(rmul(jc, 1 ./ r), r)                                    # A0_(-Mt) = int j_(-Mt) / |y|
    xfac = xs_(Mt) * s1.psi                                                  # x_(Mt) psi1   (azimuthal number m4)
    sR = (alph^2 * 0.5 * acomp) * xfac
    sind = (-4im * alph * Delta * 0.5 * A0comp) * xfac

    lmin = max(abs(2l1 - l3), abs(m4)); iseven(lmin + 2l1 + l3) || (lmin += 1)
    T = [project(src[:N], r, a, lmin, m4), project(sR, r, a, lmin, m4), project(sind, r, a, lmin, m4)]
    amps, _, _ = gf_amplitudes(erg, a, alph, lmin, m4, r, T, rmax, to_inf; h_mve=h_mve, eps_fac=eps_fac)
    return amps[2] / amps[1], amps[3] / amps[1]
end

if abspath(PROGRAM_FILE) == @__FILE__
    S1, S3, dest = ARGS[1], ARGS[2], ARGS[3]
    st(s) = (parse(Int, s[1]), parse(Int, s[2]), parse(Int, s[3]))
    outdir = ARGS[5]; mkpath(outdir)
    for als in split(ARGS[4], ",")
        al = parse(Float64, als)
        dR, dind = point_recoil(st(S1), st(S3), dest == "Inf", al)
        open(joinpath(outdir, "$(S1)_$(S1)_$(S3)_$(dest)_$(als).dat"), "w") do io
            @printf(io, "%s %.8e %.8e %.8e %.8e\n", als, real(dR), imag(dR), real(dind), imag(dind))
        end
        @printf("%s x %s -> %s %s  alpha %s:  A_R/A_N = %+.5f%+.5fi   A_ind/A_N = %+.5e%+.5ei   (/alpha^2: %+.4f)\n", S1, S1, S3, dest, als,
                real(dR), imag(dR), real(dind), imag(dind), real(dind) / al^2)
        flush(stdout)
    end
end
