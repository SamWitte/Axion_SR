# Self-gravity vs self-interaction induced 2 -> 2 level transitions  N1 x N1 -> N3 x {BH, Inf}.
#
# Both interactions are number-conserving quartic couplings of the NR field,
#     H_4 = 1/2 ∫∫ |psi(x)|^2 K(x-y) |psi(y)|^2 ,
#     K_SI = -(lambda / 8 mu^2) delta^3(x-y),   lambda = mu^2 / f_a^2   (axion cosine)
#     K_SG = -G mu^2 / |x-y|,                   G = 1 / M_pl^2          (Newtonian h_00, integrated out)
# so the amplitude for N1 + N1 -> N3 + (mode 4) is, for either kernel,
#     A ∝ ∫∫ conj(psi_4(x)) psi_N1(x) K(x-y) conj(psi_N3(y)) psi_N1(y)
# with identical combinatorics.  Mode 4 is projected out with the Green's function of gf_radial
# (solve_sr_rates.jl, copied below for both the BH and the Inf branch), so the ratio depends only on
# the source:
#     SI:  T_SI(r) = R1 R1 conj(R3) x (angular overlap, Sigma-weighted)
#     SG:  T_SG(r) = R1(r) x Phi(r),  Phi = multipole potential of the transition density
#          R1 conj(R3) Y_{l1 m1} conj(Y_{l3 m3}) = sum_L c_L Y_{L, m1-m3}.
# In units r -> r/(GM) the physical amplitude ratio is
#     A_SG / A_SI = 8 (f_a / M_pl)^2 alpha^2 * J_SG / J_SI ,
# hence  Gamma_SG / Gamma_SI = C(alpha) (f_a / M_pl)^4 with C = 64 alpha^4 sum_l|J_SG|^2 / sum_l|J_SI|^2.
# Occupation numbers enter both rates identically (N1^2 (N3 + 1)) and cancel in C.
#
# Usage:  julia sg_vs_si.jl <S1> <S3> <BH|Inf> <i_start> <i_stop> [rel|nonrel] [spin]
#   e.g.  julia sg_vs_si.jl 322 544 BH 1 19
#   indices refer to the alpha grid of rate_sve/<S1>_<S1>_<S3>_<dest>_LvrHc_.dat

const SRC = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(SRC, "Core/constants.jl"))
include(joinpath(SRC, "solve_sr_rates.jl"))
using LinearAlgebra, DelimitedFiles

# ---------------------------------------------------------------- angular integrals
function gauss_legendre(n)
    beta = [k / sqrt(4k^2 - 1) for k in 1:n-1]
    E = eigen(SymTridiagonal(zeros(n), beta))
    return E.values, 2 .* E.vectors[1, :] .^ 2
end
const GLx, GLw = gauss_legendre(100)

# ∫ dΩ  Π_i [Y_{l_i m_i} or conj] * w(θ);  factors = [(l, m, conj::Bool), ...]   (spherical harmonics)
function ang_integral(factors; cos2=false)
    msum = sum(c ? -m : m for (l, m, c) in factors)
    msum == 0 || return 0.0
    tot = 0.0
    for (x, w) in zip(GLx, GLw)
        th = acos(x)
        v = 1.0 + 0im
        for (l, m, c) in factors
            y = sphericalY(l, m, th, 0.0)
            v *= c ? conj(y) : y
        end
        tot += w * real(v) * (cos2 ? x^2 : 1.0)
    end
    return 2pi * tot
end
@assert abs(ang_integral([(3, -1, false), (3, -1, true)]) - 1) < 1e-10   # orthonormality, negative m

# ---------------------------------------------------------------- Green's function (copied from gf_radial)
# Amplitude of the (l, m) mode at the horizon (to_inf=false) or at r_max (to_inf=true) for each source.
function gf_amplitudes(erg, a, alph, l, m, rlist, sources, rmax, to_inf; h_mve=1.0, eps_fac=1e-3)
    rp = 1.0 + sqrt(1.0 - a^2)
    rmm = 1.0 - sqrt(1.0 - a^2)
    itps = [(LinearInterpolation(log10.(rlist), real.(s), extrapolation_bc=Line()),
             LinearInterpolation(log10.(rlist), imag.(s), extrapolation_bc=Line())) for s in sources]

    gam = im * a * sqrt(erg^2 - alph^2)
    r_list_map = 10 .^ LinRange(log10(rp * (1.0 + eps_fac)), log10(rmax), 100000)
    rout_star = r_list_map .+ 2.0 .* rp ./ (rp .- rmm) .* log.((r_list_map .- rp) ./ 2.0) .- 2.0 .* rmm .* log.((r_list_map .- rmm) ./ 2.0) ./ (rp .- rmm)
    itp_rrstar = LinearInterpolation(rout_star, r_list_map, extrapolation_bc=Line())
    itp_rrstar_inv = LinearInterpolation(r_list_map, rout_star, extrapolation_bc=Line())
    h_step = rp / h_mve

    LLM = l * (l + 1)
    LLM += (-1 + 2 * l * (l + 1) - 2 * m^2) * gam^2 / (-3 + 4 * l * (l + 1))
    LLM += ((l - abs(m) - 1 * (l - abs(m)) * (l + abs(m)) * (l + abs(m) - 1)) / ((-3 + 2 * l) * (2 * l - 1)^2) - (l + 1 - abs(m)) * (2 * l - abs(m)) * (l + abs(m) + 1) * (2 + l + abs(m)) / ((3 + 2 * l)^2 * (5 + 2 * l))) * gam^4 / (2 * (1 + 2 * l))
    LLM += (4 * ((-1 + 4 * m^2) * (l * (1 + l) * (121 + l * (1 + l) * (213 + 8 * l * (1 + l) * (-37 + 10 * l * (1 + l)))) - 2 * l * (1 + l) * (-137 + 56 * l * (1 + l) * (3 + 2 * l * (1 + l))) * m^2 + (705 + 8 * l * (1 + l) * (125 + 18 * l * (1 + l))) * m^4 - 15 * (1 + 46 * m^2))) * gam^6) / ((-5 + 2 * l) * (-3 + 2 * l) * (5 + 2 * l) * (7 + 2 * l) * (-3 + 4 * l * (1 + l))^5)

    Vpot(r_input) = begin
        delt = r_input^2 - 2 * r_input + a^2
        ff = r_input^2 + a^2
        delt * alph^2 / ff + delt * (LLM + a^2 * (erg^2 - alph^2)) / ff^2 + delt * (3 * r_input^2 - 4 * r_input + a^2) / ff^3 - 3 * delt^2 * r_input^2 / ff^4 + 2 * a * m * erg / ff - 2 * a * m * erg * delt / ff^2 - a^2 * m^2 / ff^2
    end

    # solution #1: fixed at r_max (decaying for BH channels, outgoing for Inf channels), integrated inwards
    rr = itp_rrstar_inv(rmax)
    kap = sqrt(alph^2 - erg^2)
    start(r) = exp(-kap * r) * sqrt(itp_rrstar(r)^2 + a^2) / itp_rrstar(r)
    while abs(start(rr)) <= 1e-200
        rr /= 2.0
    end
    outWF = ComplexF64[start(rr)]; rvals = Float64[rr]
    rr -= h_step
    push!(outWF, start(rr)); push!(rvals, rr)
    rr -= h_step
    idx = 2
    while true
        r_input = itp_rrstar(rvals[idx])
        push!(outWF, 2 * outWF[idx] - outWF[idx-1] + h_step^2 * (Vpot(r_input) - erg^2) * outWF[idx])
        push!(rvals, rr)
        rr -= h_step
        idx += 1
        itp_rrstar(rr) < rp * (1.0 + eps_fac) && break
    end
    rvals = reverse(rvals); outWF = reverse(outWF)

    # solution #2: ingoing at the horizon, integrated outwards
    rr += h_step
    omega_H = a / (rp^2 + a^2)
    hstart(r) = exp(-im * (erg - m * omega_H) * r) * sqrt(itp_rrstar(r)^2 + a^2)
    outWF_fw = ComplexF64[hstart(rr)]
    rr += h_step
    push!(outWF_fw, hstart(rr))
    idx = 2
    while length(outWF_fw) < length(rvals)
        r_input = itp_rrstar(rvals[idx])
        push!(outWF_fw, 2 * outWF_fw[idx] - outWF_fw[idx-1] + h_step^2 * (Vpot(r_input) - erg^2) * outWF_fw[idx])
        idx += 1
    end
    rphys = itp_rrstar.(rvals)
    outWF ./= sqrt.(rphys .^ 2 .+ a^2)
    outWF_fw ./= sqrt.(rphys .^ 2 .+ a^2)

    midP = to_inf ? length(rvals) - 10 : clamp(Int(round(length(rvals) / 4)), 3, length(rvals) - 3)
    wronk = (outWF_fw[midP] * (outWF[midP+1] - outWF[midP-1]) - outWF[midP] * (outWF_fw[midP+1] - outWF_fw[midP-1])) / (rphys[midP+1] - rphys[midP-1])
    wronk *= (rphys[midP]^2 - 2 * rphys[midP] + a^2) * GNew   # M = 1, as in gf_radial

    sel = rphys .> 1.01 * rp
    amps = ComplexF64[]
    for (iR, iI) in itps
        Tmm = iR.(log10.(rphys)) .+ im .* iI.(log10.(rphys))
        if to_inf
            push!(amps, outWF[end] * trapz(outWF_fw .* Tmm, rphys) / wronk)
        else
            push!(amps, trapz(outWF[sel] .* Tmm[sel], rphys[sel]) * outWF_fw[sel][1] / wronk)
        end
    end
    return amps, rvals[end], itp_rrstar
end

# ---------------------------------------------------------------- one alpha point
function compare_point(st1, st3, to_inf, alph; a=0.95, nonrel=false, rpts=4000, Npts_Bnd=2000, Ntot_safe=30000, BHlmax=2)
    M = 1.0
    mu = alph / (M * GNew)
    n1, l1, m1 = st1
    n3, l3, m3 = st3
    m4 = 2m1 - m3
    maxN = max(n1, n3)

    if nonrel
        erg1 = complex(ergL(n1, l1, m1, mu, M, a; full=true) * GNew * M)
        erg3 = complex(ergL(n3, l3, m3, mu, M, a; full=true) * GNew * M)
    else
        e1 = find_im_part(mu, M, a, n1, l1, m1; Ntot_force=Ntot_safe, for_s_rates=true, return_both=true)
        e3 = find_im_part(mu, M, a, n3, l3, m3; Ntot_force=Ntot_safe, for_s_rates=true, return_both=true)
        erg1 = ComplexF64(e1[1] + im * e1[2]); erg3 = ComplexF64(e3[1] + im * e3[2])
    end
    erg = 2 * erg1 - erg3                    # energy of mode 4
    @assert (real(erg) > alph) == to_inf "energy of mode 4 inconsistent with destination"

    rp = 1.0 + sqrt(1.0 - a^2)
    eps_fac = 1e-3
    if to_inf
        rmax = 30.0 / real(sqrt(erg^2 - alph^2))
    else
        rmax = maxN < 10 ? 2.0^(2.0 * maxN - 2 * (1 + maxN)) * gamma(2 + 2 * maxN) / alph^2 / factorial(2 * maxN - 1) * 10.0 :
                           2.0^(2.0 * maxN - 2 * (1 + maxN)) * 4 * maxN^2 / alph^2 * 10.0
    end
    rlist = 10 .^ range(log10(rp * (1.0 + eps_fac)), log10(rmax), rpts)

    # h_mve exactly as in Compute_all_rates.jl
    rmax_1 = (2.0^(2.0 * 3 - 2 * (1 + 3)) * gamma(2 + 2 * 3) / (0.03^2) / factorial(2 * 3 - 1) * 7.0)
    rmax_ratio = (2.0^(2.0 * maxN - 2 * (1 + maxN)) * gamma(2 + 2 * maxN) / alph^2 / factorial(big(2 * maxN - 1)) * 7.0) / rmax_1
    h_mve = Float64(0.2 / rmax_ratio * 10.0)

    if nonrel
        rf1 = complex.(Float64.(radial_bound_NR(n1, l1, m1, mu, M, rlist)))
        rf3 = complex.(Float64.(radial_bound_NR(n3, l3, m3, mu, M, rlist)))
    else
        interp(rl, rf) = begin
            x = log10.(Float64.(real.(rl)))
            iR = LinearInterpolation(x, Float64.(real.(rf)), extrapolation_bc=Line())
            iI = LinearInterpolation(x, Float64.(imag.(rf)), extrapolation_bc=Line())
            iR.(log10.(rlist)) .+ im .* iI.(log10.(rlist))
        end
        rl1, r1, _ = solve_radial(mu, M, a, n1, l1, m1; rpts=Npts_Bnd, return_erg=true, Ntot_safe=Ntot_safe, use_heunc=true, pre_compute_erg=erg1)
        rl3, r3, _ = solve_radial(mu, M, a, n3, l3, m3; rpts=Npts_Bnd, return_erg=true, Ntot_safe=Ntot_safe, use_heunc=true, pre_compute_erg=erg3)
        rf1 = interp(rl1, r1); rf3 = interp(rl3, r3)
        # beyond the solver's grid the bound states are exponentially small; zero them rather than extrapolate
        rf1[rlist .> maximum(Float64.(real.(rl1)))] .= 0
        rf3[rlist .> maximum(Float64.(real.(rl3)))] .= 0
    end

    preFac = 0.5                             # identical incoming states (as gf_radial)
    unitMatch = -1.0 / (2 * alph)^(3 / 2)

    # transition-density potential: Phi(r) = sum_L 4pi/(2L+1) c_L Q_L(r) Y_{L,M}
    M_ = m1 - m3
    rho = rf1 .* conj.(rf3)
    Ls = [L for L in abs(l1 - l3):(l1 + l3) if abs(M_) <= L && iseven(l1 + l3 + L)]
    cL = Dict(L => ang_integral([(l1, m1, false), (l3, m3, true), (L, M_, true)]) for L in Ls)
    function Qpot(L)
        f_in = rho .* rlist .^ (L + 2)
        f_out = rho .* rlist .^ (1 - L)
        cin = zeros(ComplexF64, rpts); cout = zeros(ComplexF64, rpts)
        for i in 2:rpts
            cin[i] = cin[i-1] + 0.5 * (f_in[i] + f_in[i-1]) * (rlist[i] - rlist[i-1])
        end
        for i in rpts-1:-1:1
            cout[i] = cout[i+1] + 0.5 * (f_out[i] + f_out[i+1]) * (rlist[i+1] - rlist[i])
        end
        return rlist .^ (-(L + 1)) .* cin .+ rlist .^ L .* cout
    end
    QL = Dict(L => Qpot(L) for L in Ls)

    # outgoing multipoles: lowest one as in gf_radial; Inf channels also get l+2 as a convergence check
    lmin = max(abs(2l1 - l3), abs(m4))
    iseven(lmin + 2l1 + l3) || (lmin += 1)
    l4s = to_inf ? [lmin, lmin + 2] : collect(lmin:2:max(BHlmax, lmin))

    res = []
    rend = 0.0
    for l4 in l4s
        CG  = ang_integral([(l1, m1, false), (l1, m1, false), (l3, m3, true), (l4, m4, true)])
        CG2 = ang_integral([(l1, m1, false), (l1, m1, false), (l3, m3, true), (l4, m4, true)]; cos2=true)
        T_SI = preFac * unitMatch .* rf1 .* rf1 .* conj.(rf3) .* (CG .* rlist .^ 2 .+ CG2 * a^2)
        T_SG = zeros(ComplexF64, rpts)
        for L in Ls
            AL = ang_integral([(l1, m1, false), (L, M_, false), (l4, m4, true)])
            BL = ang_integral([(l1, m1, false), (L, M_, false), (l4, m4, true)]; cos2=true)
            T_SG .+= (4pi / (2L + 1)) * cL[L] .* QL[L] .* (AL .* rlist .^ 2 .+ BL * a^2)
        end
        T_SG .*= preFac * unitMatch .* rf1
        (aSI, aSG), rv_end, itp_rrstar = gf_amplitudes(erg, a, alph, l4, m4, rlist, [T_SI, T_SG], rmax, to_inf; h_mve=h_mve, eps_fac=eps_fac)
        rend = itp_rrstar(rv_end)               # flux radius (matches fixed gf_radial; cancels in C)
        push!(res, (l4, aSI, aSG))
    end
    amp2_SI_low = abs2(res[1][2]); amp2_SG_low = abs2(res[1][3])
    amp2_SI = sum(abs2(r[2]) for r in res); amp2_SG = sum(abs2(r[3]) for r in res)

    # SI rate in the units of gf_radial (Gamma/mu at f_a = M_pl), for validation against the tables:
    # BH channels sum l <= BHlmax (as production); Inf channels keep only the lowest l (as production)
    lam = (mu / (M_pl * 1e9))^2
    if to_inf
        kk = real(sqrt(erg^2 - alph^2))
        rate_SI = 2 * alph * kk * amp2_SI_low * rend^2 * lam^2 / mu^2 * (GNew * M^2 * M_to_eV)^2
    else
        kfac = m4 == 0 ? real(erg) : 1.0
        rate_SI = 4 * alph * kfac * (1 + sqrt(1 - a^2)) * amp2_SI * lam^2 / mu^2 * (GNew * M^2 * M_to_eV)^2
    end

    C = 64 * alph^4 * amp2_SG / amp2_SI                 # all outgoing multipoles computed
    C_low = 64 * alph^4 * amp2_SG_low / amp2_SI_low     # lowest multipole only
    phase = angle(res[1][3] / res[1][2])
    return (; alph, a, erg1, erg3, erg, rate_SI, C, C_low, phase,
              xfrac_SI=1 - amp2_SI_low / amp2_SI, xfrac_SG=1 - amp2_SG_low / amp2_SG, l4s)
end

# ---------------------------------------------------------------- driver
if abspath(PROGRAM_FILE) == @__FILE__
    S1, S3, dest = ARGS[1], ARGS[2], ARGS[3]
    st(s) = (parse(Int, s[1]), parse(Int, s[2]), parse(Int, s[3]))
    to_inf = dest == "Inf"
    tab = readdlm(joinpath(SRC, "rate_sve", "$(S1)_$(S1)_$(S3)_$(dest)_LvrHc_.dat"))
    i0 = parse(Int, ARGS[4]); i1 = parse(Int, ARGS[5])
    nonrel = length(ARGS) >= 6 && ARGS[6] == "nonrel"
    spin = length(ARGS) >= 7 ? parse(Float64, ARGS[7]) : 0.95
    outdir = joinpath(@__DIR__, "out", "$(S1)_$(S1)_$(S3)_$(dest)"); mkpath(outdir)
    for i in i0:min(i1, size(tab, 1))
        t0 = time()
        res = try
            compare_point(st(S1), st(S3), to_inf, tab[i, 1]; a=spin, nonrel=nonrel)
        catch e
            println("i=$i alpha=$(tab[i,1]) FAILED: ", sprint(showerror, e)); flush(stdout)
            continue
        end
        stag = spin == 0.95 ? "" : "a$(spin)_"
        fn = joinpath(outdir, (nonrel ? "NR_" : "REL_") * stag * lpad(i, 2, "0") * ".dat")
        writedlm(fn, [res.alph res.a real(res.erg1) imag(res.erg1) real(res.erg3) imag(res.erg3) real(res.erg) res.rate_SI tab[i, 2] res.C res.C_low res.phase res.xfrac_SI res.xfrac_SG])
        println("i=$i alpha=$(round(res.alph, digits=4)) l4=$(res.l4s) rate_SI=$(res.rate_SI) table=$(tab[i,2]) C=$(res.C) C_low=$(res.C_low) phase=$(round(res.phase, sigdigits=3)) xfrac SI/SG=$(round(res.xfrac_SI, sigdigits=3))/$(round(res.xfrac_SG, sigdigits=3)) [$(round(time()-t0, digits=1)) s]")
        flush(stdout)
    end
end
