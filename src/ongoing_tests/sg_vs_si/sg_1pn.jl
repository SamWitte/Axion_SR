# First post-Newtonian (O(alpha^2)) corrections to the self-gravity 2 -> 2 amplitude  N1 x N1 -> N3 x {BH, Inf}.
#
# One-graviton exchange in harmonic gauge gives (G = 1, units of GM for lengths; see note_1pn/note_1pn.tex)
#     L_eff = G ∫∫ (1/r) [ 1/2 T00 T00' + 1/2 (T00 Tkk' + Tkk T00') - 2 T0i T0i' ]  + retardation  + O(v^4)
# plus the nonlinear BH x cloud (EIH three-body) terms.  For the transition 1 + 1 -> 3 + 4 each piece is written as a
# source s_X(x) for mode 4 (all derivatives moved off psi_4 by parts) and projected with the Green's function of
# sg_vs_si.jl.  We report delta_X = A_X / A_N, the fractional correction to the Newtonian self-gravity amplitude:
#   N   : s_N   = -alpha^2 psi1 Phi[rho]                                   rho = conj(psi3) psi1
#   h0i : s_GM  = -4i A·∇psi1 + 2 alpha Delta psi1 Phi[rho]                 A = Phi[j], j = -(i/2)(conj(psi3)∇psi1 - psi1∇conj(psi3))
#   hij : s_hij = 1/2 s_S   (h00 sourced by stresses gives the other, equal, half)
#         s_S   = -psi1 Phi[∇conj(psi3)·∇psi1] + ∇Phi[rho]·∇psi1 + Phi[rho]∇²psi1 - 6π psi1² conj(psi3)
#   T00 : s_T00 = 1/2 (s_S + 6π psi1² conj(psi3))                           (kinetic correction to T00, gradient form)
#   ret : s_ret = (alpha^2 Delta^2 / 2) psi1 Psi[rho]                        Psi[f] = ∫ f(y) |x - y|
#   BHc : s_BHc = alpha^2 psi1 [Phi[rho]/r + Phi[rho/r] + (1/r)∫rho/r]        (EIH G^2 M mu^2 terms)
# with Delta = (omega1 - omega3) G M.  Hydrogenic eigenfunctions are used for the corrections (they are themselves
# O(alpha^2), so the error is O(alpha^4)); mode 4 uses the Kerr Green's function as before.
#
# Usage:  julia sg_1pn.jl <S1> <S3> <BH|Inf> <i_start> <i_stop>      (alpha grid of the rate_sve table)
#         julia sg_1pn.jl test                                         (integration-by-parts / continuity checks)

include(joinpath(@__DIR__, "sg_vs_si.jl"))   # gf_amplitudes, gauss_legendre, ang_integral, SRC, solve_sr_rates.jl

const NTH = 40
const THx, THw = gauss_legendre(NTH)
const TH = acos.(THx)
const SINTH = sqrt.(1 .- THx .^ 2)

# ---------------------------------------------------------------- angular building blocks
Yvec(l, m) = [l >= abs(m) ? sphericalY(l, m, t, 0.0) : 0.0im for t in TH]
function dYvec(l, m; h=1e-6)
    l >= abs(m) || return zeros(ComplexF64, NTH)
    [(sphericalY(l, m, t + h, 0.0) - sphericalY(l, m, t - h, 0.0)) / (2h) for t in TH]
end

# a field on the (r, theta) grid with azimuthal number m:  F(r,θ) e^{i m φ}
struct Fld
    v::Matrix{ComplexF64}
    m::Int
end
Base.:*(a::Fld, b::Fld) = Fld(a.v .* b.v, a.m + b.m)
Base.:+(a::Fld, b::Fld) = (@assert a.m == b.m; Fld(a.v .+ b.v, a.m))
Base.:-(a::Fld, b::Fld) = (@assert a.m == b.m; Fld(a.v .- b.v, a.m))
Base.:*(c::Number, a::Fld) = Fld(c .* a.v, a.m)
Base.conj(a::Fld) = Fld(conj.(a.v), -a.m)
rmul(a::Fld, f::AbstractVector) = Fld(a.v .* f, a.m)          # multiply by a function of r
tmul(a::Fld, g::AbstractVector) = Fld(a.v .* g', a.m)         # multiply by a function of θ

# a state: psi = R(r) Y_lm(θ) e^{imφ};  gradient in spherical components and Laplacian
struct State
    psi::Fld; gr::Fld; gth::Fld; gph::Fld; lap::Fld
end
function make_state(R, dR, d2R, r, l, m)
    Y = Yvec(l, m); dY = dYvec(l, m)
    psi = Fld(R * transpose(Y), m)
    gr = Fld(dR * transpose(Y), m)
    gth = Fld((R ./ r) * transpose(dY), m)
    gph = Fld((im * m) .* (R ./ r) * transpose(Y ./ SINTH), m)
    lap = Fld((d2R .+ 2 .* dR ./ r .- l * (l + 1) .* R ./ r .^ 2) * transpose(Y), m)
    return State(psi, gr, gth, gph, lap)
end
cstate(s::State) = State(conj(s.psi), conj(s.gr), conj(s.gth), conj(s.gph), conj(s.lap))
graddot(a::State, b::State) = a.gr * b.gr + a.gth * b.gth + a.gph * b.gph

# Cartesian spherical components of a vector given by spherical components (vr, vth, vph), all with azimuthal M:
#   V+ = Vx + iVy  (azimuthal M+1),  V- = Vx - iVy (M-1),  Vz (M)
function cart(vr::Fld, vth::Fld, vph::Fld)
    s = SINTH; c = THx
    vp = Fld(vr.v .* s' .+ vth.v .* c' .+ im .* vph.v, vr.m + 1)
    vm = Fld(vr.v .* s' .+ vth.v .* c' .- im .* vph.v, vr.m - 1)
    vz = Fld(vr.v .* c' .- vth.v .* s', vr.m)
    return vp, vm, vz
end

# ---------------------------------------------------------------- radial kernels (multipole by multipole)
cumtrap(f, r) = (c = zeros(ComplexF64, length(r)); for i in 2:length(r); c[i] = c[i-1] + 0.5 * (f[i] + f[i-1]) * (r[i] - r[i-1]); end; c)
revtrap(f, r) = (c = zeros(ComplexF64, length(r)); for i in length(r)-1:-1:1; c[i] = c[i+1] + 0.5 * (f[i] + f[i+1]) * (r[i+1] - r[i]); end; c)
# Q(r) = r^A ∫_0^r r'^B f + r^C ∫_r^∞ r'^D f,  and dQ/dr (boundary terms cancel when A+B = C+D)
function kern(f, r, A, B, C, D)
    ci = cumtrap(f .* r .^ B, r); co = revtrap(f .* r .^ D, r)
    return r .^ A .* ci .+ r .^ C .* co, A .* r .^ (A - 1) .* ci .+ C .* r .^ (C - 1) .* co
end

function multipoles(F::Fld, Lmax)
    M = F.m
    Ls = abs(M):Lmax
    coeffs = Dict(L => (2pi) .* (F.v * (THw .* conj.(Yvec(L, M)))) for L in Ls)
    return Ls, coeffs
end

# Coulomb potential  Phi[F](x) = ∫ F(y)/|x-y| d^3y   (returns Phi and dPhi/dr, dPhi/dθ as fields)
function coulomb(F::Fld, r, Lmax)
    Ls, fL = multipoles(F, Lmax)
    P = zeros(ComplexF64, length(r), NTH); Pr = zeros(ComplexF64, length(r), NTH); Pth = zeros(ComplexF64, length(r), NTH)
    for L in Ls
        Q, dQ = kern(fL[L], r, -(L + 1), L + 2, L, 1 - L)
        c = 4pi / (2L + 1)
        Y = Yvec(L, F.m); dY = dYvec(L, F.m)
        P .+= c .* Q * transpose(Y); Pr .+= c .* dQ * transpose(Y); Pth .+= c .* (Q ./ r) * transpose(dY)
    end
    return Fld(P, F.m), Fld(Pr, F.m), Fld(Pth, F.m)
end

# Psi[F](x) = ∫ F(y) |x-y| d^3y,  |x-y| = sum_L P_L [ r<^{L+2}/((2L+3) r>^{L+1}) - r<^L/((2L-1) r>^{L-1}) ]
function linear_kernel(F::Fld, r, Lmax)
    Ls, fL = multipoles(F, Lmax)
    P = zeros(ComplexF64, length(r), NTH)
    for L in Ls
        Q1, _ = kern(fL[L], r, -(L + 1), L + 4, L + 2, 1 - L)
        Q2, _ = kern(fL[L], r, -(L - 1), L + 2, L, 3 - L)
        P .+= (4pi / (2L + 1)) .* (Q1 ./ (2L + 3) .- Q2 ./ (2L - 1)) * transpose(Yvec(L, F.m))
    end
    return Fld(P, F.m)
end

# projection onto conj(Y_{l4 m4}) with the Kerr Sigma weight  -> radial source T(r)
project(S::Fld, r, a, l4, m4) = (@assert S.m == m4; (2pi) .* ((S.v .* (r .^ 2 .+ a^2 .* (THx .^ 2)')) * (THw .* conj.(Yvec(l4, m4)))))
# full volume integral ∫ F d^3x (flat measure), F must have azimuthal number 0
volint(F::Fld, r) = (@assert F.m == 0; 2pi * real(trapz(r .^ 2 .* (F.v * THw), r)))

# ---------------------------------------------------------------- hydrogenic states with derivatives on a log grid
function hyd_state(n, l, m, alph, r)
    mu = alph / GNew
    R = complex.(Float64.(radial_bound_NR(n, l, m, mu, 1.0, r)))
    u = log.(r); h = u[2] - u[1]
    d(f) = (g = similar(f); g[3:end-2] = (-f[5:end] .+ 8f[4:end-1] .- 8f[2:end-3] .+ f[1:end-4]) ./ (12h);
            g[1:2] = (f[2:3] .- f[1:2]) ./ h; g[end-1:end] = (f[end-1:end] .- f[end-2:end-1]) ./ h; g)
    Ru = d(R); Ruu = d(Ru)
    dR = Ru ./ r; d2R = (Ruu .- Ru) ./ r .^ 2
    return make_state(R, dR, d2R, r, l, m)
end
bohr(n, alph) = alph - alph^3 / (2n^2)          # omega G M at Bohr order

# ---------------------------------------------------------------- sources
function all_sources(s1::State, s3::State, alph, Delta, r, Lmax)
    c3 = cstate(s3)
    rho = c3.psi * s1.psi                                          # conj(psi3) psi1
    Phi, Phir, Phith = coulomb(rho, r, Lmax)
    # Newtonian
    sN = (-alph^2) * (s1.psi * Phi)
    # gravitomagnetic: j = -(i/2)(conj(psi3) ∇psi1 - psi1 ∇conj(psi3))
    jr = (-0.5im) * (c3.psi * s1.gr - s1.psi * c3.gr)
    jth = (-0.5im) * (c3.psi * s1.gth - s1.psi * c3.gth)
    jph = (-0.5im) * (c3.psi * s1.gph - s1.psi * c3.gph)
    jp, jm, jz = cart(jr, jth, jph)
    Ap, _, _ = coulomb(jp, r, Lmax + 1); Am, _, _ = coulomb(jm, r, Lmax + 1); Az, _, _ = coulomb(jz, r, Lmax + 1)
    gp, gm, gz = cart(s1.gr, s1.gth, s1.gph)                      # ∇+ psi1, ∇- psi1, ∇z psi1
    Adotg = 0.5 * (Ap * gm) + 0.5 * (Am * gp) + Az * gz
    sGM = (-4im) * Adotg + (2 * alph * Delta) * (s1.psi * Phi)
    # stresses (h_ij + stress-sourced h_00)
    gdg = graddot(c3, s1)
    PhiK, _, _ = coulomb(gdg, r, Lmax + 2)
    gradPhi_dot_g1 = Phir * s1.gr + Phith * s1.gth + Fld((im * Phi.m) .* (Phi.v ./ r) ./ SINTH', Phi.m) * s1.gph
    contact = s1.psi * s1.psi * c3.psi
    sS = (-1.0) * (s1.psi * PhiK) + gradPhi_dot_g1 + Phi * s1.lap - (6pi) * contact
    sT00 = 0.5 * (sS + (6pi) * contact)
    # retardation
    Psi = linear_kernel(rho, r, Lmax)
    sret = (alph^2 * Delta^2 / 2) * (s1.psi * Psi)
    # EIH three-body (BH x cloud)
    Phi_rho_r, _, _ = coulomb(rmul(rho, 1 ./ r), r, Lmax)
    s3b = alph^2 * (rmul(s1.psi * Phi, 1 ./ r) + s1.psi * Phi_rho_r)
    if rho.m == 0   # separable 1/(r_x r_y) piece survives only if rho has a monopole (never for our channels)
        s3b = s3b + alph^2 * rmul(volint(rmul(rho, 1 ./ r), r) * s1.psi, 1 ./ r)
    end
    return Dict(:N => sN, :GM => sGM, :hij => 0.5 * sS, :S => sS, :T00 => sT00, :ret => sret, :BHc => s3b),
           (; rho, Phi, jp, jm, jz, Ap, Am, Az, c3, gdg)
end

# ---------------------------------------------------------------- one alpha point
function point_1pn(st1, st3, to_inf, alph; a=0.95, rpts=3000)
    n1, l1, m1 = st1; n3, l3, m3 = st3
    mu = alph / GNew; M = 1.0
    m4 = 2m1 - m3; maxN = max(n1, n3)
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
    src, _ = all_sources(s1, s3, alph, Delta, r, Lmax)

    lmin = max(abs(2l1 - l3), abs(m4)); iseven(lmin + 2l1 + l3) || (lmin += 1)
    keys_ = [:N, :GM, :hij, :T00, :ret, :BHc]
    T = [project(src[k], r, a, lmin, m4) for k in keys_]
    amps, _, _ = gf_amplitudes(erg, a, alph, lmin, m4, r, T, rmax, to_inf; h_mve=h_mve, eps_fac=eps_fac)
    return Dict(k => amps[i] / amps[1] for (i, k) in enumerate(keys_))
end

# ---------------------------------------------------------------- self-tests (all-bound, hydrogenic)
# Replace mode 4 by a bound hydrogenic state and compare  ∫ conj(psi4) s_X d^3x  with the symmetric double integral
# evaluated directly (derivatives acting on psi4), i.e. check every integration by parts and the continuity step.
function selftest(; alph=0.2)
    r = collect(10 .^ range(log10(1e-3), log10(150 / alph^2), 4000))
    s1 = hyd_state(2, 1, 1, alph, r); s3 = hyd_state(3, 2, 2, alph, r); s4 = hyd_state(3, 2, 0, alph, r)
    Delta = bohr(2, alph) - bohr(3, alph); Lmax = 1 + 2 + 3
    src, aux = all_sources(s1, s3, alph, Delta, r, Lmax)
    c4 = cstate(s4)
    via_source(k) = volint(c4.psi * src[k], r)
    # direct evaluations
    rho41 = c4.psi * s1.psi
    N_dir = -alph^2 * volint(rho41 * aux.Phi, r)
    jr = (-0.5im) * (c4.psi * s1.gr - s1.psi * c4.gr); jth = (-0.5im) * (c4.psi * s1.gth - s1.psi * c4.gth); jph = (-0.5im) * (c4.psi * s1.gph - s1.psi * c4.gph)
    jp4, jm4, jz4 = cart(jr, jth, jph)
    GM_dir = 4 * volint(0.5 * (jp4 * aux.Am) + 0.5 * (jm4 * aux.Ap) + jz4 * aux.Az, r)
    # T_kk^{41} = ∇conj(psi4)·∇psi1 - 3/4 ∇²(conj(psi4) psi1),  ∇²(fg) = f∇²g + g∇²f + 2∇f·∇g
    g41 = graddot(c4, s1)
    lap41 = c4.psi * s1.lap + s1.psi * c4.lap + 2.0 * g41
    Tkk41 = g41 - 0.75 * lap41
    Tkk31_pot, _, _ = coulomb(aux.gdg - 0.75 * (aux.c3.psi * s1.lap + s1.psi * aux.c3.lap + 2.0 * aux.gdg), r, Lmax + 2)
    S_dir = -(volint(rho41 * Tkk31_pot, r) + volint(Tkk41 * aux.Phi, r))
    # continuity check: ∇·j31 = i alpha Delta rho31  ->  compare ∫ conj(psi4) psi1 ∇·A  two ways is implicit in GM_dir
    println("selftest alpha=$alph (mode 4 -> bound 320):")
    for (k, dir) in [(:N, N_dir), (:GM, GM_dir), (:S, S_dir)]
        v = via_source(k)
        println("  $k  via source = $(round(v, sigdigits=8))   direct = $(round(dir, sigdigits=8))   rel.diff = $(round(abs(v - dir) / abs(dir), sigdigits=3))")
    end
    lapchk = maximum(abs.(s1.lap.v .- ((alph^4 / 4) .- 2 * alph^2 ./ r) .* s1.psi.v)[r .> 1, :]) / maximum(abs.(s1.lap.v))
    println("  hydrogenic ∇²psi check (Schrödinger): max rel. deviation = $(round(lapchk, sigdigits=3))")
end

# ---------------------------------------------------------------- driver
if abspath(PROGRAM_FILE) == @__FILE__
    if ARGS[1] == "test"
        selftest(alph=0.2); selftest(alph=0.05)
    else
        S1, S3, dest = ARGS[1], ARGS[2], ARGS[3]
        st(s) = (parse(Int, s[1]), parse(Int, s[2]), parse(Int, s[3]))
        tab = readdlm(joinpath(SRC, "rate_sve", "$(S1)_$(S1)_$(S3)_$(dest)_LvrHc_.dat"))
        outdir = joinpath(@__DIR__, "out_1pn"); mkpath(outdir)
        for i in parse(Int, ARGS[4]):min(parse(Int, ARGS[5]), size(tab, 1))
            t0 = time()
            d = point_1pn(st(S1), st(S3), dest == "Inf", tab[i, 1])
            row = [tab[i, 1]; vcat([[real(d[k]), imag(d[k])] for k in (:GM, :hij, :T00, :ret, :BHc)]...)]
            writedlm(joinpath(outdir, "$(S1)_$(S1)_$(S3)_$(dest)_" * lpad(i, 2, "0") * ".dat"), permutedims(row))
            println("i=$i alpha=$(round(tab[i,1], digits=4))  GM=$(round(d[:GM], sigdigits=3))  hij=$(round(d[:hij], sigdigits=3))  T00=$(round(d[:T00], sigdigits=3))  ret=$(round(d[:ret], sigdigits=3))  3b=$(round(d[:BHc], sigdigits=3))  [$(round(time()-t0, digits=1)) s]")
            flush(stdout)
        end
    end
end
