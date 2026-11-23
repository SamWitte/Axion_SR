# Fast version of the full scattering-rate calculation for channels whose fourth quantum is bound (horizon absorption).
# Same ingredients as gf_radial in code/solve_sr_rates.jl - relativistic (Leaver/Heun) cloud modes for n <= rel_nmax from the
# shared cache, hydrogenic ones above, spheroidal angular functions, Klein-Gordon normalisation, contact (SI) and
# self-gravity (SG, with recoil) sources, the same source grid and interpolation, the same rate formula - but the two
# homogeneous radial solutions of the fourth quantum are integrated in double precision with O(1) memory, which makes a
# fine radial step affordable (gf_radial does this recurrence in BigFloat: ~0.5 ms per step).
#   julia spin_fast.jl <channel list> <alphas a,b,..> <spins a,b,..> <levels dir> <cache dir> <out dir> <worker i> <workers N> [key=value ..]
#     rel_nmax=8   levels with n <= rel_nmax relativistic (0: all hydrogenic)
#     hmax=0.5     cap of the step in r*, on top of min(default, 0.12/alpha) and of 0.1/|w4 - m4 Omega_H| for m4 != 0
#     hfac=1       divide the step further (convergence checks)
#     step=fine    "prod": exactly the production step min(default, 0.12/alpha) (for like-for-like tests)
#     kg=1         Klein-Gordon normalisation (0: as before 19 Nov 2026)
#     tag=fast     prefix of the output file
#     wfsub=FILE   list 'n l m spin alpha' of cached relativistic radial functions that are broken (wf_check.jl): for these the
#                  function of the same level at spin asub (default 0.99) is used, with the level energy of the actual spin
#                  (the radial function of a cloud level hardly depends on spin: overlap 0.95 vs 0.99 >= 0.999994)
#     only_sub=1   compute only rows in which a function is substituted (the others are skipped, nothing written)
#     only=FILE    compute only rows that contain a listed 'n l m spin alpha' at one of the requested spins (recomputation
#                  after a correction of level energies: lev_bad.txt, levels_out/lev_zfix.txt)
# Output (out dir/<tag>_<list name>_<i>of<N>.txt), one row per channel and alpha:
#   S1 S2 S3 S4 alpha seconds  then per spin:  status G_SI G_SG G_x     (units of the rate tables, f_a = M_pl; without kH for
#   m4 != 0, as in the tables)
#   status: 1 fourth quantum bound, computed; 2 fourth quantum not bound at this spin (nothing computed);
#           0 an eigenvalue search failed (gf_radial then uses hydrogenic energies; nothing computed here)
#   A row is not written at all (FAILED in the log) while a level energy is still missing from the levels dir.
const CODE = get(ENV, "SPIN_TAB_CODE", normpath(joinpath(@__DIR__, "..", "..")))   # repository src/ (the campaign ran from a frozen copy of it)
include(joinpath(CODE, "solve_sr_rates.jl"))
include(joinpath(CODE, "Core/constants.jl"))
include(joinpath(CODE, "state_utils.jl"))
using Printf

const EPSF = 1e-3
const RPTS = 1500          # source grid, as Compute_sg_points.jl
const NPTS_BND = 2000
const NTOT_SAFE = 30000

struct Src
    lr0::Float64
    dlr::Float64
    T::Matrix{ComplexF64}      # T[k, i]: source k at grid point i
end
@inline function src_at(s::Src, k::Int, r::Float64)
    x = (log10(r) - s.lr0) / s.dlr
    i = clamp(floor(Int, x), 0, size(s.T, 2) - 2); t = x - i
    return (1 - t) * s.T[k, i+1] + t * s.T[k, i+2]
end

# amplitudes amp_k = R_H(r_cut) * Int_{r > 1.01 r+} R_inf T_k dr / W  for the sources k of s
function hom_amp(alph::Float64, a::Float64, w4::Float64, l::Int, m::Int, rmax::Float64, h::Float64, s::Src)
    rp = 1 + sqrt(1 - a^2); rmm = 1 - sqrt(1 - a^2)
    rstar(r) = r + 2rp / (rp - rmm) * log((r - rp) / 2) - 2rmm * log((r - rmm) / 2) / (rp - rmm)
    function rofrs(x, y0)                    # invert r*(r): Newton on y = log(r - rp), start y0
        local yy = y0
        for _ in 1:80
            local rr = rp + exp(yy)
            local dy = (rstar(rr) - x) * (rr^2 - 2rr + a^2) / ((rr^2 + a^2) * exp(yy))
            yy -= clamp(dy, -2.0, 2.0)
            abs(dy) < 1e-13 && break
        end
        return yy
    end
    ycold(x) = x > 10 ? log(max(x - rp, 1e-3)) : (x - rp) * (rp - rmm) / (2rp) + log(2.0) - 1
    gam2 = -a^2 * (w4^2 - alph^2)
    LLM = l * (l + 1) + (-1 + 2l * (l + 1) - 2m^2) * gam2 / (-3 + 4l * (l + 1))
    function VmE(r)
        local delt = r^2 - 2r + a^2; local ff = r^2 + a^2
        delt * alph^2 / ff + delt * (LLM + a^2 * (w4^2 - alph^2)) / ff^2 + delt * (3r^2 - 4r + a^2) / ff^3 - 3delt^2 * r^2 / ff^4 +
            2a * m * w4 / ff - 2a * m * w4 * delt / ff^2 - a^2 * m^2 / ff^2 - w4^2
    end
    kap = sqrt(alph^2 - w4^2)
    rsmax = rstar(rmax)
    while exp(-kap * rsmax) <= 1e-200; rsmax /= 2; end
    rsmin = rstar(rp * (1 + EPSF))
    N = floor(Int, (rsmax - rsmin) / h) + 1          # grid rv_j = rsmax - (j - 1) h, j = 1..N (descending); ascending index i = N + 1 - j
    N < 16 && return (amp=fill(NaN + 0im, size(s.T, 1)), n=N)
    mid = clamp(round(Int, N / 4), 3, N - 3)         # ascending index of the Wronskian point (as in gf_radial)
    ns = size(s.T, 1)
    acc = zeros(ComplexF64, ns); yprev = zeros(ComplexF64, ns)
    Rm = zeros(Float64, 3); rm = zeros(Float64, 3)   # R_inf and r at ascending mid-1, mid, mid+1
    rcut = 1.01 * rp
    # ---- inward pass: decaying solution
    y = rofrs(rsmax, ycold(rsmax)); r = rp + exp(y)
    u0 = exp(-kap * rsmax) * sqrt(r^2 + a^2) / r
    y1 = rofrs(rsmax - h, y); r1 = rp + exp(y1)
    u1 = exp(-kap * (rsmax - h)) * sqrt(r1^2 + a^2) / r1
    rprev = r; R = u0 / sqrt(r^2 + a^2)
    for k in 1:ns; yprev[k] = R * src_at(s, k, r); end
    ih_r = r                                          # smallest r > rcut seen (radius of R_H(r_cut))
    uj, ujm = u1, u0; rj = r1; yj = y1
    for j in 2:N
        # here uj = u_j at radius rj
        R = uj / sqrt(rj^2 + a^2)
        i = N + 1 - j
        if rj > rcut
            for k in 1:ns
                yk = R * src_at(s, k, rj)
                acc[k] += (yk + yprev[k]) * (rprev - rj) / 2
                yprev[k] = yk
            end
            ih_r = rj
        end
        (mid - 1 <= i <= mid + 1) && (Rm[i-mid+2] = R; rm[i-mid+2] = rj)
        j == N && break
        un = 2uj - ujm + h^2 * VmE(rj) * uj
        ujm, uj = uj, un
        rprev = rj
        yj = rofrs(rsmax - j * h, yj - h * (rj^2 - 2rj + a^2) / ((rj^2 + a^2) * (rj - rp)))
        rj = rp + exp(yj)
    end
    # ---- outward pass: ingoing solution at the horizon, up to max(mid + 1, index of r_cut)
    omH = a / (rp^2 + a^2); kH = w4 - m * omH
    x1 = rsmax - (N - 1) * h
    ya = rofrs(x1, ycold(x1)); ra = rp + exp(ya)
    yb = rofrs(x1 + h, ya); rb = rp + exp(yb)
    fa = exp(-im * kH * x1) * sqrt(ra^2 + a^2); fb = exp(-im * kH * (x1 + h)) * sqrt(rb^2 + a^2)
    Rf = zeros(ComplexF64, 3); RH = 0.0 + 0.0im; gotH = false
    i = 2; fi, fim = fb, fa; ri = rb; yi = yb
    (ra > rcut) && (RH = fa / sqrt(ra^2 + a^2); gotH = true)
    while true
        Rv = fi / sqrt(ri^2 + a^2)
        (!gotH && ri > rcut) && (RH = Rv; gotH = true)
        (mid - 1 <= i <= mid + 1) && (Rf[i-mid+2] = Rv)
        (i >= mid + 1 && gotH) && break
        i >= N && break
        fn = 2fi - fim + h^2 * VmE(ri) * fi
        fim, fi = fi, fn
        i += 1
        yi = rofrs(x1 + (i - 1) * h, yi + h * (ri^2 - 2ri + a^2) / ((ri^2 + a^2) * (ri - rp)))
        ri = rp + exp(yi)
    end
    W = (Rf[2] * (Rm[3] - Rm[1]) - Rm[2] * (Rf[3] - Rf[1])) / (rm[3] - rm[1]) * (rm[2]^2 - 2rm[2] + a^2)
    W *= GNew                 # gf_radial: wronk *= Delta * (GNew * M), M = 1
    return (amp=acc .* RH ./ W, n=N)
end

lev(s) = parse_state_from_args(String(s))
lkey(n, l, m, a, al) = @sprintf("%d %d %d %.6g %.6g", n, l, m, a, al)

function load_levels(dir)
    tab = Dict{String,Tuple{Float64,Int}}()
    for f in readdir(dir; join=true)
        startswith(basename(f), "lev_") || continue
        for ln in eachline(f)
            p = split(ln)
            length(p) == 9 || continue
            tab[join(p[1:5], " ")] = (parse(Float64, p[6]), parse(Int, p[8]))
        end
    end
    return tab
end

# one (channel, alpha, spin) point: "status G_SI G_SG G_x"
function fast_point(L, alpp::Float64, a::Float64, tab, P; BHlmax=2)
    (n1, l1, m1), (n2, l2, m2), (n3, l3, m3) = L
    M = 1.0; mu = alpp / (M * GNew); alph = alpp
    eG = [get(tab, lkey(n, l, m, a, alpp), nothing) for (n, l, m) in L]
    any(isnothing, eG) && error("level energy missing at spin $a (row not written; rerun when the levels are complete)")
    any(x -> x[2] == 0 || !isfinite(x[1]), eG) && return "0 0 0 0"
    (eG[1][1] + eG[2][1] - eG[3][1] < alpp) || return "2 0 0 0"
    maxN = max(n1, n2, n3)
    q = (2 * maxN) * (2 * maxN + 1)
    rmax = maxN < 10 ? 2.5 * q / alpp^2 : 10.0 * maxN^2 / alpp^2
    h_def = 0.2 / ((q / alpp^2) / (42 / 0.03^2))
    rpB = BigFloat(1.0 + sqrt(1.0 - a^2)); rp = Float64(rpB)
    rlist = 10 .^ range(log10(rpB * (1.0 + EPSF)), log10(rmax), RPTS)
    rlF = Float64.(rlist)
    rf = Vector{Vector{ComplexF64}}(undef, 3); erg = Vector{Float64}(undef, 3)
    for (i, (n, l, m)) in enumerate(L)
        if i == 2 && L[2] == L[1]
            rf[2] = rf[1]; erg[2] = erg[1]; continue
        end
        done = false
        if n <= P.rel_nmax
            try
                aw = (lkey(n, l, m, a, alpp) in P.wfsub) ? P.asub : a          # broken cache entry: function of spin asub
                rl, r1, e = gf_cached("wf", (n, l, m, mu, M, aw, NPTS_BND, NTOT_SAFE, true)) do
                    solve_radial(mu, M, aw, n, l, m; rpts=NPTS_BND, return_erg=true, Ntot_safe=NTOT_SAFE, use_heunc=true, debug=false, pre_compute_erg=nothing)
                end
                if !(isnan(real(e)) || isnan(imag(e)) || (abs(real(e) - alph) / alph > 1.0) || (abs(imag(e)) > 0.1))
                    itp = LinearInterpolation(log10.(rl), r1, extrapolation_bc=Line())
                    rf[i] = ComplexF64.(itp(log10.(rlist))); erg[i] = aw == a ? Float64(real(e)) : eG[i][1]; done = true
                end
            catch err
                println("solve_radial failed for ", (n, l, m), " a ", a, " alpha ", alpp, ": hydrogenic fallback (as gf_radial)")
            end
        end
        if !done
            rf[i] = ComplexF64.(radial_bound_NR(n, l, m, mu, M, rlist)); erg[i] = eG[i][1]
        end
    end
    w4 = erg[1] + erg[2] - erg[3]
    (abs(w4) < alpp) || return "2 0 0 0"
    m4 = m1 + m2 - m3
    kg_fac = 1.0
    if P.kg
        kg_fac = 1.0 / (kg_norm_ratio(rf[1], rlist, l1, m1, a, alph, erg[1]) * kg_norm_ratio(rf[2], rlist, l2, m2, a, alph, erg[2]) *
                        kg_norm_ratio(rf[3], rlist, l3, m3, a, alph, erg[3]))
    end
    Z1p = sg_theta_profile(spheroidals(l1, m1, a, erg[1], alph)); Z2p = sg_theta_profile(spheroidals(l2, m2, a, erg[2], alph))
    Z3p = sg_theta_profile(spheroidals(l3, m3, a, erg[3], alph))
    ident = L[1] == L[2]
    preFac = ident ? 0.5 : 1.0
    unitMatch = -1.0 / (2 * alpp)^(3 / 2)
    gammaT = preFac * unitMatch .* rf[1] .* rf[2] .* conj.(rf[3])
    tr23 = sg_transition(rf[2], rf[3], Z2p, Z3p, m2, m3, l2, l3, rlF; recoil=true)
    tr13 = ident ? tr23 : sg_transition(rf[1], rf[3], Z1p, Z3p, m1, m3, l1, l3, rlF; recoil=true)
    sg_pref = preFac * unitMatch * 0.5 * 8.0 * alpp^2
    kH = w4 - m4 * a / (2rp)
    h = min(rp / h_def, 0.12 / alpp)
    if P.step != "prod"
        h = min(h, P.hmax)
        m4 != 0 && (h = min(h, 0.1 / abs(kH)))
    end
    h /= P.hfac
    lr0 = log10(rlF[1]); dlr = (log10(rlF[end]) - lr0) / (RPTS - 1)
    l4 = max(abs(l1 + l2 - l3), abs(m4))
    vsi = 0.0; vsg = 0.0; vx = 0.0
    for l in l4:2:max(BHlmax, l4)
        if l != l4
            _, ang_iszero = integral4(l1, m1, l2, m2, l3, -m3, l, -m4)
            ang_iszero && continue
        end
        Z4p = sg_theta_profile(spheroidals(l, m4, a, w4, alph))
        CG = sg_ang([Z1p, Z2p, conj.(Z3p), conj.(Z4p)], [m1, m2, -m3, -m4])
        CG2 = sg_ang([Z1p, Z2p, conj.(Z3p), conj.(Z4p)], [m1, m2, -m3, -m4]; cos2=true)
        T = Matrix{ComplexF64}(undef, 2, RPTS)
        T[1, :] .= gammaT .* (CG .* rlF .^ 2 .+ CG2 * a^2)
        T[2, :] .= sg_pref .* (sg_radial_source(tr23, Z1p, m1, rf[1], Z4p, m4, rlF, a) .+ sg_radial_source(tr13, Z2p, m2, rf[2], Z4p, m4, rlF, a))
        res = hom_amp(alpp, a, w4, l, m4, rmax, h, Src(lr0, dlr, T))
        (isfinite(abs2(res.amp[1])) && isfinite(abs2(res.amp[2]))) || continue
        vsi += abs2(res.amp[1]); vsg += abs2(res.amp[2]); vx += 2 * real(res.amp[1] * conj(res.amp[2]))
    end
    lam = (mu / (M_pl * 1e9))^2
    pref = kg_fac * 4 * alpp * (m4 == 0 ? kH : 1.0) * rp * lam^2 / mu^2 * (GNew * M^2 * M_to_eV)^2
    return @sprintf("1 %.10e %.10e %.10e", pref * vsi, pref * vsg, pref * vx)
end

function main()
    list, levdir, cache, outdir = ARGS[1], ARGS[4], ARGS[5], ARGS[6]
    alphas = parse.(Float64, split(ARGS[2], ",")); spins = parse.(Float64, split(ARGS[3], ","))
    iw, nw = parse(Int, ARGS[7]), parse(Int, ARGS[8])
    kv = Dict(String(split(x, "=")[1]) => String(split(x, "=")[2]) for x in ARGS[9:end])
    P = (rel_nmax=parse(Int, get(kv, "rel_nmax", "8")), hmax=parse(Float64, get(kv, "hmax", "0.5")), hfac=parse(Float64, get(kv, "hfac", "1")),
         step=get(kv, "step", "fine"), kg=get(kv, "kg", "1") == "1",
         wfsub=haskey(kv, "wfsub") ? Set(join(split(l)[1:5], " ") for l in eachline(kv["wfsub"]) if !startswith(l, "#") && length(split(l)) >= 5) : Set{String}(),
         asub=parse(Float64, get(kv, "asub", "0.99")), only_sub=get(kv, "only_sub", "0") == "1",
         only=haskey(kv, "only") ? Set(join(split(l)[1:5], " ") for l in eachline(kv["only"]) if !startswith(l, "#") && length(split(l)) >= 5) : nothing)
    GF_CACHE_DIR[] = cache
    mkpath(outdir)
    tab = load_levels(levdir)
    tag = get(kv, "tag", "fast") * "_" * replace(basename(list), ".txt" => "")
    done = Set{String}()
    for f in readdir(outdir; join=true)
        startswith(basename(f), tag * "_") || continue
        for ln in eachline(f)
            p = split(ln)
            (length(p) == 6 + 4 * length(spins)) && push!(done, join(p[1:5], " "))
        end
    end
    fout = joinpath(outdir, tag * "_$(iw)of$(nw).txt")
    isfile(fout) || open(fout, "w") do f
        println(f, "# spins ", join(spins, " "), "   hmax ", P.hmax, " hfac ", P.hfac, " rel_nmax ", P.rel_nmax, " step ", P.step, " kg ", P.kg,
                isempty(P.wfsub) ? "" : " wfsub $(length(P.wfsub)) entries -> spin $(P.asub)")
    end
    for (k, ln) in enumerate(eachline(list))
        (k - iw) % nw == 0 || continue
        p = split(ln)
        (length(p) < 4 || startswith(p[1], "#")) && continue
        L = (lev(p[1]), lev(p[2]), lev(p[3]))
        for al in alphas
            key = join((p[1:4]..., @sprintf("%.6g", al)), " ")
            key in done && continue
            (P.only_sub && !any(lkey(n, l, m, a, al) in P.wfsub for (n, l, m) in L for a in spins if n <= P.rel_nmax)) && continue
            (P.only !== nothing && !any(lkey(n, l, m, a, al) in P.only for (n, l, m) in L for a in spins)) && continue
            t0 = time()
            try
                r = [fast_point(L, al, a, tab, P) for a in spins]
                open(fout, "a") do f
                    println(f, key, " ", round(time() - t0, digits=2), " ", join(r, " "))
                end
            catch e
                println("FAILED ", key, ": ", first(sprint(showerror, e), 300))
            end
            length(GF_CACHE_MEM) > 400 && empty!(GF_CACHE_MEM)
        end
        flush(stdout)
    end
    println("FAST_DONE ", iw, "/", nw)
end
main()
