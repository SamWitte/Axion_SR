# Superradiance tables of one level from the mixed-precision Leaver solver (Numerics/kerr_modes.jl), in the layout
# that grid_sr_growth.jl writes and Numerics/rate_interpolation.jl reads:
#
#   Imag_zeroC_<nlm>.dat          threshold spin a_c(alpha), where Re w = m Omega_H (two columns: alpha, a_c)
#   Imag_ergC_pos_<nlm>.npz       array (alpha, spin, 3) = [alpha, a, Im w] where Im w > 0, else [alpha, a, 1e-100]
#   Imag_ergC_neg_<nlm>.npz       the same with -Im w where Im w < 0
#
# Im w in units of 1/(G M) (M = 1, mu = alpha). The level is followed continuously in alpha at each spin
# (km_eigenfrequency, gw_kerr_track), so the overtone stays the same at large alpha where the hydrogenic expansion
# no longer identifies it. Im w is then refined to the precision floor of the continued fraction with two settings
# (BigFloat head 256 bits / 600 terms and 512 bits / 1500 terms, Float64 tail of leaver_terms and twice as many);
# their difference is the error sigma. As in the tables of 10 Nov 2026 (commit bac01b2):
#   - Leaver values where resolved: |Im w| >= 300 sigma and >= 1e-31 Re w;
#   - elsewhere (small alpha, and next to the threshold where Im w passes through 0) the perturbative rate
#     (Detweiler, as sr_rates in solve_sr_rates.jl, at the Leaver Re w) times the ratio Leaver/perturbative of the
#     resolved points: interpolated in ln alpha between them, and below the lowest one joined C1 to a ratio that
#     tends to 1 as alpha -> 0 (ln r = ln r_0 (alpha/alpha_0)^p, p from the slope at the join, clamped to [1, 3]);
#   - above the last alpha at which the level could be followed (broad resonances at large alpha below
#     threshold): Im w / (m Omega_H - Re w) held fixed, Re w from the expansion plus the last offset.
#
#   julia -t 1 scripts/build_sr_tables_leaver.jl --level 4-3-3 --outdir DIR [--what zero,erg] [--spins_from NPZ]
#                                                 [--nalpha 60] [--nzero 220] [--spin_index 1:49]
# Writes the three files and a diagnostics table sr_leaver_<nlm>.tsv (one row per grid point) into DIR.
const SRC = joinpath(@__DIR__, "..")
using ArgParse, NPZ, Printf, Suppressor, DelimitedFiles
include(joinpath(SRC, "Core", "constants.jl"))
include(joinpath(SRC, "state_utils.jl"))
@suppress include(joinpath(SRC, "Numerics", "gw_rates.jl"))
include(joinpath(SRC, "Numerics", "gw_kerr_build.jl"))

const A_MAX = 0.998
sr_OmegaH(a) = a / (2 * (1 + sqrt(1 - a^2)))
sr_alpha_max(m) = A_MAX * m / (2 * (1 + sqrt(1 - A_MAX^2))) * 1.3

"""Detweiler's Im w (units 1/(G M)) at Re w = w, exactly as sr_rates in solve_sr_rates.jl (without its clamping)."""
function sr_detweiler(n::Int, l::Int, m::Int, α, a, w)
    rp = 1 + sqrt(1 - a^2)
    Anl = 2.0^(4l + 1) * Float64(factorial(big(l + n)) / (big(n)^(2l + 4) * factorial(big(n - l - 1))))
    Anl *= Float64(factorial(big(l)) / (factorial(big(2l)) * factorial(big(2l + 1))))^2
    χ = prod((k^2 * (1 - a^2) + (a * m - 2rp * α)^2 for k in 1:l); init=1.0)
    return 2α * rp * (m * sr_OmegaH(a) - w) * α^(4l + 4) * Anl * χ
end

"""Im w refined to the precision floor of the continued fraction, from a located root w0; returns (w, sigma)."""
function sr_refine(s::NTuple{3,Int}, α, a, w0)
    n, l, m = s
    gap = km_level_gap(n, l, α)
    Kmax = leaver_terms(α, w0, a)
    # spheroidal eigenvalue as a smooth function of w in extended precision: Λ(c²) = Λ(c0) + (c² - c0) Λ'(c0) around
    # c0 = a²(w0² - α²) (the refinement moves c² by ~1e-16 relative, so the neglected term is ~1e-32). Recomputing
    # Λ from w rounded to Float64 at every step makes the condition jagged at the 1e-17 level and let the secant stop
    # on a jag (isolated Im w wrong by orders of magnitude where Im w ~ 1e-30).
    c0 = a^2 * (ComplexF64(w0)^2 - α^2)
    Λ0 = real(km_spheroidal(l, m, real(c0))[1])
    ε = 1e-12 * max(abs(c0), 1e-6)
    dΛ = imag(km_spheroidal(l, m, real(c0) + im * ε)[1]) / ε       # = -<S|cos²θ|S> (first-order branch: exact)
    res = ComplexF64[]
    for (bits, J, K) in ((256, 600, Kmax), (512, 1500, min(2Kmax, 8_000_000)))
        push!(res, setprecision(BigFloat, bits) do
            Λof(w) = Λ0 + (a^2 * (w^2 - big(α)^2) - real(c0)) * dΛ
            f(w) = km_cf_eval(w, α, a, m, Λof(w), K, min(J, K))[1]
            w1 = Complex{BigFloat}(w0); w2 = w1 + 1e-9 * gap
            f1 = f(w1); f2 = f(w2)
            best = w2; bestf = abs(f2)
            for it in 1:60
                (f2 == f1 || !isfinite(abs(f2))) && break
                w3 = w2 - f2 * (w2 - w1) / (f2 - f1)
                step = abs(w3 - w2)
                w1, f1 = w2, f2
                w2 = w3; f2 = f(w2)
                abs(f2) < bestf && (bestf = abs(f2); best = w2)     # keep the iterate with the smallest residual
                step < 1e-60 * gap && break
            end
            ComplexF64(best)
        end)
    end
    return res[2], abs(imag(res[2]) - imag(res[1]))
end

"""One spin column: the level followed through the alpha grid. Returns per alpha (Re w, Im w Leaver, sigma, ok)."""
function sr_column(s::NTuple{3,Int}, a, αs)
    n, l, m = s
    N = length(αs)
    wr = fill(NaN, N); wi = fill(NaN, N); sg = fill(Inf, N); followed = falses(N)
    hist = nothing
    for (i, α) in enumerate(αs)
        w = nothing
        if hist === nothing
            α <= gw_kerr_alpha_free(l) || break
            e = try km_eigenfrequency(n, l, m, α, a) catch; nothing end
            (e !== nothing && e.converged && e.err < 1e-8 && real(e.ω) < α) || continue
            w = e.ω
            hist = gw_kerr_push(nothing, α, w)
        else
            h = try gw_kerr_track(s, a, hist, α) catch; nothing end
            h === nothing && break                # the level is lost: nothing above this alpha at this spin
            hist = h; w = h[2]
        end
        wref, σ = try sr_refine(s, α, a, w) catch; (w, Inf) end
        wr[i] = real(wref); wi[i] = imag(wref); sg[i] = σ; followed[i] = true
    end
    return wr, wi, sg, followed
end

"""Assemble the column: Leaver where resolved, the perturbative rate times the joined ratio elsewhere, the held
G/(m Omega_H - Re w) above the last followed alpha. Returns (Im w, source) per alpha; source 1 Leaver, 2 joined,
3 held, 0 none."""
function sr_assemble(s::NTuple{3,Int}, a, αs, wr, wi, sg, followed)
    n, l, m = s
    N = length(αs)
    imw = fill(NaN, N); src = zeros(Int, N)
    resolved = [followed[i] && isfinite(wi[i]) && abs(wi[i]) >= 300sg[i] && abs(wi[i]) >= 1e-31 * wr[i] for i in 1:N]
    pert(i) = sr_detweiler(n, l, m, αs[i], a, followed[i] ? wr[i] : real(km_omega_guess(n, l, m, αs[i], a)))
    lr = fill(NaN, N)                            # ln(Leaver/perturbative) at resolved points
    for i in 1:N
        resolved[i] || continue
        p = pert(i)
        if p != 0 && sign(p) == sign(wi[i])
            lr[i] = log(wi[i] / p); imw[i] = wi[i]; src[i] = 1
        else
            resolved[i] = false                  # sign disagreement next to the threshold: treat as unresolved
        end
    end
    R = findall(isfinite, lr)
    lα = log.(αs)
    if !isempty(R)
        k0 = R[1]
        if length(R) >= 2 && lr[k0] != 0 && lr[R[2]] * lr[k0] > 0
            p = clamp((log(abs(lr[R[2]])) - log(abs(lr[k0]))) / (lα[R[2]] - lα[k0]), 1.0, 3.0)
        else
            p = 1.0
        end
        for i in 1:N
            (followed[i] && src[i] == 0) || continue
            if i < k0
                l_i = lr[k0] * exp(p * (lα[i] - lα[k0]))
            else
                j1 = findlast(<(i), R); j2 = findfirst(>(i), R)
                if j1 !== nothing && j2 !== nothing
                    a1, a2 = R[j1], R[j2]
                    l_i = lr[a1] + (lr[a2] - lr[a1]) * (lα[i] - lα[a1]) / (lα[a2] - lα[a1])
                else
                    l_i = lr[R[j1 === nothing ? 1 : j1]]
                end
            end
            imw[i] = pert(i) * exp(l_i); src[i] = 2
        end
    else
        for i in 1:N                              # nothing resolved (all far below 1e-31): perturbative as it is
            followed[i] || continue
            imw[i] = pert(i); src[i] = 2
        end
    end
    last = findlast(followed)
    if last !== nothing && last < N
        g = km_level_gap(n, l, αs[last])
        x = (wr[last] - real(km_omega_guess(n, l, m, αs[last], a))) / g
        A = imw[last] / (m * sr_OmegaH(a) - wr[last])
        for i in (last + 1):N
            w = real(km_omega_guess(n, l, m, αs[i], a)) + x * km_level_gap(n, l, αs[i])
            imw[i] = A * (m * sr_OmegaH(a) - w); src[i] = 3
        end
    end
    return imw, src
end

"""Threshold spin at alpha: bisection on m Omega_H(a) - Re w(alpha, a), the level followed from alpha 0.5 (or alpha)."""
function sr_threshold(s::NTuple{3,Int}, α; it::Int=40)
    n, l, m = s
    function w_at(a)
        α0 = min(α, 0.5, gw_kerr_alpha_free(l))
        e = km_eigenfrequency(n, l, m, α0, a)
        α <= α0 && return real(e.ω)
        h = gw_kerr_track(s, a, gw_kerr_push(nothing, α0, e.ω), α)
        return h === nothing ? NaN : real(h[2])
    end
    f(a) = m * sr_OmegaH(a) - w_at(a)
    lo, hi = 1e-4, A_MAX
    flo = f(lo); fhi = f(hi)
    (isfinite(flo) && isfinite(fhi) && flo < 0 < fhi) || return NaN
    for _ in 1:it
        mid = (lo + hi) / 2; fm = f(mid)
        isfinite(fm) || return NaN
        fm < 0 ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end

function main()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--level";      arg_type = String; required = true
        "--outdir";     arg_type = String; required = true
        "--what";       arg_type = String; default = "zero,erg"
        "--spins_from"; arg_type = String; default = joinpath(SRC, "rate_sve", "Imag_ergC_pos_544.npz")
        "--nalpha";     arg_type = Int; default = 60
        "--nzero";      arg_type = Int; default = 220
        "--spin_index"; arg_type = String; default = ""; help = "range of spin nodes (e.g. 1:10) for a partial run; output gets _part<i>_<j>"
    end
    p = parse_args(s)
    lev = Tuple(parse.(Int, split(p["level"], "-")))
    n, l, m = lev
    nlm = "$(n)$(l)$(m)"
    mkpath(p["outdir"])
    what = split(p["what"], ",")
    if "zero" in what
        αz = exp10.(range(log10(0.03), log10(sr_alpha_max(m)), length=p["nzero"]))
        rows = Tuple{Float64,Float64}[]
        for α in αz
            ac = try sr_threshold(lev, α) catch; NaN end
            isfinite(ac) || continue
            push!(rows, (α, ac))
            @printf("zero %s alpha %.5f a_c %.8f\n", nlm, α, ac); flush(stdout)
        end
        open(joinpath(p["outdir"], "Imag_zeroC_$(nlm).dat"), "w") do io
            for (α, ac) in rows
                @printf(io, "%.18e %.18e\n", α, ac)
            end
        end
    end
    if "erg" in what
        spins = sort(unique(round.(vec(npzread(p["spins_from"])[:, :, 2]), digits=12)))
        js = isempty(p["spin_index"]) ? (1:length(spins)) : eval(Meta.parse(p["spin_index"]))
        αs = exp10.(range(log10(0.03), log10(sr_alpha_max(m)), length=p["nalpha"]))
        P = zeros(length(αs), length(spins), 3); M = zeros(length(αs), length(spins), 3)
        tag = isempty(p["spin_index"]) ? "" : "_part$(first(js))_$(last(js))"
        diag = open(joinpath(p["outdir"], "sr_leaver_$(nlm)$(tag).tsv"), "w")
        println(diag, "alpha\ta\tRe_w\tIm_w_leaver\tsigma\tIm_w_table\tsource")
        for j in js
            a = spins[j]
            t0 = time()
            wr, wi, sg, fol = sr_column(lev, a, αs)
            imw, src = sr_assemble(lev, a, αs, wr, wi, sg, fol)
            for (i, α) in enumerate(αs)
                v = imw[i]
                P[i, j, :] .= (α, a, isfinite(v) && v > 0 ? v : 1e-100)
                M[i, j, :] .= (α, a, isfinite(v) && v < 0 ? -v : 1e-100)
                @printf(diag, "%.12e\t%.10f\t%.15e\t%.6e\t%.3e\t%.6e\t%d\n", α, a, wr[i], wi[i], sg[i], v, src[i])
            end
            flush(diag)
            @printf("erg %s spin %.4f: %d followed, %d Leaver, %d joined, %d held  (%.0f s)\n", nlm, a, count(fol), count(==(1), src),
                    count(==(2), src), count(==(3), src), time() - t0)
            flush(stdout)
        end
        close(diag)
        npzwrite(joinpath(p["outdir"], "Imag_ergC_pos_$(nlm)$(tag).npz"), P)
        npzwrite(joinpath(p["outdir"], "Imag_ergC_neg_$(nlm)$(tag).npz"), M)
    end
end

main()
