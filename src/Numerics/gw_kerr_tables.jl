"""
    gw_kerr_tables.jl

Relativistic GW rate tables (gw_model = :rel; rate_sve/gw_kerr/{tr,ann}_<a>_<b>.dat, columns GW_KERR_COLS, built by
gw_kerr_build.jl): ratio = Γ_kerr/Γ_nr per channel, Re ω/μ per level, bilinear in (ln α, a), clamped outside the grid.
"""

using LinearAlgebra: eigen, Symmetric

const GW_KERR_TABLE_DIR = joinpath(@__DIR__, "..", "rate_sve", "gw_kerr")
gw_kerr_table_dir() = get(ENV, "NL_GW_KERR_DIR", GW_KERR_TABLE_DIR)

const GW_KERR_COLS = ("alpha", "a", "Gamma_kerr", "Gamma_nr", "ratio", "omega_tilde", "Re_omega_a", "Re_omega_b",
                      "Im_omega_a", "Im_omega_b", "Gamma_H", "nl", "secs", "status", "chk", "lfrac")

gw_kerr_level_str(s) = "$(s[1])-$(s[2])-$(s[3])"

"""Scale of the asinh mapping of the horizon fraction Γ_H/Γ_∞ in the interpolants (gw_hfrac inverts it)."""
const GW_HFRAC_SCALE = 1e-6
gw_kerr_table_name(kind::Symbol, sa, sb) = "$(kind)_$(gw_kerr_level_str(sa))_$(gw_kerr_level_str(sb)).dat"
gw_kerr_parse_level(s::AbstractString) = (p = parse.(Int, split(s, "-")); (p[1], p[2], p[3]))

# ----------------------------------------------------------------------------
# Bilinear interpolation in (ln α, a) with edge clamping
# ----------------------------------------------------------------------------

struct GWInterp2D
    x::Vector{Float64}        # ln α, increasing
    y::Vector{Float64}        # a, increasing
    z::Matrix{Float64}        # size (length(x), length(y)), finite everywhere
end

# index j and weight s with v[j] <= t <= v[j+1]; clamps outside; NaN/∞ go to the lower edge
@inline function _gw_bracket(v::Vector{Float64}, t::Float64)
    n = length(v)
    (n == 1 || !(t > v[1])) && return 1, 0.0
    t >= v[n] && return n - 1, 1.0
    j = clamp(searchsortedlast(v, t), 1, n - 1)
    return j, (t - v[j]) / (v[j + 1] - v[j])
end

"""Value of the table at (α, a); constant outside the grid."""
@inline function gw_interp2(it::GWInterp2D, α::Real, a::Real)
    lx = α > 0 ? log(Float64(α)) : -Inf
    i, s = _gw_bracket(it.x, lx)
    j, t = _gw_bracket(it.y, Float64(a))
    z = it.z
    i1 = min(i + 1, length(it.x)); j1 = min(j + 1, length(it.y))
    @inbounds return (1 - s) * ((1 - t) * z[i, j] + t * z[i, j1]) + s * ((1 - t) * z[i1, j] + t * z[i1, j1])
end

gw_interp_alpha_range(it::GWInterp2D) = (exp(it.x[1]), exp(it.x[end]))

"""Noisy rows (direct and kernel forms of the source disagree, gw_kerr_table_point) whose disagreement |chk| is at
most this are used: the rate is then good to about that fraction."""
const GW_KERR_NOISY_ACCEPT = 0.05
"""Interior gaps along α of at most this many grid nodes are interpolated (gw_kerr_fill_grid!)."""
const GW_KERR_MAX_GAP = 2
"""At the top of a table's α range at most this many nodes are extrapolated."""
const GW_KERR_MAX_EXTRAP = 1
"""Largest change of ln(ratio) an extrapolation over one node may make (0.3 dex)."""
const GW_KERR_MAX_EXTRAP_STEP = 0.3 * log(10)

"""How a grid value of a ratio interpolant was obtained (gw_kerr_fill_grid!)."""
const GW_FILL_NAMES = Dict(Int8(0) => "computed", Int8(1) => "interpolated", Int8(2) => "limit", Int8(3) => "extrapolated",
                           Int8(4) => "spin-interpolated", Int8(5) => "held")

"""Does a table row carry a ratio the interpolants use? ok; broad (computed, with a level wider than half its
spacing: such a level is emptied by the hole within a few beat periods, so its rates matter little, and the
computed value is the best estimate there is); noisy with |chk| <= GW_KERR_NOISY_ACCEPT."""
gw_kerr_usable(r) = isfinite(r.ratio) && r.ratio > 0 &&
    (r.status == "ok" || r.status == "broad" || (r.status == "noisy" && isfinite(r.chk) && abs(r.chk) <= GW_KERR_NOISY_ACCEPT))

"""+1 if level s = (n, l, m) is superradiant at (α, a) (m Ω_H > ω, hydrogenic ω), -1 otherwise."""
function gw_kerr_sr_sign(s::NTuple{3,Int}, α, a)
    ΩH = a / (2 * (1 + sqrt(1 - a^2)))
    return s[3] * ΩH > α * gw_omega(s..., α, a) ? 1 : -1
end

"""
    gw_kerr_fill_grid!(z, αs, as, kind, sa, sb) -> how::Matrix{Int8}
Fill NaNs of z = ln(ratio) (rows α, columns a) by steps 1-5, never across a superradiance threshold (GW_FILL_NAMES).
"""
function gw_kerr_fill_grid!(z::Matrix{Float64}, αs::AbstractVector, as::AbstractVector, kind::Symbol, sa, sb)
    nα, na = size(z)
    how = fill(Int8(-1), nα, na)
    any(isfinite, z) || return how
    how[isfinite.(z)] .= 0
    lα = log.(αs)
    sr = [(gw_kerr_sr_sign(sa, αs[i], as[j]), gw_kerr_sr_sign(sb, αs[i], as[j])) for i in 1:nα, j in 1:na]
    for j in 1:na
        v = findall(i -> how[i, j] == 0, 1:nα)
        isempty(v) && continue
        for (k1, k2) in zip(v[1:end-1], v[2:end])                   # 1. interior gaps
            L = k2 - k1 - 1
            (L == 0 || L > GW_KERR_MAX_GAP) && continue
            all(i -> sr[i, j] == sr[k1, j], k1:k2) || continue
            for i in (k1 + 1):(k2 - 1)
                t = (lα[i] - lα[k1]) / (lα[k2] - lα[k1])
                z[i, j] = (1 - t) * z[k1, j] + t * z[k2, j]; how[i, j] = 1
            end
        end
        k = v[1]                                                     # 2. below the lowest computed α
        if k > 1
            flat = length(v) >= 2 && abs(z[v[2], j] - z[k, j]) <= 0.05 * (lα[v[2]] - lα[k])
            pfit = length(v) >= 2 && z[k, j] * z[v[2], j] > 0 ? log(z[v[2], j] / z[k, j]) / (lα[v[2]] - lα[k]) : 1.0
            if kind == :tr && (pfit >= 0.5 || !flat)
                # approaching ratio 1 (or too little information to say otherwise): ln ratio ∝ α^p
                p = clamp(pfit, 1.0, 2.0)
                for i in 1:(k - 1)
                    z[i, j] = z[k, j] * exp(p * (lα[i] - lα[k])); how[i, j] = pfit >= 0.5 ? 2 : 5
                end
            elseif kind == :tr
                # a transition whose ratio levelled off away from 1: hold it, as for annihilations
                for i in 1:(k - 1)
                    z[i, j] = z[k, j]; how[i, j] = 2
                end
            else
                for i in 1:(k - 1)
                    z[i, j] = z[k, j]; how[i, j] = flat ? 2 : 5
                end
            end
        end
        k = v[end]                                                   # 3. above the highest computed α
        if k < nα && length(v) >= 2
            k0 = v[end - 1]
            slope = (z[k, j] - z[k0, j]) / (lα[k] - lα[k0])
            for i in (k + 1):min(nα, k + GW_KERR_MAX_EXTRAP)
                step = slope * (lα[i] - lα[k])
                (abs(step) <= GW_KERR_MAX_EXTRAP_STEP && all(q -> sr[q, j] == sr[k, j], (k0, i))) || break
                z[i, j] = z[k, j] + step; how[i, j] = 3
            end
        end
    end
    for i in 1:nα, j in 1:na                                         # 4. along spin
        how[i, j] == -1 || continue
        jl = findlast(q -> 0 <= how[i, q] <= 3, 1:(j - 1))
        jr = findfirst(q -> 0 <= how[i, q] <= 3, (j + 1):na)
        (jl === nothing || jr === nothing) && continue
        jr += j
        all(q -> sr[i, q] == sr[i, j], (jl, jr)) || continue
        t = (as[j] - as[jl]) / (as[jr] - as[jl])
        z[i, j] = (1 - t) * z[i, jl] + t * z[i, jr]; how[i, j] = 4
    end
    if any(==(Int8(-1)), how)                                        # 5. held
        miss = how .== -1
        gw_fill_nan!(z)
        how[miss] .= 5
    end
    return how
end

"""
    gw_kerr_ratio_interp(kind, sa, sb, rows) -> (GWInterp2D or nothing, αs, as, how)
Interpolant of ln(ratio) on the full grid of the usable rows, gaps filled by gw_kerr_fill_grid!.
"""
function gw_kerr_ratio_interp(kind::Symbol, sa, sb, rows)
    αs = sort(unique(round.([r.alpha for r in rows], sigdigits=10)))
    as = sort(unique(round.([r.a for r in rows], sigdigits=10)))
    z = fill(NaN, length(αs), length(as))
    for r in rows
        gw_kerr_usable(r) || continue
        z[searchsortedfirst(αs, round(r.alpha, sigdigits=10)), searchsortedfirst(as, round(r.a, sigdigits=10))] = log(r.ratio)
    end
    how = gw_kerr_fill_grid!(z, αs, as, kind, sa, sb)
    any(isfinite, z) || return nothing, αs, as, how
    return GWInterp2D(log.(αs), as, z), αs, as, how
end

"""Fill NaNs of z (rows α, columns a) with the nearest finite value, along α first (lower α preferred), then whole
columns; false if nothing is finite."""
function gw_fill_nan!(z::Matrix{Float64})
    any(isfinite, z) || return false
    nx, ny = size(z)
    for j in 1:ny
        col = view(z, :, j)
        any(isfinite, col) || continue
        last = NaN
        for i in 1:nx                     # forward (from low α)
            isfinite(col[i]) ? (last = col[i]) : (isfinite(last) && (col[i] = last))
        end
        last = NaN
        for i in nx:-1:1                  # backward for leading NaNs
            isfinite(col[i]) ? (last = col[i]) : (col[i] = last)
        end
    end
    good = [all(isfinite, view(z, :, j)) for j in 1:ny]
    for j in 1:ny
        good[j] && continue
        k = argmin([good[q] ? abs(q - j) : typemax(Int) for q in 1:ny])
        z[:, j] .= z[:, k]
    end
    return true
end

"""Interpolant from scattered (α, a, value) points on a rectangular grid (missing points -> filled)."""
function gw_interp_from_points(pts::Vector{NTuple{3,Float64}})
    isempty(pts) && return nothing
    αs = sort(unique(round.(first.(pts), sigdigits=10)))
    as = sort(unique(round.(getindex.(pts, 2), sigdigits=10)))
    z = fill(NaN, length(αs), length(as))
    for (α, a, v) in pts
        i = searchsortedfirst(αs, round(α, sigdigits=10)); j = searchsortedfirst(as, round(a, sigdigits=10))
        isfinite(v) && (z[i, j] = v)
    end
    gw_fill_nan!(z) || return nothing
    return GWInterp2D(log.(αs), as, z)
end

# ----------------------------------------------------------------------------
# Reading the channel files
# ----------------------------------------------------------------------------

"""
    gw_kerr_read_table(path) -> (kind, sa, sb, rows)
rows: NamedTuples with the columns of GW_KERR_COLS (numbers parsed, status a String).
"""
function gw_kerr_read_table(path::AbstractString)
    kind = :none; sa = sb = (0, 0, 0)
    rows = NamedTuple[]
    for line in eachline(path)
        if startswith(line, "#")
            m = match(r"kind=(\w+)\s+a=([\d-]+)\s+b=([\d-]+)", line)
            m === nothing || (kind = Symbol(m[1]); sa = gw_kerr_parse_level(m[2]); sb = gw_kerr_parse_level(m[3]))
            continue
        end
        f = split(line)
        length(f) >= 12 || continue
        vals = [tryparse(Float64, f[k]) for k in 1:11]
        any(x -> x === nothing, vals) && continue
        nl = tryparse(Int, f[12]); nl === nothing && (nl = 0)
        secs = length(f) >= 13 ? something(tryparse(Float64, f[13]), NaN) : NaN
        status = length(f) >= 14 ? String(f[14]) : "ok"
        chk = length(f) >= 15 ? something(tryparse(Float64, f[15]), NaN) : NaN
        lfrac = length(f) >= 16 ? String(f[16]) : "-"
        push!(rows, (alpha=vals[1], a=vals[2], Gamma_kerr=vals[3], Gamma_nr=vals[4], ratio=vals[5],
                     omega_tilde=vals[6], Re_omega_a=vals[7], Re_omega_b=vals[8], Im_omega_a=vals[9],
                     Im_omega_b=vals[10], Gamma_H=vals[11], nl=nl, secs=secs, status=status, chk=chk, lfrac=lfrac))
    end
    kind == :none && error("gw_kerr_read_table: no '# kind= a= b=' header in $(path)")
    return kind, sa, sb, rows
end

"""
    GWKerrTables
Interpolants: ratio_* = ln(Γ_kerr/Γ_nr), omega = Re ω/μ - gw_omega, hfrac, lfrac per ℓ̃; fill codes as GW_FILL_NAMES.
"""
struct GWKerrTables
    ratio_ann::Dict{Tuple{NTuple{3,Int},NTuple{3,Int}}, GWInterp2D}
    ratio_tr::Dict{Tuple{NTuple{3,Int},NTuple{3,Int}}, GWInterp2D}
    omega::Dict{NTuple{3,Int}, GWInterp2D}
    hfrac::Dict{Tuple{Symbol,NTuple{3,Int},NTuple{3,Int}}, GWInterp2D}
    lfrac::Dict{Tuple{Symbol,NTuple{3,Int},NTuple{3,Int}}, Dict{Int,GWInterp2D}}
    files::Vector{String}
    fill::Dict{Tuple{Symbol,NTuple{3,Int},NTuple{3,Int}}, Tuple{Vector{Float64},Vector{Float64},Matrix{Int8}}}
end

"""Parse an lfrac column entry ("2:0.98;3:0.02", or "-") into ℓ̃ => fraction."""
function gw_kerr_parse_lfrac(s::AbstractString)
    d = Dict{Int,Float64}()
    (isempty(s) || s == "-") && return d
    for p in split(s, ';')
        kv = split(p, ':')
        length(kv) == 2 || continue
        l = tryparse(Int, kv[1]); f = tryparse(Float64, kv[2])
        (l === nothing || f === nothing) || (d[l] = f)
    end
    return d
end

gw_pair_key(sa, sb) = sa <= sb ? (sa, sb) : (sb, sa)

"""Level pairs (gw_pair_key) whose tabulated horizon fraction of the transition reaches h < -1 somewhere: the
transitions the hole can reverse (same-n pairs, and pairs with an absorbed level at large α)."""
function gw_kerr_reversible_pairs(tabs)
    zc = asinh(-1 / GW_HFRAC_SCALE)
    return Set((k[2], k[3]) for (k, it) in tabs.hfrac if k[1] == :tr && minimum(it.z) < zc)
end

"""
    gw_kerr_load_tables(dir=gw_kerr_table_dir(); verbose=false) -> GWKerrTables or nothing
Read every tr_/ann_ table in `dir` (cached per directory); nothing if none is usable.
"""
function gw_kerr_load_tables(dir=gw_kerr_table_dir(); verbose::Bool=false)
    key = "kerr:" * abspath(dir)
    haskey(GW_TABLE_CACHE, key) && return GW_TABLE_CACHE[key]
    isdir(dir) || return nothing
    files = sort(filter(f -> endswith(f, ".dat") && (startswith(f, "tr_") || startswith(f, "ann_")), readdir(dir)))
    isempty(files) && return nothing
    rann = Dict{Tuple{NTuple{3,Int},NTuple{3,Int}}, GWInterp2D}()
    rtr = Dict{Tuple{NTuple{3,Int},NTuple{3,Int}}, GWInterp2D}()
    hfr = Dict{Tuple{Symbol,NTuple{3,Int},NTuple{3,Int}}, GWInterp2D}()
    lfr = Dict{Tuple{Symbol,NTuple{3,Int},NTuple{3,Int}}, Dict{Int,GWInterp2D}}()
    ωpts = Dict{NTuple{3,Int}, Dict{Tuple{Float64,Float64}, Float64}}()
    fil = Dict{Tuple{Symbol,NTuple{3,Int},NTuple{3,Int}}, Tuple{Vector{Float64},Vector{Float64},Matrix{Int8}}}()
    used = String[]
    for f in files
        kind, sa, sb, rows = gw_kerr_read_table(joinpath(dir, f))
        ok = [r for r in rows if gw_kerr_usable(r)]
        isempty(ok) && (verbose && println("  gw_kerr: $(f) has no valid row; skipped"); continue)
        it, αg, ag, how = gw_kerr_ratio_interp(kind, sa, sb, rows)
        it === nothing && continue
        (kind == :ann ? rann : rtr)[gw_pair_key(sa, sb)] = it
        fil[(kind, gw_pair_key(sa, sb)...)] = (αg, ag, how)
        push!(used, f)
        # per-ℓ̃ power fractions (angular pattern of the line); a row without an ℓ̃ has 0 there
        fr = [gw_kerr_parse_lfrac(r.lfrac) for r in ok]
        ls = sort(unique(vcat([collect(keys(d)) for d in fr]...)))
        if !isempty(ls)
            dl = Dict{Int,GWInterp2D}()
            for l in ls
                itl = gw_interp_from_points([(r.alpha, r.a, get(d, l, 0.0)) for (r, d) in zip(ok, fr)])
                itl === nothing || (dl[l] = itl)
            end
            isempty(dl) || (lfr[(kind, gw_pair_key(sa, sb)...)] = dl)
        end
        # horizon fraction h interpolated as asinh(h/GW_HFRAC_SCALE) (changes sign, spans decades)
        hp = [(r.alpha, r.a, asinh(r.Gamma_H / r.Gamma_kerr / GW_HFRAC_SCALE)) for r in ok if isfinite(r.Gamma_H) && r.Gamma_kerr != 0]
        if !isempty(hp)
            ith = gw_interp_from_points(hp)
            ith === nothing || (hfr[(kind, gw_pair_key(sa, sb)...)] = ith)
        end
        for r in ok
            for (lv, w) in ((sa, r.Re_omega_a), (sb, r.Re_omega_b))
                isfinite(w) && w > 0 || continue
                get!(get!(ωpts, lv, Dict{Tuple{Float64,Float64}, Float64}()), (round(r.alpha, sigdigits=10), round(r.a, sigdigits=10)),
                     w - gw_omega(lv..., r.alpha, r.a))
            end
        end
    end
    isempty(used) && return nothing
    ω = Dict{NTuple{3,Int}, GWInterp2D}()
    for (lv, d) in ωpts
        it = gw_interp_from_points([(k[1], k[2], v) for (k, v) in d])
        it === nothing || (ω[lv] = it)
    end
    tabs = GWKerrTables(rann, rtr, ω, hfr, lfr, used, fil)
    if verbose
        println("gw_kerr: $(length(used)) tables from $(dir): $(length(rann)) annihilation, $(length(rtr)) transition channels, $(length(ω)) levels with Leaver frequencies")
        cnt = Dict{Int8,Int}()
        for (_, _, h) in values(fil), c in h
            cnt[c] = get(cnt, c, 0) + 1
        end
        tot = sum(values(cnt))
        println("  ratio grid values: " * join(["$(GW_FILL_NAMES[c]) $(round(100 * cnt[c] / tot, digits=2)) %" for c in sort(collect(keys(cnt)))], ", "))
    end
    GW_TABLE_CACHE[key] = tabs
    return tabs
end

# ----------------------------------------------------------------------------
# Angular pattern of a line: spin-weighted spheroidal harmonics (spectral method)
# ----------------------------------------------------------------------------

"""Wigner small-d matrix element d^j_{m'm}(β), with cos β = c."""
function wigner_d(j::Int, mp::Int, m::Int, c::Float64)
    (abs(m) > j || abs(mp) > j) && return 0.0
    ch, sh = sqrt((1 + c) / 2), sqrt(max(1 - c, 0.0) / 2)
    lf(n) = loggamma(n + 1.0)
    pre = 0.5 * (lf(j + m) + lf(j - m) + lf(j + mp) + lf(j - mp))
    s = 0.0
    for k in max(0, m - mp):min(j + m, j - mp)
        term = exp(pre - lf(j + m - k) - lf(k) - lf(j - k - mp) - lf(k - m + mp))
        s += (isodd(k - m + mp) ? -1 : 1) * term * ch^(2j - 2k + m - mp) * sh^(2k - m + mp)
    end
    return s
end

"""Spin-weighted spherical harmonic ₛY_lm(θ, φ = 0) = sqrt((2l+1)/4π) d^l_{m,-s}(θ) (real; orthonormal on the sphere)."""
gw_sYlm(s::Int, l::Int, m::Int, cθ::Float64) = sqrt((2l + 1) / (4π)) * wigner_d(l, m, -s, cθ)

const GW_SWSH_CACHE = Dict{Tuple{Int,Int,Int,Float64}, Tuple{Float64,Int,Vector{Float64}}}()

"""
    gw_swsh(s, l, m, c) -> (A, lmin, b)
Spin-weighted spheroidal harmonic (c = aω): S = Σ_k b[k] ₛY_{lmin+k-1, m}, unit norm, A → (l-s)(l+s+1) at c = 0.
"""
function gw_swsh(s::Int, l::Int, m::Int, c::Real)
    cF = Float64(c)
    key = (s, l, m, round(cF, sigdigits=8))
    haskey(GW_SWSH_CACHE, key) && return GW_SWSH_CACHE[key]
    lmin = max(abs(m), abs(s))
    l >= lmin || error("gw_swsh: l = $l below lmin = $lmin")
    nb = (l - lmin + 1) + max(6, ceil(Int, 2abs(cF)) + 4)
    ls = lmin:lmin+nb-1
    x, w = gauss(ls[end] + 3)                  # exact for polynomials of degree <= 2 l_max + 5
    Y = [gw_sYlm(s, lk, m, xi) for xi in x, lk in ls]          # nodes × basis
    H = zeros(nb, nb)
    for i in 1:nb, j in i:nb
        c1 = 2π * sum(w[q] * Y[q, i] * Y[q, j] * x[q] for q in eachindex(x))
        c2 = 2π * sum(w[q] * Y[q, i] * Y[q, j] * x[q]^2 for q in eachindex(x))
        H[i, j] = -(cF^2 * c2 - 2cF * s * c1)
        H[j, i] = H[i, j]
    end
    for i in 1:nb
        H[i, i] += (ls[i] - s) * (ls[i] + s + 1)
    end
    F = eigen(Symmetric(H))
    k = l - lmin + 1
    j = argmax(abs.(F.vectors[k, :]))        # the eigenvector that continues ₛY_lm
    b = F.vectors[:, j]
    b[k] < 0 && (b .*= -1)
    res = (F.values[j], lmin, b)
    GW_SWSH_CACHE[key] = res
    return res
end

"""ₛS_lm(θ; c) at cos θ, from gw_swsh."""
function gw_swsh_eval(s::Int, l::Int, m::Int, c::Real, cθ::Real)
    _, lmin, b = gw_swsh(s, l, m, c)
    x = clamp(Float64(cθ), -1.0, 1.0)
    return sum(b[k] * gw_sYlm(s, lmin + k - 1, m, x) for k in eachindex(b))
end

"""
    gw_rel_pattern(tabs, kind, sa, sb, α, a, cos_iota) -> 4π (dP/dΩ)/P, or nothing
Line pattern Σ_ℓ̃ f_ℓ̃ 2π [S(θ)² + S(π - θ)²], S = ₋₂S_{ℓ̃m̃}(θ; aω̃), ℓ̃ summed incoherently; nothing without a table.
"""
function gw_rel_pattern(tabs::GWKerrTables, kind::Symbol, sa, sb, α, a, cι)
    fr = get(tabs.lfrac, (kind, gw_pair_key(sa, sb)...), nothing)
    fr === nothing && return nothing
    m̃ = kind == :ann ? sa[3] + sb[3] : sa[3] - sb[3]
    ωl(s) = gw_omega(s..., α, a) + (haskey(tabs.omega, s) ? gw_interp2(tabs.omega[s], α, a) : 0.0)
    ω̃ = α * (kind == :ann ? ωl(sa) + ωl(sb) : abs(ωl(sa) - ωl(sb)))     # M = 1 units
    c = a * ω̃
    x = clamp(Float64(cι), -1.0, 1.0)
    num = den = 0.0
    for (l, it) in fr
        f = max(gw_interp2(it, α, a), 0.0)
        f > 0 || continue
        num += f * 2π * (gw_swsh_eval(-2, l, m̃, c, x)^2 + gw_swsh_eval(-2, l, m̃, c, -x)^2)
        den += f
    end
    return den > 0 ? num / den : nothing
end

# ----------------------------------------------------------------------------
# Default channel set for the tables
# ----------------------------------------------------------------------------

"""
    gw_kerr_default_channels(; lmax=5, overtones=1) -> Vector{(kind, sa, sb)}
Ground-pair and overtone-self annihilations, transitions with an NR multipole, among l = m levels with l ≤ lmax.
"""
function gw_kerr_default_channels(; lmax::Int=5, overtones::Int=1)
    ground = [(l + 1, l, l) for l in 1:lmax]
    over = [(l + 1 + k, l, l) for l in 1:lmax for k in 1:overtones]
    levels = vcat(ground, over)
    ch = Tuple{Symbol,NTuple{3,Int},NTuple{3,Int}}[]
    for (i, sa) in enumerate(ground), sb in ground[i:end]
        push!(ch, (:ann, sa, sb))
    end
    for sa in over
        push!(ch, (:ann, sa, sa))
    end
    higher(sa, sb) = sa[1] > sb[1] || (sa[1] == sb[1] && sa[2] > sb[2])
    for sa in levels, sb in levels
        higher(sa, sb) || continue
        isempty(gw_transition_channels(sa, sb)) && continue
        push!(ch, (:tr, sa, sb))
    end
    return ch
end
