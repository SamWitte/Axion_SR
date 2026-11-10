"""
    rate_interpolation.jl

SR rate Γ(α, a) of each level from the Leaver tables in rate_sve/ (Imag_ergC_{pos,neg}_<nlm>.npz,
growth and damping rates on an (α, a) grid; floored values ~1e-100 mark "no rate on this side").

Default (SR_RATE_INTERP[] = :anchored): each α row of the tables is taken as a function of the distance
from threshold, u = Ω_H(a) - Ω_H(a_c(α_row)), with a_c(α) from the zero-crossing table
(Imag_zeroC_<nlm>.dat). Writing the signed rate as G = u A (A > 0 on both sides), each row holds ln A
at its table points plus a value at u = 0 (the slope of G at threshold, see anchored_rows), joined by a
monotone cubic (PCHIP) and continued smoothly (capped) beyond its outermost points. At (α, a),
u = Ω_H(a) - Ω_H(a_c(α)) and the four rows around α are combined at that same u, ln A as a PCHIP in
ln α. The rate vanishes exactly at a_c(α), is C1 in a and in α, reproduces the table points, and rows
are never mixed across their thresholds. Near threshold at large α, ω_I(a) is lopsided (a steep rise to
a near-plateau above a_c, an orders-of-magnitude drop below it); :factorised could not follow this and
overshot up to 18x just above a_c (544 at α = 1.49).

Tables for 322, 544, l = 5 (651-855) and l = 6 with m < 6 (761-865) recomputed on 2026-11-10 with
double-double Leaver, each mode tracked continuously from small α (the earlier tables jumped to the
n+1 overtone at large α, and l = 5, 6 were the perturbative formula); α <= 2, 49 spins (9 extra
nodes in [0.94, 0.9965]), so the spin grid need not be uniform. Against direct Leaver at α between
table rows, more than 0.01 from a_c: median error 0.01-0.3%, at most ~5-18% at large α (653 at
α = 1.33, from the row spacing); a few 10% at the node closest to a_c. a_c(α) matches Leaver to
< 1e-6 (up to 3e-4 where a_c -> 0.998).

:factorised (default on 2026-11-09 and 10): Γ = m [Ω_H(a) - Ω_H(a_c(α))] A(α, a), ln A linear in
(ln α, a), the two rows blended at fixed a. Same zero crossing, but kinks at the spin nodes, and an
overshoot near a_c when a bracketing row is already past its threshold.

:legacy interpolates the signed rate linearly in α between table nodes (~7% apart) and in a. Since
Γ ∝ α^(4l+5) the larger node dominates and the zero crossing lands up to ~0.02 above a_c, as a
staircase in α (flat between nodes, a jump just above each one); kept for comparisons.
"""

module RateInterpolation

using Interpolations
using NPZ
using DelimitedFiles

# Import constants and utilities
include(joinpath(@__DIR__, "..", "Core/constants.jl"))
include(joinpath(@__DIR__, "..", "state_utils.jl"))

export compute_sr_rates_smooth, pre_computed_sr_rates_unified
export symlog, inverse_symlog

# Floor value used in the data files - anything at or below this is artificial
const RATE_FLOOR = 1e-90

# Rate tables are read once per session (accretion runs rebuild the SR
# interpolators on a grid of BH masses).
const NPZ_CACHE = Dict{String, Any}()
cached_npzread(fn) = get!(() -> npzread(fn), NPZ_CACHE, fn)

# :anchored (default), :factorised or :legacy, see the module docstring
const SR_RATE_INTERP = Ref(:anchored)

# :anchored: table points closer than this (in spin) to a_c are left out (A = G/u there is dominated by
# the uncertainty of a_c itself)
const ANCHOR_GAP = 0.004
# ln A at threshold from the knots either side: :harmonic, :linear or :pchip (see anchored_rows)
const ANCHOR_RULE = Ref(:linear)
# largest drop of ln A from threshold to the first knot on either side (see anchored_rows)
const ANCHOR_MAX_DLNA = 2.0
# (state, cheby) -> AnchoredTables, built once per session
const ANCHORED_CACHE = Dict{Tuple{String, Bool}, Any}()

# Table points closer than this to a_c are left out of ln A (Γ and the prefactor both -> 0 there)
const SPIN_GAP = 0.004
# Zero-crossing points with ω/μ = m Ω_H(a_c)/α above 1 + ZERO_W_TOL are not bound states (outliers of
# the root finder); the tolerance covers the ~1e-4 resolution of a_c at small α
const ZERO_W_TOL = 2e-3

# (state, cheby) -> (a_c interpolation, ln A tables), built once per session
const FACTORISED_CACHE = Dict{Tuple{String, Bool}, Any}()

"""Horizon angular velocity in units of 1/(G M)."""
horizon_omega(a) = a / (2 * (1 + sqrt(1 - a^2)))

"""
    symlog(x, linthresh)

Symmetric logarithm transformation for smooth sign changes.
"""
function symlog(x::Real, linthresh::Real=1e-30)
    return sign(x) * log10(1.0 + abs(x) / linthresh)
end

function inverse_symlog(y::Real, linthresh::Real=1e-30)
    return sign(y) * linthresh * (10.0^abs(y) - 1.0)
end

symlog(x::AbstractArray, linthresh::Real=1e-30) = symlog.(x, linthresh)
inverse_symlog(y::AbstractArray, linthresh::Real=1e-30) = inverse_symlog.(y, linthresh)

"""
    load_zero_crossing(state_str, rate_sve_dir; cheby=true)

Load the zero-crossing data: alpha -> spin where rate = 0
"""
function load_zero_crossing(state_str::String, rate_sve_dir::String; cheby::Bool=true)
    prefix = cheby ? "Imag_zeroC" : "Imag_zero"
    fn = joinpath(rate_sve_dir, "$(prefix)_$(state_str).dat")

    if !isfile(fn)
        return nothing
    end

    zerolist = readdlm(fn)
    itp = LinearInterpolation(zerolist[:, 1], zerolist[:, 2], extrapolation_bc=Line())
    return itp
end

"""
    extract_valid_rates_at_alpha(data_3d, alph, M; floor_threshold=RATE_FLOOR)

Extract rate vs spin at a specific alpha, filtering out floored values.
Returns (spin_array, rate_array, alpha_min) where alpha_min is the smallest
alpha in the dataset. If alph < alpha_min the interpolation is evaluated at
alpha_min; the caller is responsible for applying any power-law rescaling.
"""
function extract_valid_rates_at_alpha(data_3d, alph, M; floor_threshold=RATE_FLOOR)
    alpha_vals = sort(unique(data_3d[:, :, 1]))
    spin_vals  = sort(unique(data_3d[:, :, 2]))
    n_alpha    = length(alpha_vals)
    n_spin     = length(spin_vals)

    rates = reshape(data_3d[:, :, 3], n_alpha, n_spin) .* 2.0 ./ (GNew * M)

    alpha_min = alpha_vals[1]
    alpha_eff = max(alph, alpha_min)

    itp_2d = LinearInterpolation((alpha_vals, spin_vals), rates, extrapolation_bc=Line())
    rates_at_alpha = [itp_2d(alpha_eff, a) for a in spin_vals]

    valid_mask = abs.(rates_at_alpha) .> floor_threshold
    return spin_vals[valid_mask], rates_at_alpha[valid_mask], alpha_min
end

"""
    load_threshold(state_str, m, rate_sve_dir; cheby=true)

Zero-crossing table α -> a_c as a monotone C1 curve (PCHIP; continued linearly with the end slopes),
without the root finder's failures: points that are not bound (ω/μ > 1 + ZERO_W_TOL, e.g. 511 at
α = 0.231), points at the a = maxSpin clamp ("no root"), and points where a_c stops increasing with α
(the 432 and 831 tails, where a_c -> 1). A piecewise-linear a_c(α) puts a slope jump into the rate at
every zeroC point near threshold (kinks along α at fixed spin) and, where a_c bends towards 1,
misplaces the threshold between points. Returns nothing if the file is missing or fewer than two
points remain.
"""
function load_threshold(state_str::String, m::Int, rate_sve_dir::String; cheby::Bool=true)
    prefix = cheby ? "Imag_zeroC" : "Imag_zero"
    fn = joinpath(rate_sve_dir, "$(prefix)_$(state_str).dat")
    isfile(fn) || return nothing
    z = readdlm(fn)
    al, ac = z[:, 1], z[:, 2]
    keep = falses(length(al))
    a_run = -Inf
    for i in eachindex(al)
        w = m * horizon_omega(min(ac[i], 1.0)) / al[i]
        if w < 1 + ZERO_W_TOL && ac[i] < maxSpin - 1e-6 && ac[i] > a_run
            keep[i] = true
            a_run = ac[i]
        end
    end
    count(keep) < 2 && return nothing
    x = Float64.(al[keep]); y = Float64.(ac[keep])
    d = pchip_slopes(x, y)
    # end slopes = the end secants: the extrapolation continues the trend of the last points, as the linear
    # interpolation did (PCHIP's own end slope can be 0, e.g. 543 after the far last point at α = 1.29,
    # which would freeze a_c below 1 and open a spurious growth region at larger α)
    d[1] = (y[2] - y[1]) / (x[2] - x[1]); d[end] = (y[end] - y[end - 1]) / (x[end] - x[end - 1])
    return PchipCurve(x, y, d)
end

"""Monotone C1 curve through (x, y) with PCHIP slopes d; beyond the ends, linear with the end slopes."""
struct PchipCurve
    x::Vector{Float64}
    y::Vector{Float64}
    d::Vector{Float64}
end
function (c::PchipCurve)(t::Real)
    x = c.x; n = length(x)
    t <= x[1] && return c.y[1] + c.d[1] * (t - x[1])
    t >= x[n] && return c.y[n] + c.d[n] * (t - x[n])
    k = clamp(searchsortedlast(x, t), 1, n - 1)
    h = x[k + 1] - x[k]; s = (t - x[k]) / h
    return (1 + 2s) * (1 - s)^2 * c.y[k] + s * (1 - s)^2 * h * c.d[k] + s^2 * (3 - 2s) * c.y[k + 1] +
           s^2 * (s - 1) * h * c.d[k + 1]
end

"""ln A = ln(|Γ| / |m (Ω_H(a) - Ω_H(a_c))|) of one rate table (sgn = +1 growth, -1 damping); each row is
filled along a (see factorised_tables), rows without any usable point are NaN."""
struct LnATable
    sgn::Float64
    lnalpha::Vector{Float64}
    spins::Vector{Float64}
    lnA::Matrix{Float64}
end

"""
    factorised_tables(state_str, l, m, rate_sve_dir; cheby=true) -> (a_c, [LnATable...]) or nothing

Left out of ln A: floored points, points within SPIN_GAP of a_c, and points on the wrong side of a_c
for their table (some rows label the first grid points above threshold as damping, e.g. 433 at
α = 0.121 and α = 0.793, while a_c agrees with the Kerr spectrum there). The gaps are then filled
along a, linearly between usable points and flat beyond them, so that both α nodes around any α have
a value at every spin and ln A is continuous in α.
"""
function factorised_tables(state_str::String, l::Int, m::Int, rate_sve_dir::String; cheby::Bool=true)
    get!(FACTORISED_CACHE, (state_str, cheby)) do
        a_c = load_threshold(state_str, m, rate_sve_dir; cheby=cheby)
        a_c === nothing && return nothing
        prefix = cheby ? "Imag_ergC" : "Imag_erg"
        tabs = LnATable[]
        for (side, sgn) in (("pos", 1.0), ("neg", -1.0))
            fn = joinpath(rate_sve_dir, "$(prefix)_$(side)_$(state_str).npz")
            isfile(fn) || continue
            d = cached_npzread(fn)
            alphas = sort(unique(d[:, :, 1]))
            spins = sort(unique(d[:, :, 2]))
            (d[:, 1, 1] ≈ alphas && d[1, :, 2] ≈ spins) || error("$(fn): expected a sorted (alpha, spin) grid")
            lnA = fill(NaN, length(alphas), length(spins))
            for i in eachindex(alphas)
                ac = min(a_c(alphas[i]), 1.0)
                for k in eachindex(spins)
                    v = d[i, k, 3]
                    (v > RATE_FLOOR && abs(spins[k] - ac) > SPIN_GAP) || continue
                    D = m * (horizon_omega(spins[k]) - horizon_omega(ac))
                    sign(D) == sgn || continue
                    lnA[i, k] = log(v / abs(D))
                end
                ok = findall(!isnan, lnA[i, :])
                isempty(ok) && continue
                lnA[i, :] = length(ok) == 1 ? fill(lnA[i, ok[1]], length(spins)) :
                    LinearInterpolation(spins[ok], lnA[i, ok], extrapolation_bc=Flat()).(spins)
            end
            # rows without any usable point (all floored, or a whole row on the wrong side of a_c, e.g.
            # 433 at α = 0.105, 0.113): linear in ln α between filled rows, α^(4l+5) below the first, flat above
            filled = findall(i -> !isnan(lnA[i, 1]), eachindex(alphas))
            isempty(filled) && continue
            lnal = log.(alphas)
            for i in eachindex(alphas)
                isnan(lnA[i, 1]) || continue
                if i < filled[1]
                    lnA[i, :] = lnA[filled[1], :] .+ (4l + 5) * (lnal[i] - lnal[filled[1]])
                elseif i > filled[end]
                    lnA[i, :] = lnA[filled[end], :]
                else
                    j = searchsortedlast(filled, i)
                    i1, i2 = filled[j], filled[j + 1]
                    t = (lnal[i] - lnal[i1]) / (lnal[i2] - lnal[i1])
                    lnA[i, :] = (1 - t) .* lnA[i1, :] .+ t .* lnA[i2, :]
                end
            end
            push!(tabs, LnATable(sgn, lnal, spins, lnA))
        end
        isempty(tabs) ? nothing : (a_c, tabs)
    end
end

"""
    factorised_rate_at_alpha(a_c, tabs, l, m, alph, M) -> (a_mid, rate(a), spin range) or nothing

ln A at α: each table on its own side of a_c(α) (growth above, damping below), linear in ln α
between its two nodes (only an empty node row leaves the other one alone), then linear in a over
the points of both tables; flat beyond the spin and α range of the tables, and ∝ α^(4l+5) below it
(as the legacy rescaling).
"""
function factorised_rate_at_alpha(a_c, tabs::Vector{LnATable}, l, m, alph, M)
    ac = min(a_c(alph), 1.0)
    acc = Dict{Float64, Tuple{Float64, Int}}()    # spin -> (Σ ln A, count) over the two tables
    lna = log(alph)
    for T in tabs
        x = clamp(lna, T.lnalpha[1], T.lnalpha[end])
        i = clamp(searchsortedlast(T.lnalpha, x), 1, length(T.lnalpha) - 1)
        t = (x - T.lnalpha[i]) / (T.lnalpha[i + 1] - T.lnalpha[i])
        shift = lna < T.lnalpha[1] ? (4l + 5) * (lna - T.lnalpha[1]) : 0.0
        for k in eachindex(T.spins)
            T.sgn * (T.spins[k] - ac) > 0 || continue
            v1, v2 = T.lnA[i, k], T.lnA[i + 1, k]
            v = isnan(v1) ? v2 : isnan(v2) ? v1 : (1 - t) * v1 + t * v2
            isnan(v) && continue
            s = round(T.spins[k], digits=9)
            sv, c = get(acc, s, (0.0, 0))
            acc[s] = (sv + v + shift, c + 1)
        end
    end
    length(acc) < 2 && return nothing
    spins = sort!(collect(keys(acc)))
    itp = LinearInterpolation(spins, [acc[s][1] / acc[s][2] for s in spins], extrapolation_bc=Flat())
    W_c = horizon_omega(ac)
    pref = 2.0 / (GNew * M)
    rate(aspin) = m * (horizon_omega(clamp(aspin, 0.0, 1.0)) - W_c) * exp(itp(aspin)) * pref
    return clamp(ac, minSpin, maxSpin), rate, (spins[1], spins[end])
end

# ----------------------------------------------------------------------------------------------------
# :anchored
# ----------------------------------------------------------------------------------------------------

"""Monotone cubic Hermite slopes at the knots (x, y), as scipy's PchipInterpolator (Fritsch-Carlson):
zero at local extrema, weighted harmonic means of the secants inside, one-sided three-point ends."""
function pchip_slopes(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    n == 1 && return [0.0]
    h = diff(x); δ = diff(y) ./ h
    n == 2 && return [δ[1], δ[1]]
    d = zeros(n)
    for k in 2:(n - 1)
        if δ[k - 1] * δ[k] > 0
            w1 = 2h[k] + h[k - 1]; w2 = h[k] + 2h[k - 1]
            d[k] = (w1 + w2) / (w1 / δ[k - 1] + w2 / δ[k])
        end
    end
    function edge(h0, h1, δ0, δ1)
        dd = ((2h0 + h1) * δ0 - h0 * δ1) / (h0 + h1)
        sign(dd) != sign(δ0) && return 0.0
        (sign(δ0) != sign(δ1) && abs(dd) > 3abs(δ0)) && return 3δ0
        return dd
    end
    d[1] = edge(h[1], h[2], δ[1], δ[2])
    d[n] = edge(h[n - 1], h[n - 2], δ[n - 1], δ[n - 2])
    return d
end

"""One α row at u = Ω_H(a) - Ω_H(a_c(α_row)): knots sorted in u, holding ln A = ln(G/u) with G = G M ω_I
the signed rate (> 0 growth, < 0 damping), so A > 0 on both sides; the anchor u = 0 (index i0) holds
ln A_0, the rate's slope at threshold; d are the PCHIP slopes of ln A."""
struct RateRow
    u::Vector{Float64}
    lnA::Vector{Float64}
    d::Vector{Float64}
    i0::Int
end
has_side(r::RateRow, pos::Bool) = pos ? r.i0 < length(r.u) : r.i0 > 1

# beyond a row's outermost knot ln A continues with the knot's own slope, saturating (tanh) at a factor
# e^EXTRAP_CAP, so the rate is C1 there and cannot run away (top end: a near 0.998 at α between rows)
const EXTRAP_CAP = 0.5

"""ln A(u) of a row: PCHIP inside the knots; beyond the outermost knot, continued with the knot's
slope (C1), capped at ± EXTRAP_CAP."""
function row_lnA(r::RateRow, u::Float64)
    U = r.u; n = length(U)
    n == 1 && return r.lnA[1]
    u >= U[n] && return r.lnA[n] + EXTRAP_CAP * tanh(r.d[n] * (u - U[n]) / EXTRAP_CAP)
    u <= U[1] && return r.lnA[1] + EXTRAP_CAP * tanh(r.d[1] * (u - U[1]) / EXTRAP_CAP)
    k = clamp(searchsortedlast(U, u), 1, n - 1)
    h = U[k + 1] - U[k]; s = (u - U[k]) / h
    h00 = (1 + 2s) * (1 - s)^2; h10 = s * (1 - s)^2; h01 = s^2 * (3 - 2s); h11 = s^2 * (s - 1)
    return h00 * r.lnA[k] + h10 * h * r.d[k] + h01 * r.lnA[k + 1] + h11 * h * r.d[k + 1]
end

"""ln of the rate's slope at the row's threshold (A at u = 0)."""
anchor_lnA(r::RateRow) = r.lnA[r.i0]

"""|u| of the row's first knot on side `pos` (Inf if none): the scale over which the data constrain
the rate near threshold."""
first_knot(r::RateRow, pos::Bool) = pos ? (r.i0 < length(r.u) ? r.u[r.i0 + 1] : Inf) :
                                          (r.i0 > 1 ? -r.u[r.i0 - 1] : Inf)

"""Rows used on one side of threshold (growth or damping): ln α of the rows, the rows, and for each row
the index of the nearest row with points on this side (0 if none), used when a row has none."""
struct RowSet
    lnalpha::Vector{Float64}
    rows::Vector{RateRow}
    src::Vector{Int}
end

struct AnchoredTables
    a_c::Any
    pos::RowSet
    neg::RowSet
    spin_range::Tuple{Float64, Float64}
end

"""Rows at each α node of the tables in `sides` ((pos::Bool, data) pairs on one common grid). The anchor's
ln A_0 is the weighted harmonic mean of A at the first knots on either side (the PCHIP slope of G
there: the shallower side dominates), or A at the one first knot of a row with points on one side."""
function anchored_rows(a_c, alphas, spins, sides)
    lnal = Float64[]; rows = RateRow[]
    for i in eachindex(alphas)
        ac = min(a_c(alphas[i]), 1.0)
        Wc = horizon_omega(ac)
        us = Float64[]; lA = Float64[]
        for (pos, d) in sides, k in eachindex(spins)
            v = d[i, k, 3]
            (v > RATE_FLOOR && (pos ? spins[k] - ac : ac - spins[k]) > ANCHOR_GAP) || continue
            u = horizon_omega(spins[k]) - Wc
            push!(us, u); push!(lA, log(v / (abs(u))))
        end
        isempty(us) && continue
        o = sortperm(us); us = us[o]; lA = lA[o]
        k1 = searchsortedfirst(us, 0.0)              # first knot with u > 0 (us has no zeros)
        hl = k1 > 1 ? -us[k1 - 1] : NaN; hr = k1 <= length(us) ? us[k1] : NaN
        lnA0 = if k1 > 1 && k1 <= length(us)
            if ANCHOR_RULE[] === :harmonic
                Al, Ar = exp(lA[k1 - 1]), exp(lA[k1]); w1 = 2hr + hl; w2 = hr + 2hl
                log((w1 + w2) / (w1 / Al + w2 / Ar))
            elseif ANCHOR_RULE[] === :linear
                (hr * lA[k1 - 1] + hl * lA[k1]) / (hl + hr)
            else   # :pchip: the PCHIP of ln A through the knots without the anchor, at u = 0
                dd = pchip_slopes(us, lA); h = hl + hr; t = hl / h
                (1 + 2t) * (1 - t)^2 * lA[k1 - 1] + t * (1 - t)^2 * h * dd[k1 - 1] + t^2 * (3 - 2t) * lA[k1] + t^2 * (t - 1) * h * dd[k1]
            end
        else
            k1 > 1 ? lA[k1 - 1] : lA[k1]
        end
        # with ln A falling by Δ from threshold to the first knot, |G| = |u| A peaks between them at
        # e^(Δ-1)/Δ times the knot's value: Δ <= ANCHOR_MAX_DLNA keeps that <= 1.36 (a very lopsided row,
        # e.g. 322 at α = 0.95, would otherwise put a 10x bump next to threshold)
        k1 <= length(us) && (lnA0 = min(lnA0, lA[k1] + ANCHOR_MAX_DLNA))
        k1 > 1 && (lnA0 = min(lnA0, lA[k1 - 1] + ANCHOR_MAX_DLNA))
        insert!(us, k1, 0.0); insert!(lA, k1, lnA0)
        push!(lnal, log(alphas[i]))
        push!(rows, RateRow(us, lA, pchip_slopes(us, lA), k1))
    end
    return lnal, rows
end

function row_set(lnal, rows, pos::Bool)
    c = [j for j in eachindex(rows) if has_side(rows[j], pos)]
    src = [isempty(c) ? 0 : c[argmin(abs.(c .- i))] for i in eachindex(rows)]
    return RowSet(lnal, rows, src)
end

"""
    anchored_tables(state_str, m, rate_sve_dir; cheby=true) -> AnchoredTables or nothing

Rows from the pos (growth) and neg (damping) tables at each α node, keeping only points above the
floor, on the right side of a_c(α), and further than ANCHOR_GAP from it (some rows label the first
points above threshold as damping, e.g. 433 at α = 0.121 and 0.793). When the two tables share a grid
each row holds both sides (so the slope at the anchor sees both); otherwise (421, 431, 521, 621, 631,
632, 642, 841, 842, 844: growth on 40 x 30, damping on 60 x 40) each side has its own rows.
"""
function anchored_tables(state_str::String, m::Int, rate_sve_dir::String; cheby::Bool=true)
    get!(ANCHORED_CACHE, (state_str, cheby)) do
        a_c = load_threshold(state_str, m, rate_sve_dir; cheby=cheby)
        a_c === nothing && return nothing
        prefix = cheby ? "Imag_ergC" : "Imag_erg"
        tabs = Dict{Bool, Any}()
        for (side, pos) in (("pos", true), ("neg", false))
            fn = joinpath(rate_sve_dir, "$(prefix)_$(side)_$(state_str).npz")
            isfile(fn) || continue
            d = cached_npzread(fn)
            al = sort(unique(d[:, :, 1])); sp = sort(unique(d[:, :, 2]))
            (d[:, 1, 1] ≈ al && d[1, :, 2] ≈ sp) || error("$(fn): expected a sorted (alpha, spin) grid")
            tabs[pos] = (al, sp, d)
        end
        isempty(tabs) && return nothing
        empty = RowSet(Float64[], RateRow[], Int[])
        if length(tabs) == 2 && size(tabs[true][3]) == size(tabs[false][3]) &&
           tabs[true][1] ≈ tabs[false][1] && tabs[true][2] ≈ tabs[false][2]
            al, sp = tabs[true][1], tabs[true][2]
            lnal, rows = anchored_rows(a_c, al, sp, ((true, tabs[true][3]), (false, tabs[false][3])))
            isempty(rows) && return nothing
            P, N = row_set(lnal, rows, true), row_set(lnal, rows, false)
        else
            sets = Dict{Bool, RowSet}()
            for (pos, (al, sp, d)) in tabs
                lnal, rows = anchored_rows(a_c, al, sp, ((pos, d),))
                sets[pos] = isempty(rows) ? empty : row_set(lnal, rows, pos)
            end
            P, N = get(sets, true, empty), get(sets, false, empty)
            isempty(P.rows) && isempty(N.rows) && return nothing
        end
        sps = vcat([t[2] for t in values(tabs)]...)
        AnchoredTables(a_c, P, N, (minimum(sps), maximum(sps)))
    end
end

"""PCHIP (Fritsch-Carlson) slope at an interior point between secants δ0 (width h0) and δ1 (width h1)."""
function pchip_mid(h0, h1, δ0, δ1)
    δ0 * δ1 > 0 || return 0.0
    w1 = 2h1 + h0; w2 = h1 + 2h0
    return (w1 + w2) / (w1 / δ0 + w2 / δ1)
end
"""PCHIP slope at an end point next to secant δ0 (width h0), the next one being δ1 (width h1)."""
function pchip_end(h0, h1, δ0, δ1)
    dd = ((2h0 + h1) * δ0 - h0 * δ1) / (h0 + h1)
    sign(dd) != sign(δ0) && return 0.0
    (sign(δ0) != sign(δ1) && abs(dd) > 3abs(δ0)) && return 3δ0
    return dd
end

"""ln α stencil of a RowSet at ln α: rows (as src indices, 0 = absent) at nodes i-1, i, i+1, i+2 around
the interval [L_i, L_{i+1}] holding x = clamp(ln α), their ln α, s in [0, 1], and the power-law shift
(4l+5) Δln α below the first row."""
struct Stencil
    k::NTuple{4, Int}
    x::NTuple{4, Float64}
    s::Float64
    shift::Float64
end

function stencil(S::RowSet, la, l)
    L = S.lnalpha; nr = length(L)
    nr == 0 && return Stencil((0, 0, 0, 0), (0.0, 0.0, 0.0, 0.0), 0.0, 0.0)
    shift = la < L[1] ? (4l + 5) * (la - L[1]) : 0.0
    nr == 1 && return Stencil((0, S.src[1], 0, 0), (0.0, L[1], 0.0, 0.0), 0.0, shift)
    x = clamp(la, L[1], L[end])
    i = clamp(searchsortedlast(L, x), 1, nr - 1)
    idx = (i - 1, i, i + 1, i + 2)
    k = ntuple(q -> 1 <= idx[q] <= nr ? S.src[idx[q]] : 0, 4)
    xs = ntuple(q -> 1 <= idx[q] <= nr ? L[idx[q]] : 0.0, 4)
    return Stencil(k, xs, (x - L[i]) / (L[i + 1] - L[i]), shift)
end

"""Cubic Hermite in ln α through the stencil values v (v[1], v[4] may be absent, then NaN), PCHIP
slopes; linear if only the two middle values exist, the one value if only one does."""
function stencil_eval(st::Stencil, v::NTuple{4, Float64})
    vi, vj = v[2], v[3]
    isfinite(vj) || return vi
    isfinite(vi) || return vj
    st.s == 0.0 && return vi
    H = st.x[3] - st.x[2]; δ = (vj - vi) / H
    has_m = st.k[1] != 0 && isfinite(v[1]); has_p = st.k[4] != 0 && isfinite(v[4])
    δm = has_m ? (vi - v[1]) / (st.x[2] - st.x[1]) : 0.0
    δp = has_p ? (v[4] - vj) / (st.x[4] - st.x[3]) : 0.0
    di = has_m ? pchip_mid(st.x[2] - st.x[1], H, δm, δ) : has_p ? pchip_end(H, st.x[4] - st.x[3], δ, δp) : δ
    dj = has_p ? pchip_mid(H, st.x[4] - st.x[3], δ, δp) : has_m ? pchip_end(H, st.x[2] - st.x[1], δ, δm) : δ
    t = st.s
    return (1 + 2t) * (1 - t)^2 * vi + t * (1 - t)^2 * H * di + t^2 * (3 - 2t) * vj + t^2 * (t - 1) * H * dj
end

"""ln|G| of one side at u: ln|u| + the stencil rows' ln A at that u, cubic in ln α."""
function side_lnG(S::RowSet, st::Stencil, u)
    v = ntuple(q -> st.k[q] == 0 ? NaN : row_lnA(S.rows[st.k[q]], u), 4)
    return log(abs(u)) + stencil_eval(st, v) + st.shift
end

"""
    anchored_rate_at_alpha(T, l, m, alph, M) -> (a_mid, rate(a), spin range)

Rate at α: u = Ω_H(a) - Ω_H(a_c(α)); on the side of u, the rows' G at that u (a row without points on
that side uses the nearest row that has some), ln|G| cubic (PCHIP) in ln α through the four rows around α.
Below the first row the rate scales as α^(4l+5); beyond the last row the last row is used. A side
without any table rows continues the other side's rows linearly through the anchor.

The two sides can come from different rows (a row without points on one side, or the separate growth
and damping grids of some levels); their slopes at the anchor then differ, while ω_I(a) has one slope
there. Each side's ln A is then moved to the common A_0 (harmonic mean of the two, so the shallower
side dominates, as at a PCHIP knot) with weight (1 - (u/u_s)^2)^2 for |u| < u_s, u_s the distance of
that side's first knot (linear in ln α between the two middle rows): the rate is C1 through a_c, the
weight and its slope vanish at u_s, and table points are reproduced. With the same rows on both sides
the two slopes are equal and nothing changes.
"""
function anchored_rate_at_alpha(T::AnchoredTables, l, m, alph, M)
    ac = min(T.a_c(alph), 1.0)
    Wc = horizon_omega(ac)
    la = log(alph)
    SP = isempty(T.pos.rows) ? T.neg : T.pos
    SN = isempty(T.neg.rows) ? T.pos : T.neg
    stP, stN = stencil(SP, la, l), stencil(SN, la, l)
    lnA0(S, st) = stencil_eval(st, ntuple(q -> st.k[q] == 0 ? NaN : anchor_lnA(S.rows[st.k[q]]), 4)) + st.shift
    function us(S, st, pos)                     # first-knot distance, linear in ln α between the middle rows
        f(q) = st.k[q] == 0 ? Inf : first_knot(S.rows[st.k[q]], pos)
        fi, fj = f(2), f(3)
        isfinite(fi) && isfinite(fj) && return (1 - st.s) * fi + st.s * fj
        return isfinite(fi) ? fi : fj
    end
    aP, aN = lnA0(SP, stP), lnA0(SN, stN)
    lnA0c = log(2.0) + aP + aN - log(exp(aP) + exp(aN))            # harmonic mean of the two
    uP, uN = us(SP, stP, true), us(SN, stN, false)
    # no knot on a side (its table is missing): that side continues the other's rows, the slopes agree
    cP = isfinite(uP) && isfinite(aP) && isfinite(aN) ? lnA0c - aP : 0.0
    cN = isfinite(uN) && isfinite(aP) && isfinite(aN) ? lnA0c - aN : 0.0
    wt(u, us) = abs(u) < us ? (1 - (u / us)^2)^2 : 0.0
    pref = 2.0 / (GNew * M)
    function rate(aspin)
        u = horizon_omega(clamp(aspin, 0.0, 1.0)) - Wc
        u == 0.0 && return 0.0
        if u > 0
            (stP.k[2] == 0 && stP.k[3] == 0) && return 0.0
            lg = side_lnG(SP, stP, u) + (cP == 0.0 ? 0.0 : cP * wt(u, uP))
            return exp(lg) * pref
        else
            (stN.k[2] == 0 && stN.k[3] == 0) && return 0.0
            lg = side_lnG(SN, stN, u) + (cN == 0.0 ? 0.0 : cN * wt(u, uN))
            return -exp(lg) * pref
        end
    end
    return clamp(ac, minSpin, maxSpin), rate, T.spin_range
end

"""
    pre_computed_sr_rates_unified(n, l, m, alph, M; cheby=true) -> (a_mid, rate(a), spin range)

SR rate of (n, l, m) at α as a function of spin (eV), with the interpolation chosen by SR_RATE_INTERP.
Levels without rate tables get rate 0.
"""
function pre_computed_sr_rates_unified(n, l, m, alph, M; cheby::Bool=true)
    SR_RATE_INTERP[] === :legacy && return sr_rates_linear_legacy(n, l, m, alph, M; cheby=cheby)
    state_str = format_state_for_filename(n, l, m)
    rate_sve_dir = joinpath(@__DIR__, "..", "rate_sve")
    if SR_RATE_INTERP[] === :anchored
        T = anchored_tables(state_str, m, rate_sve_dir; cheby=cheby)
        # no zero-crossing file (or no usable points): the legacy interpolation, which gives 0 without tables
        T === nothing && return sr_rates_linear_legacy(n, l, m, alph, M; cheby=cheby)
        return anchored_rate_at_alpha(T, l, m, alph, M)
    end
    SR_RATE_INTERP[] === :factorised || error("SR_RATE_INTERP must be :anchored, :factorised or :legacy")
    ft = factorised_tables(state_str, l, m, rate_sve_dir; cheby=cheby)
    # no zero-crossing file (or no usable points): the legacy interpolation, which gives 0 without tables
    ft === nothing && return sr_rates_linear_legacy(n, l, m, alph, M; cheby=cheby)
    res = factorised_rate_at_alpha(ft[1], ft[2], l, m, alph, M)
    res === nothing && return clamp(ft[1](alph), minSpin, maxSpin), (aspin -> 0.0), (minSpin, maxSpin)
    return res
end

"""
    sr_rates_linear_legacy(n, l, m, alph, M; cheby=true)

The interpolation used before 2026-11-09 (SR_RATE_INTERP[] = :legacy): valid (non-floored) growth and
damping rates linear in (α, a), averaged where both tables have a point, then linear in a.
"""
function sr_rates_linear_legacy(n, l, m, alph, M; cheby::Bool=true)
    state_str = format_state_for_filename(n, l, m)
    rate_sve_dir = joinpath(@__DIR__, "..", "rate_sve")

    # Get zero-crossing point (for reference, but we won't force it)
    zero_itp = load_zero_crossing(state_str, rate_sve_dir; cheby=cheby)
    if isnothing(zero_itp)
        a_mid = 0.5
    else
        a_mid = zero_itp(alph)
        a_mid = clamp(a_mid, minSpin, maxSpin)
    end

    # Load data files
    prefix = cheby ? "Imag_ergC" : "Imag_erg"
    pos_file = joinpath(rate_sve_dir, "$(prefix)_pos_$(state_str).npz")
    neg_file = joinpath(rate_sve_dir, "$(prefix)_neg_$(state_str).npz")

    # Collect all valid data points
    all_spins = Float64[]
    all_rates = Float64[]

    # Add positive rate data (these are growth rates, keep as positive)
    if isfile(pos_file)
        pos_data = cached_npzread(pos_file)
        spins_pos, rates_pos, alpha_min_pos = extract_valid_rates_at_alpha(pos_data, alph, M)
        if alph < alpha_min_pos
            rates_pos = rates_pos .* (alph / alpha_min_pos)^(4*l + 5)
        end

        for (a, r) in zip(spins_pos, rates_pos)
            if r > 0
                push!(all_spins, a)
                push!(all_rates, r)
            end
        end
    end

    # Add negative rate data (these are damping rates, make negative)
    if isfile(neg_file)
        neg_data = cached_npzread(neg_file)
        spins_neg, rates_neg, alpha_min_neg = extract_valid_rates_at_alpha(neg_data, alph, M)
        if alph < alpha_min_neg
            rates_neg = rates_neg .* (alph / alpha_min_neg)^(4*l + 5)
        end

        for (a, r) in zip(spins_neg, rates_neg)
            if r > 0
                push!(all_spins, a)
                push!(all_rates, -r)
            end
        end
    end

    # Sort by spin
    sort_idx = sortperm(all_spins)
    all_spins = all_spins[sort_idx]
    all_rates = all_rates[sort_idx]

    # Handle edge cases
    if length(all_spins) < 3
        interp_func = aspin -> 0.0
        return a_mid, interp_func, (minSpin, maxSpin)
    end

    # Remove duplicates (average rates at same spin)
    unique_spins = Float64[]
    unique_rates = Float64[]
    i = 1
    while i <= length(all_spins)
        current_spin = all_spins[i]
        sum_rate = all_rates[i]
        count = 1

        # Collect all points at approximately the same spin
        while i + count <= length(all_spins) && abs(all_spins[i + count] - current_spin) < 1e-6
            sum_rate += all_rates[i + count]
            count += 1
        end

        push!(unique_spins, current_spin)
        push!(unique_rates, sum_rate / count)
        i += count
    end

    if length(unique_spins) < 3
        interp_func = aspin -> 0.0
        return a_mid, interp_func, (minSpin, maxSpin)
    end
    
    
    # Create linear interpolation directly on rates (no symlog transformation)
    itp = LinearInterpolation(unique_spins, unique_rates, extrapolation_bc=Line())

    # Create the final interpolation function
    function interp_func(aspin)
        return itp(aspin)
    end

    return a_mid, interp_func, (minimum(unique_spins), maximum(unique_spins))
end

"""
    compute_sr_rates_smooth(qtm_cfigs, M_BH, aBH, alph; cheby=true)

SR rates of every level at (α, aBH) and their spin interpolants (pre_computed_sr_rates_unified).
"""
function compute_sr_rates_smooth(qtm_cfigs, M_BH, aBH, alph; cheby::Bool=true)
    SR_rates = zeros(length(qtm_cfigs))
    interp_functions = []

    for (idx, (n, l, m, alph_threshold)) in enumerate(qtm_cfigs)

        a_mid, interp_func, a_range = pre_computed_sr_rates_unified(n, l, m, alph, M_BH; cheby=cheby)
        push!(interp_functions, interp_func)
        SR_rates[idx] = interp_func(aBH)
    end

    interp_dict = Dict{Tuple{Int,Int,Int}, Function}()
    for (i, (n, l, m, _)) in enumerate(qtm_cfigs)
        interp_dict[(n, l, m)] = interp_functions[i]
    end

    return SR_rates, interp_functions, interp_dict
end

end # module RateInterpolation
