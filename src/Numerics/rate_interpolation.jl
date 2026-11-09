"""
    rate_interpolation.jl

SR rate Γ(α, a) of each level from the Leaver tables in rate_sve/ (Imag_ergC_{pos,neg}_<nlm>.npz,
growth and damping rates on an (α, a) grid; floored values ~1e-100 mark "no rate on this side").

Default (SR_RATE_INTERP[] = :factorised): near threshold Γ = m [Ω_H(a) - Ω_H(a_c(α))] A(α, a) with
A > 0 and smooth across a_c, so Γ is evaluated as that product, with a_c(α) from the zero-crossing
table (Imag_zeroC_<nlm>.dat) and ln A interpolated linearly in (ln α, a). Γ then vanishes exactly at
a_c(α), and a_c moves continuously as α = μ M(t) changes.

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

# :factorised (default) or :legacy, see the module docstring
const SR_RATE_INTERP = Ref(:factorised)

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

Zero-crossing table α -> a_c, linear in α (linear extrapolation), without the root finder's
failures: points that are not bound (ω/μ > 1 + ZERO_W_TOL, e.g. 511 at α = 0.231), points at the
a = maxSpin clamp ("no root"), and points where a_c stops increasing with α (the 432 and 831 tails,
where a_c -> 1). Returns nothing if the file is missing or fewer than two points remain.
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
    return LinearInterpolation(al[keep], ac[keep], extrapolation_bc=Line())
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

"""
    pre_computed_sr_rates_unified(n, l, m, alph, M; cheby=true) -> (a_mid, rate(a), spin range)

SR rate of (n, l, m) at α as a function of spin (eV), with the interpolation chosen by SR_RATE_INTERP.
Levels without rate tables get rate 0.
"""
function pre_computed_sr_rates_unified(n, l, m, alph, M; cheby::Bool=true)
    SR_RATE_INTERP[] === :legacy && return sr_rates_linear_legacy(n, l, m, alph, M; cheby=cheby)
    SR_RATE_INTERP[] === :factorised || error("SR_RATE_INTERP must be :factorised or :legacy")
    state_str = format_state_for_filename(n, l, m)
    rate_sve_dir = joinpath(@__DIR__, "..", "rate_sve")
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
